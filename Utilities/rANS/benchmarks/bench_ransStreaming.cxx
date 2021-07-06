// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   bench_ransCombinedIterator.cxx
/// @author Michael Lettrich
/// @since  2021-05-03
/// @brief

#include <vector>
#include <random>

#include <benchmark/benchmark.h>

#include "rANS/utils.h"
#include "rANS/rans.h"
#include "rANS/internal/backend/simd/kernel.h"

#ifdef ENABLE_VTUNE_PROFILER
#include <ittnotify.h>
#endif

using source_t = uint16_t;
using count_t = uint32_t;
using ransState_t = uint64_t;
using stream_t = uint32_t;

using namespace o2::rans;
using namespace o2::rans::internal;

using epi64_t = simd::simdepi64_t<simd::SSE>;
using epi32_t = simd::simdepi32_t<simd::SSE>;
using pd_t = simd::simdpd_t<simd::SSE>;

template <typename T>
inline constexpr size_t toBits()
{
  return sizeof(T) * 8;
};

inline constexpr size_t LOWER_BOUND_BITS = 31;
inline constexpr size_t SYMBOL_TABLE_PRECISION_BITS = 22;
inline constexpr size_t STREAM_BITS = toBits<stream_t>();

struct Fixture : public benchmark::Fixture {
 public:
  void SetUp(const ::benchmark::State& state) override
  {
    std::mt19937 mt(0); // same seed we want always the same distrubution of random numbers;
    std::binomial_distribution<source_t> dist(tries, probability);

    for (size_t i = 0; i < nRepetitions; ++i) {
      frequencies.push_back(dist(mt));
    }
  };

  void
    TearDown(const ::benchmark::State& state) override{};

  static constexpr double probability = 0.3;
  static constexpr size_t tries = 1 << (SYMBOL_TABLE_PRECISION - toBits<source_t>());
  static constexpr size_t nRepetitions = 1ul << 24;
  static constexpr uint64_t initialState = ((LOWER_BOUND >> SYMBOL_TABLE_PRECISION) << STREAM_BITS) * (tries * probability);
  std::vector<count_t> frequencies;
};

template <typename output_IT>
inline auto renorm(ransState_t state, count_t frequency, output_IT outputIter) noexcept
{
  ransState_t newState = state;
  output_IT newOutputIT = outputIter;

  ransState_t maxState = frequency << (LOWER_BOUND_BITS - SYMBOL_TABLE_PRECISION_BITS + STREAM_BITS);
  if (newState >= maxState) {
    ++newOutputIT;
    *newOutputIT = static_cast<stream_t>(state);
    newState >>= STREAM_BITS;
    assert(newState < maxState);
  }
  return std::make_tuple(newOutputIT, newState);
};

template <typename output_IT>
inline auto renorm(const simd::simdepi64_t<simd::SSE> states&, const simd::simdepi32_t<simd::SSE> frequencies&, output_IT outputIter) noexcept
{
  __m128i stateVec = _mm_load_si128(reinterpret_cast<__m128i const*>(states.data()));
  __m128i frequencyVec = _mm_load_si128(reinterpret_cast<__m128i const*>(freqencies.data()));
  constexpr uint8_t shift = LOWER_BOUND_BITS - SYMBOL_TABLE_PRECISION_BITS + STREAM_BITS;
  _m128i maxStateVec = _mm_slli_epi64(frequencyVec, shift);

  //cmpVec = (maxState >= state)
  __m128i cmpGreater = _mm_cmpgt_epi64(stateVec, maxStateVec);
  __m128i cmpEqual = _mm_cmpeq_epi64(stateVec, maxStateVec);
  __m128i cmpVec = _mm_or_si128(cmpGreater, cmpEqual);

  //store new state
  __m128i newStateVec = _mm_srli_epi64(stateVec, STREAM_BITS);
  simd::simdepi64_t<simd::SSE> newState = states;
  _mm_maskmoveu_si128(newStateVec, cmpVec, newState.data());

  //stream out
  constexpr simd::simdepi64_t<simd::AVX> extractionMask{0xFFFFFFFF, 0x0, 0xFFFFFFFF, 0x0, 0xFFFFFFFF, 0x0, 0xFFFFFFFF, 0x0};
  __m128i extractionMaskVec = _mm_load_si128(reinterpret_cast<__m128i const*>(extractionMask.data()));
  __m128i streamOutMask = _mm_and_si128(cmp, extractionMaskVec);
  const uint32_t id = _mm_movemask_epi8(streamOutMask);

  __m128i streamOutVec = _mm_and_si128(stateVec, streamOutMask);
  streamOutVec = _mm_shuffle_epi32(streamOutVec, lut[id]);

  simd::simdepi32_t<simd::AVX> tmpStream;
  _mm_store_si128(reinterpret_cast<__m128i*>(tmpStream.data()), streamOutVec);

  output_IT newOutputIT = outputIter;
  for (size_t i = 0; i < lut[id]; ++i) {
    ++newOutputIT;
    *newOutputIT = tmpStream[i];
  }

  return std::make_tuple(newOutputIT, newState);
};

BENCHMARK_F(Fixture, renorm)
(benchmark::State& s)
{
  std::vector<stream_t> out(nRepetitions, {});
  auto outputIter = out.begin();

  for (auto _ : s) {
#ifdef ENABLE_VTUNE_PROFILER
    __itt_resume();
#endif
    for (size_t i = 0; i < nRepetitions; ++i) {
      const count_t frequency = frequencies[i];
      ransState_t state = initialState;

      std::tie(outputIter, state) = renorm(state, frequency, outputIter);
      benchmark::DoNotOptimize(state);
    }
    benchmark::DoNotOptimize(out.data());
    benchmark::ClobberMemory();
    outputIter = out.begin();
#ifdef ENABLE_VTUNE_PROFILER
    __itt_pause();
#endif
  }

  s.SetItemsProcessed(int64_t(s.iterations()) * nRepetitions);
  s.SetBytesProcessed(int64_t(s.iterations()) * sizeof(stream_t) * nRepetitions);
};

BENCHMARK_MAIN();
