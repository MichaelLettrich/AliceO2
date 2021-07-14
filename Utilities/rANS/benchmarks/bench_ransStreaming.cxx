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
using namespace o2::rans::internal::simd;

template <typename T>
inline constexpr size_t toBits()
{
  return sizeof(T) * 8;
};

inline constexpr size_t LOWER_BOUND_BITS = 31;
inline constexpr size_t SYMBOL_TABLE_PRECISION_BITS = 22;
inline constexpr size_t STREAM_BITS = toBits<stream_t>();
inline constexpr size_t MAX_SHIFT = LOWER_BOUND_BITS - SYMBOL_TABLE_PRECISION_BITS + STREAM_BITS;

alignas(16) inline constexpr std::array<std::array<uint8_t, 16>, 6> permutationLUT{{
  {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}, //0b0000 valid
  {0x00, 0x01, 0x02, 0x03, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}, //0b0001 valid
  {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}, //0b0010 invalid
  {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}, //0b0011 invalid
  {0x08, 0x09, 0x0A, 0x0B, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}, //0b0100 valid
  {0x08, 0x09, 0x0A, 0x0B, 0x00, 0x01, 0x02, 0x03, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}  //0b0101  valid
}};

inline constexpr std::array<int, 6> nStreamElemsLUT{
  0, //0b0000 valid
  1, //0b0001 valid
  0, //0b0010 invalid
  0, //0b0011 invalid
  1, //0b0100 valid
  2  //0b0101  valid
};

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
  static constexpr size_t tries = 1 << (SYMBOL_TABLE_PRECISION_BITS - toBits<source_t>());
  static constexpr size_t nRepetitions = 1ul << 24;
  static constexpr size_t expectationValue = probability * nRepetitions;
  static constexpr ransState_t initialState = expectationValue << MAX_SHIFT;
  std::vector<count_t> frequencies;
};

struct SSEFixture : public Fixture {
 public:
  void SetUp(const ::benchmark::State& state) override
  {
    std::mt19937 mt(0); // same seed we want always the same distrubution of random numbers;
    std::binomial_distribution<source_t> dist(tries, probability);

    for (size_t i = 0; i < nRepetitions; i += 2) {
      frequencies.emplace_back(dist(mt), dist(mt), 0u, 0u);
    }
  };

  void TearDown(const ::benchmark::State& state) override{};

  static constexpr epi64_t<SIMDWidth::SSE> initialState = {expectationValue << MAX_SHIFT, expectationValue << MAX_SHIFT};
  std::vector<epi32_t<SIMDWidth::SSE>> frequencies;
};

template <typename output_IT>
inline auto renorm(ransState_t state, count_t frequency, output_IT outputIter) noexcept
{
  ransState_t newState = state;
  output_IT newOutputIT = outputIter;

  ransState_t maxState = static_cast<ransState_t>(frequency) << MAX_SHIFT;
  if (newState >= maxState) {
    ++newOutputIT;
    *newOutputIT = static_cast<stream_t>(state);
    newState >>= STREAM_BITS;
    assert(newState < maxState);
  }
  return std::make_tuple(newOutputIT, newState);
};

template <typename output_IT>
inline auto renorm(const epi64_t<SIMDWidth::SSE>& states, const epi32_t<SIMDWidth::SSE>& frequencies, output_IT outputIter) noexcept -> std::tuple<output_IT, epi64_t<SIMDWidth::SSE>>
{
  __m128i stateVec = _mm_load_si128(reinterpret_cast<__m128i const*>(states.data()));
  __m128i frequencyVec = _mm_load_si128(reinterpret_cast<__m128i const*>(frequencies.data()));

  // calculate maximum state
  constexpr uint8_t shift = LOWER_BOUND_BITS - SYMBOL_TABLE_PRECISION_BITS + STREAM_BITS;
  frequencyVec = _mm_cvtepi32_epi64(frequencyVec);
  __m128i maxStateVec = _mm_slli_epi64(frequencyVec, shift);

  //cmpVec = (state >= maxState)
  __m128i cmpGreater = _mm_cmpgt_epi64(stateVec, maxStateVec);
  __m128i cmpEqual = _mm_cmpeq_epi64(stateVec, maxStateVec);
  __m128i cmpVec = _mm_or_si128(cmpGreater, cmpEqual);

  //store new state
  __m128i newStateVec = _mm_srli_epi64(stateVec, STREAM_BITS);
  newStateVec = _mm_blendv_epi8(stateVec, newStateVec, cmpVec);
  epi64_t<SIMDWidth::SSE> newState;
  _mm_store_si128(reinterpret_cast<__m128i*>(newState.data()), newStateVec);

  //stream out
  constexpr epi32_t<SIMDWidth::SSE> extractionMask{0xFFFFFFFFu, 0x0u, 0xFFFFFFFFu, 0x0u};
  __m128i extractionMaskVec = _mm_load_si128(reinterpret_cast<__m128i const*>(extractionMask.data()));
  __m128i streamOutMaskVec = _mm_and_si128(cmpVec, extractionMaskVec);
  const uint32_t id = _mm_movemask_ps(_mm_castsi128_ps(streamOutMaskVec));
  __m128i shuffleMaskVec = _mm_load_si128(reinterpret_cast<__m128i const*>(permutationLUT[id].data()));

  __m128i streamOutVec = _mm_and_si128(stateVec, streamOutMaskVec);
  streamOutVec = _mm_shuffle_epi8(streamOutVec, shuffleMaskVec);
  epi32_t<SIMDWidth::SSE> tmpStream;
  _mm_store_si128(reinterpret_cast<__m128i*>(tmpStream.data()), streamOutVec);

  output_IT newOutputIT = outputIter;
  for (size_t i = 0; i < nStreamElemsLUT[id]; ++i) {
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
    // #ifdef ENABLE_VTUNE_PROFILER
    //     __itt_resume();
    // #endif
    for (size_t i = 0; i < nRepetitions; ++i) {
      const count_t frequency = frequencies[i];
      ransState_t state = initialState;

      std::tie(outputIter, state) = renorm(state, frequency, outputIter);
      benchmark::DoNotOptimize(state);
    }
    benchmark::DoNotOptimize(out.data());
    benchmark::ClobberMemory();
    outputIter = out.begin();
    // #ifdef ENABLE_VTUNE_PROFILER
    //     __itt_pause();
    // #endif
  }

  s.SetItemsProcessed(int64_t(s.iterations()) * nRepetitions);
  s.SetBytesProcessed(int64_t(s.iterations()) * sizeof(stream_t) * nRepetitions);
};

BENCHMARK_F(SSEFixture, renormSSE)
(benchmark::State& s)
{
  std::vector<stream_t> out(nRepetitions, {});
  auto outputIter = out.begin();

  for (auto _ : s) {
#ifdef ENABLE_VTUNE_PROFILER
    __itt_resume();
#endif
    for (size_t i = 0; i < nRepetitions / 2; ++i) {
      const epi32_t<SIMDWidth::SSE> frequency = frequencies[i];
      epi64_t<SIMDWidth::SSE> state = initialState;

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
