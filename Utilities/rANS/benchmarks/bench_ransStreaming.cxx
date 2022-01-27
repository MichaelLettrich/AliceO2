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
      std::tie(outputIter, state) = ransRenorm<decltype(outputIter), (1ull << LOWER_BOUND_BITS), STREAM_BITS>(state, frequency, static_cast<uint8_t>(SYMBOL_TABLE_PRECISION_BITS), outputIter);
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
