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
#include <cmath>

#include <benchmark/benchmark.h>
#include <boost/align/aligned_allocator.hpp>

#include "rANS/utils.h"
#include "rANS/rans.h"
#include "rANS/internal/backend/simd/kernel.h"

#ifdef ENABLE_VTUNE_PROFILER
#include <ittnotify.h>
#endif

class RANSData
{
 public:
  using source_T = uint16_t;

  RANSData()
  {
    mSourceMessage = [this] {
      std::mt19937 mt(0); // same seed we want always the same distrubution of random numbers;
      std::binomial_distribution<source_T> dist(this->mMax, this->mProbability);

      std::vector<source_T> symbols(this->mMax);
      for (size_t i = 0; i < mNSymbols; ++i) {
        symbols.push_back(dist(mt));
      }
      return symbols;
    }();
  };

  auto buildSymbolStats() const
  {
    using namespace o2::rans;

    // build a symbol table
    FrequencyTable ft{};
    ft.addSamples(std::begin(mSourceMessage), std::end(mSourceMessage), 0, mMax);
    return internal::SymbolStatistics{ft, mRescalingBits};
  };

  inline const auto& getSourceMessage() const noexcept { return mSourceMessage; };

  inline source_T getMin() const noexcept { return mMin; };
  inline source_T getMax() const noexcept { return mMax; };

 protected:
  static constexpr source_T mMin{0};
  static constexpr source_T mMax{std::numeric_limits<source_T>::max()};
  static constexpr size_t mRescalingBits = 24;
  static constexpr size_t mNSymbols = 1ul << 24;
  static constexpr double mProbability = 0.3;
  std::vector<source_T> mSourceMessage{};
};

class RANSDataFixture : public RANSData, public benchmark::Fixture
{
 public:
  using source_T = uint16_t;

  void SetUp(const ::benchmark::State& state) override{};

  void TearDown(const ::benchmark::State& state) override{};
};

BENCHMARK_F(RANSDataFixture, accessSOA)
(benchmark::State& s)
{
  using namespace o2::rans;
  using coder_T = uint64_t;
  const auto symbolStats = buildSymbolStats();
  std::vector<uint32_t> frequencies;
  std::vector<uint32_t> cumulative;

  for (const auto [freq, cumul] : symbolStats) {
    frequencies.push_back(freq);
    cumulative.push_back(cumul);
  }

  std::array<uint64_t, 4> states{0, 0, 0, 0};

  for (auto _ : s) {

#pragma GCC unroll(4)
    for (size_t i = 0; i < mSourceMessage.size(); i += 4) {
#pragma omp simd
      for (size_t j = 0; j < 4; ++j) {
        const auto sourceSymbol = mSourceMessage[i + j];
        states[j] += frequencies[sourceSymbol - mMin] + cumulative[sourceSymbol - mMin];
      }
    }
    benchmark::DoNotOptimize(states.data());
  }

  s.SetItemsProcessed(int64_t(s.iterations()) * mSourceMessage.size());
  s.SetBytesProcessed(int64_t(s.iterations()) * mSourceMessage.size() * sizeof(source_T));
};

BENCHMARK_F(RANSDataFixture, accessAOS)
(benchmark::State& s)
{
  using namespace o2::rans;
  using symbol_T = internal::simd::EncoderSymbol;
  using coder_T = uint64_t;
  const auto symbolStats = buildSymbolStats();
  std::vector<symbol_T> symbolTable;

  for (const auto [freq, cumul] : symbolStats) {
    symbolTable.emplace_back(freq, cumul);
  }

  std::array<uint64_t, 4> states{0, 0, 0, 0};

  for (auto _ : s) {

#pragma GCC unroll(4)
    for (size_t i = 0; i < mSourceMessage.size(); i += 4) {
#pragma omp simd
      for (size_t j = 0; j < 4; ++j) {
        const auto sourceSymbol = mSourceMessage[i + j];
        const auto symbol = symbolTable[sourceSymbol - mMin];
        states[j] += symbol.getFrequency() + symbol.getCumulative();
      }
    }
    benchmark::DoNotOptimize(states.data());
  }

  s.SetItemsProcessed(int64_t(s.iterations()) * mSourceMessage.size());
  s.SetBytesProcessed(int64_t(s.iterations()) * mSourceMessage.size() * sizeof(source_T));
};

BENCHMARK_F(RANSDataFixture, accessLightSymbolTable)
(benchmark::State& s)
{
  using namespace o2::rans;
  using symbol_T = internal::simd::EncoderSymbol;
  using coder_T = uint64_t;
  const internal::SymbolTable<symbol_T> symbolTable{buildSymbolStats()};

  std::array<uint64_t, 4> states{0, 0, 0, 0};

  for (auto _ : s) {

#pragma GCC unroll(4)
    for (size_t i = 0; i < mSourceMessage.size(); i += 4) {
#pragma omp simd
      for (size_t j = 0; j < 4; ++j) {
        const auto sourceSymbol = mSourceMessage[i + j];
        const auto encoderSymbol = symbolTable[sourceSymbol];
        states[j] += encoderSymbol.getFrequency() + encoderSymbol.getCumulative();
      }
    }
    benchmark::DoNotOptimize(states.data());
  }

  s.SetItemsProcessed(int64_t(s.iterations()) * mSourceMessage.size());
  s.SetBytesProcessed(int64_t(s.iterations()) * mSourceMessage.size() * sizeof(source_T));
};

BENCHMARK_F(RANSDataFixture, accessHeavySymbolTable)
(benchmark::State& s)
{
  using namespace o2::rans;
  using symbol_T = internal::cpp::EncoderSymbol<uint64_t>;
  using coder_T = uint64_t;
  const internal::SymbolTable<symbol_T> symbolTable{buildSymbolStats()};

  std::array<uint64_t, 4> states{0, 0, 0, 0};

  for (auto _ : s) {

#pragma GCC unroll(4)
    for (size_t i = 0; i < mSourceMessage.size(); i += 4) {
#pragma omp simd
      for (size_t j = 0; j < 4; ++j) {
        const auto sourceSymbol = mSourceMessage[i + j];
        const auto encoderSymbol = symbolTable[sourceSymbol];
        states[j] += encoderSymbol.freq + encoderSymbol.bias;
      }
    }
    benchmark::DoNotOptimize(states.data());
  }

  s.SetItemsProcessed(int64_t(s.iterations()) * mSourceMessage.size());
  s.SetBytesProcessed(int64_t(s.iterations()) * mSourceMessage.size() * sizeof(source_T));
};

template <o2::rans::internal::simd::SIMDWidth SIMDWidth_V>
void accessSOASIMD(benchmark::State& s)
{
  using namespace o2::rans;
  using namespace o2::rans::internal;

  using epi64_t = simd::simdepi64_t<SIMDWidth_V>;
  using epi32_t = simd::simdepi32_t<SIMDWidth_V>;
  using pd_t = simd::simdpd_t<SIMDWidth_V>;

  using source_T = typename RANSData::source_T;
  using coder_T = uint64_t;

  RANSData data;
  const auto symbolStats = data.buildSymbolStats();
  const auto& sourceMessage = data.getSourceMessage();
  std::vector<uint32_t> frequencies;
  std::vector<uint32_t> cumulative;

  for (const auto [freq, cumul] : symbolStats) {
    frequencies.push_back(freq);
    cumulative.push_back(cumul);
  }

  epi64_t states{0};

  for (auto _ : s) {

    for (size_t i = 0; i < sourceMessage.size(); i += SIMDWidth_V) {

      //convert first
      pd_t freqs;
      pd_t cumuls;
      for (size_t j = 0; j < SIMDWidth_V; ++j) {
        const auto sourceSymbol = sourceMessage[i + j];
        const auto tableIndex = sourceSymbol - data.getMin();
        freqs[j] = frequencies[tableIndex];
        cumuls[j] = cumulative[tableIndex];
      }

      //       const auto [freqs, cumuls] = [&frequencies, &cumulative, &data, i]() {
      //         epi32_t freqs_epi32;
      //         epi32_t cumuls_epi32;
      //         const auto& sourceMessage = data.getSourceMessage();
      // #pragma omp simd
      //         for (size_t j = 0; j < SIMDWidth_V; ++j) {
      //           const auto sourceSymbol = sourceMessage[i + j];
      //           const auto tableIndex = sourceSymbol - data.getMin();
      //           freqs_epi32[j] = frequencies[tableIndex];
      //           cumuls_epi32[j] = cumulative[tableIndex];
      //         }

      //         return std::pair(internal::simd::int32ToDouble(freqs_epi32), internal::simd::int32ToDouble(cumuls_epi32));
      //       }();

#pragma omp simd
      for (size_t j = 0; j < SIMDWidth_V; ++j) {
        states[j] += freqs[j] + cumuls[j];
      }
    }
    benchmark::DoNotOptimize(states.data());
  }

  s.SetItemsProcessed(int64_t(s.iterations()) * sourceMessage.size());
  s.SetBytesProcessed(int64_t(s.iterations()) * sourceMessage.size() * sizeof(source_T));
};

BENCHMARK_TEMPLATE(accessSOASIMD, o2::rans::internal::simd::SIMDWidth::SSE);
BENCHMARK_TEMPLATE(accessSOASIMD, o2::rans::internal::simd::SIMDWidth::AVX);

BENCHMARK_MAIN();
