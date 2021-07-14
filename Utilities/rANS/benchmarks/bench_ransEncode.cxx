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
  }

  template <class EncoderSymbol_T>
  auto buildSymbolTable() const -> o2::rans::internal::SymbolTable<EncoderSymbol_T>
  {
    using namespace o2::rans;

    // build a symbol table
    RenormedFrequencyTable renormedFrequencyTable = renorm(makeFrequencyTableFromSamples(std::begin(mSourceMessage), std::end(mSourceMessage), 0, mMax), mRescalingBits);
    return internal::SymbolTable<EncoderSymbol_T>{renormedFrequencyTable};
  }

  const auto& getSourceMessage() const noexcept
  {
    return mSourceMessage;
  }

 protected:
  static constexpr size_t mNSymbols = 1ul << 24;
  static constexpr double mProbability = 0.3;
  static constexpr source_T mMin{0};
  static constexpr source_T mMax{std::numeric_limits<source_T>::max()};
  std::vector<uint16_t> mSourceMessage{};
  static constexpr size_t mRescalingBits = 24;
};

class RANSDataFixture : public RANSData, public benchmark::Fixture
{
 public:
  void SetUp(const ::benchmark::State& state)
  {
  }

  void TearDown(const ::benchmark::State& state)
  {
  }
};

BENCHMARK_F(RANSDataFixture, encodeSlow)
(benchmark::State& s)
{
  using namespace o2::rans;
  using coder_T = uint64_t;
  using symbol_T = internal::simd::EncoderSymbol;

  const auto symbolTable = buildSymbolTable<symbol_T>();
  const auto symbols = [this, &symbolTable]() {
    std::vector<symbol_T> symbols;
    symbols.reserve(mNSymbols);
    for (auto i : mSourceMessage) {
      symbols.emplace_back(symbolTable[i]);
    }
    return symbols;
  }();

  std::vector<coder_T> states(symbols.size(), 1ul << 31);

  for (auto _ : s) {

#pragma omp simd
#pragma GCC unroll(4)
    for (size_t i = 0; i < symbols.size(); ++i) {

      const symbol_T& symbol = symbols[i];
      const coder_T x = states[i];
      const coder_T div = x / symbol.getFrequency();
      const coder_T mod = x - (div * symbol.getFrequency());

      states[i] = (div << mRescalingBits) + mod * symbol.getCumulative();
    }
  }

  s.SetItemsProcessed(int64_t(s.iterations()) * symbols.size());
  s.SetBytesProcessed(int64_t(s.iterations()) * symbols.size() * sizeof(source_T));
};

BENCHMARK_F(RANSDataFixture, encodeFast)
(benchmark::State& s)
{
  __extension__ typedef unsigned __int128 uint128;

  using namespace o2::rans;
  using coder_T = uint64_t;
  using symbol_T = internal::cpp::EncoderSymbol<coder_T>;

  const auto symbolTable = buildSymbolTable<symbol_T>();
  const auto symbols = [this, &symbolTable]() {
    std::vector<symbol_T> symbols;
    symbols.reserve(mNSymbols);
    for (auto i : mSourceMessage) {
      symbols.emplace_back(symbolTable[i]);
    }
    return symbols;
  }();

  std::vector<coder_T> states(symbols.size(), 1ul << 31);

  for (auto _ : s) {
#pragma omp simd
#pragma GCC unroll(2)
    for (size_t i = 0; i < symbols.size(); ++i) {
      const symbol_T& symbol = symbols[i];
      const coder_T q = static_cast<coder_T>((static_cast<uint128>(states[i]) * symbol.getReciprocalFrequency()) >> 64);
      states[i] = states[i] + symbol.getBias() + (q >> symbol.getReciprocalShift()) * symbol.getFrequencyComplement();
    }
  }

  s.SetItemsProcessed(int64_t(s.iterations()) * symbols.size());
  s.SetBytesProcessed(int64_t(s.iterations()) * symbols.size() * sizeof(source_T));
};

template <o2::rans::internal::simd::SIMDWidth SIMDWidth_V>
void encodeSIMD(benchmark::State& s)
{
  using namespace o2::rans;
  using namespace o2::rans::internal;
  using namespace o2::rans::internal::simd;

  RANSData data;
  const auto symbolTable = data.buildSymbolTable<simd::EncoderSymbol>();

  auto [frequencies, cumulative, states] = [&data, &symbolTable](size_t initialState) {
    const auto& sourceMessage = data.getSourceMessage();
    const size_t nSIMDIterations = sourceMessage.size() / elementCount_v<epi64_t<SIMDWidth_V>>;
    std::vector<epi32_t<SIMDWidth_V>> frequencies;
    std::vector<epi32_t<SIMDWidth_V>> cumulative;
    std::vector<epi64_t<SIMDWidth_V>> states;
    frequencies.reserve(nSIMDIterations);
    cumulative.reserve(nSIMDIterations);
    states.reserve(nSIMDIterations);

    for (size_t i = 0; i < sourceMessage.size(); i += elementCount_v<simd::epi64_t<SIMDWidth_V>>) {
      epi32_t<SIMDWidth_V> frequency(0);
      epi32_t<SIMDWidth_V> cumul(0);

      for (size_t j = 0; j < elementCount_v<epi64_t<SIMDWidth_V>>; ++j) {
        auto& symbol = symbolTable[sourceMessage[i + j]];
        frequency[j] = symbol.getFrequency();
        cumul[j] = symbol.getCumulative();
      }
      frequencies.push_back(frequency);
      cumulative.push_back(cumul);

      epi64_t<SIMDWidth_V> state{initialState};
      states.push_back(state);
    }
    return std::make_tuple(frequencies, cumulative, states);
  }(1ul << 37);

  const size_t messageLength = states.size() * elementCount_v<epi64_t<SIMDWidth_V>>;

  for (auto _ : s) {
#ifdef ENABLE_VTUNE_PROFILER
    __itt_resume();
#endif
    [&frequencies = frequencies, &cumulative = cumulative, &states = states]() {
      for (size_t i = 0; i < frequencies.size(); ++i) {
        const pd_t<SIMDWidth_V> frequency = simd::int32ToDouble(frequencies[i]);
        const pd_t<SIMDWidth_V> cumul = simd::int32ToDouble(cumulative[i]);
        const double normalization = 1 << 24;
        auto newState = simd::ransEncode(states[i], frequency, cumul, normalization);
        states[i] = newState;
      }
#ifdef ENABLE_VTUNE_PROFILER
      __itt_pause();
#endif
    }();
  }

  s.SetItemsProcessed(int64_t(s.iterations()) * messageLength);
  s.SetBytesProcessed(int64_t(s.iterations()) * messageLength * sizeof(RANSData::source_T));
};

BENCHMARK_TEMPLATE(encodeSIMD, o2::rans::internal::simd::SIMDWidth::SSE);
BENCHMARK_TEMPLATE(encodeSIMD, o2::rans::internal::simd::SIMDWidth::AVX);

// BENCHMARK_F(RANSData, encodeSIMD)
// (benchmark::State& state)
// {
//   using namespace o2::rans;
//   using coder_T = uint64_t;
//   using symbol_T = internal::simd::EncoderSymbol;

//   const auto symbolTable = buildSymbolTable<symbol_T>();

//   std::vector < internal::simd::simdepi32_t<internal::simd>

//     const auto symbols = [this, &symbolTable]() {
//     std::vector<symbol_T> symbols;
//     symbols.reserve(mNSymbols);
//     for (auto i : mSourceMessage) {
//       symbols.emplace_back(symbolTable[i]);
//     }
//     return symbols;
//   }();

//   std::vector<coder_T> states(symbols.size(), 1ul << 31);

//   for (auto _ : state) {

//     for (size_t i = 0; i < symbols.size(); ++i) {
//       const symbol_T& symbol = symbols[i];
//       const coder_T q = static_cast<coder_T>((static_cast<uint128>(states[i]) * symbol.rcp_freq) >> 64);
//       states[i] = states[i] + symbol.bias + (q >> symbol.rcp_shift) * symbol.cmpl_freq;
//     }
//   }

//   state.counters["Symbols"] = symbols.size();
//   state.counters["Symbols/sec"] = benchmark::Counter(symbols.size(), benchmark::Counter::kIsRate);
// };

BENCHMARK_MAIN();
