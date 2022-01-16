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
#include <utility>

#include <benchmark/benchmark.h>

#include "rANS/utils.h"
#include "rANS/rans.h"
#include "rANS/internal/backend/simd/kernel.h"

#ifdef ENABLE_VTUNE_PROFILER
#include <ittnotify.h>
#endif

using namespace o2::rans::internal;
using namespace o2::rans::internal::simd;

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

  auto buildFrequencyTable() const
  {
    using namespace o2::rans;

    // build a symbol table
    FrequencyTable frequencyTable{};
    frequencyTable.addSamples(std::begin(mSourceMessage), std::end(mSourceMessage), 0, mMax);
    return renorm(std::move(frequencyTable), mRescalingBits);
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
  const auto frequencyTable = buildFrequencyTable();
  std::vector<uint32_t> frequencies;
  std::vector<uint32_t> cumulative;

  uint32_t cumulatedFrequency = 0;
  for (auto frequency : frequencyTable) {
    frequencies.push_back(frequency);
    cumulative.push_back(cumulatedFrequency);
    cumulatedFrequency += frequency;
  }

  std::array<uint64_t, 4> states{0, 0, 0, 0};

  for (auto _ : s) {

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
  using symbol_T = internal::simd::Symbol;
  using coder_T = uint64_t;
  const auto frequencyTable = buildFrequencyTable();
  std::vector<symbol_T> symbolTable;

  uint32_t cumulatedFrequency = 0;
  for (auto frequency : frequencyTable) {
    symbolTable.emplace_back(frequency, cumulatedFrequency);
    cumulatedFrequency += frequency;
  }

  std::array<uint64_t, 4> states{0, 0, 0, 0};

  for (auto _ : s) {

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
  using symbol_T = internal::simd::Symbol;
  using coder_T = uint64_t;
  const internal::SymbolTable<symbol_T> symbolTable{buildFrequencyTable()};

  std::array<uint64_t, 4> states{0, 0, 0, 0};

  for (auto _ : s) {

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
  const internal::SymbolTable<symbol_T> symbolTable{buildFrequencyTable()};

  std::array<uint64_t, 4> states{0, 0, 0, 0};

  for (auto _ : s) {

    for (size_t i = 0; i < mSourceMessage.size(); i += 4) {
#pragma omp simd
      for (size_t j = 0; j < 4; ++j) {
        const auto sourceSymbol = mSourceMessage[i + j];
        const auto encoderSymbol = symbolTable[sourceSymbol];
        states[j] += encoderSymbol.getFrequency() + encoderSymbol.getBias();
      }
    }
    benchmark::DoNotOptimize(states.data());
  }

  s.SetItemsProcessed(int64_t(s.iterations()) * mSourceMessage.size());
  s.SetBytesProcessed(int64_t(s.iterations()) * mSourceMessage.size() * sizeof(source_T));
};

inline std::pair<simd::epi32_t<simd::SIMDWidth::SSE>, simd::epi32_t<simd::SIMDWidth::SSE>> aosToSoaSSE(const uint32_t* in0, const uint32_t* in1) noexcept
{

  __m128i in0vec = _mm_load_si128(reinterpret_cast<__m128i const*>(in0));
  __m128i in1vec = _mm_load_si128(reinterpret_cast<__m128i const*>(in1));

  __m128i resvec0 = _mm_unpacklo_epi32(in0vec, in1vec);
  __m128i resvec1 = _mm_shuffle_epi32(resvec0, _MM_SHUFFLE(0, 0, 3, 2));

  simd::epi32_t<simd::SIMDWidth::SSE> outvec0;
  simd::epi32_t<simd::SIMDWidth::SSE> outvec1;
  _mm_store_si128(reinterpret_cast<__m128i*>(outvec0.data()), resvec0);
  _mm_store_si128(reinterpret_cast<__m128i*>(outvec1.data()), resvec1);

  return {outvec0, outvec1};
};

inline simd::epi64_t<simd::SIMDWidth::SSE> add(const simd::epi32_t<simd::SIMDWidth::SSE>& in0, const simd::epi32_t<simd::SIMDWidth::SSE>& in1, const simd::epi64_t<simd::SIMDWidth::SSE>& in2) noexcept
{
  __m128i in0vec = _mm_load_si128(reinterpret_cast<__m128i const*>(in0.data()));
  __m128i in1vec = _mm_load_si128(reinterpret_cast<__m128i const*>(in1.data()));
  __m128i in2vec = _mm_load_si128(reinterpret_cast<__m128i const*>(in2.data()));

  __m128i result = _mm_add_epi32(in0vec, in1vec);
  result = _mm_cvtepi32_epi64(result);
  result = _mm_add_epi64(result, in2vec);

  simd::epi64_t<simd::SIMDWidth::SSE> resultvec;
  _mm_store_si128(reinterpret_cast<__m128i*>(resultvec.data()), result);
  return resultvec;
};

void accessAOSsse(benchmark::State& s)
{
  using namespace o2::rans;
  using namespace o2::rans::internal;

  using epi64_t = simd::epi64_t<simd::SIMDWidth::SSE>;
  using epi32_t = simd::epi32_t<simd::SIMDWidth::SSE>;
  using pd_t = simd::pd_t<simd::SIMDWidth::SSE>;

  using source_T = typename RANSData::source_T;
  using coder_T = uint64_t;

  RANSData data;
  const auto& sourceMessage = data.getSourceMessage();
  const auto frequencyTable = data.buildFrequencyTable();
  std::vector<simd::Symbol> symbolTable;

  uint32_t cumulatedFrequency = 0;
  for (auto frequency : frequencyTable) {
    symbolTable.emplace_back(frequency, cumulatedFrequency);
    cumulatedFrequency += frequency;
  }
  epi64_t states{0};

  for (auto _ : s) {
#ifdef ENABLE_VTUNE_PROFILER
    __itt_resume();
#endif
    for (size_t i = 0; i < sourceMessage.size(); i += elementCount_v<epi64_t>) {

      std::array<const uint32_t*, 2> symbols;
      for (size_t j = 0; j < elementCount_v<epi64_t>; ++j) {
        const auto sourceSymbol = sourceMessage[i + j];
        const auto tableIndex = sourceSymbol - data.getMin();
        symbols[j] = symbolTable[tableIndex].data();
      }

      const auto [freqs, cumuls] = aosToSoaSSE(symbols[0], symbols[1]);
      states = add(freqs, cumuls, states);
    }
    benchmark::DoNotOptimize(states.data());
#ifdef ENABLE_VTUNE_PROFILER
    __itt_pause();
#endif
  }

  s.SetItemsProcessed(int64_t(s.iterations()) * sourceMessage.size());
  s.SetBytesProcessed(int64_t(s.iterations()) * sourceMessage.size() * sizeof(source_T));
};

// inline std::pair<simdepi32_t<AVX>, simdepi32_t<AVX>> aosToSOA4(const Symbol* base, const uint64_t* index) noexcept
// {
//   __m256i indexVec = _mm256_load_si256(reinterpret_cast<__m256i const*>(index));
//   __m256i data = _mm256_i64gather_epi64(reinterpret_cast<const long long int*>(base), indexVec, 8);

//   //simdepi32_t<AVX512>shuffleMask{0u,2u,4u,6u,1u,3u,5u,7u};
//   // __m256i shuffleMaskVec = _mm256_load_si256(reinterpret_cast<__m256i const*>(shuffleMask.data()));

//   data = _mm256_shuffle_epi32(data, _MM_SHUFFLE(3, 1, 2, 0));
//   //data = _mm256_permutevar8x32_epi32(data,shuffleMaskVec);
//   __m128i lo = _mm256_extractf128_si256(data, 0);
//   __m128i hi = _mm256_extractf128_si256(data, 1);

//   __m128i resvec0 = _mm_unpacklo_epi64(lo, hi);
//   __m128i resvec1 = _mm_unpackhi_epi64(lo, hi);

//   simdepi32_t<AVX> outvec0;
//   simdepi32_t<AVX> outvec1;
//   _mm_store_si128(reinterpret_cast<__m128i*>(outvec0.data()), resvec0);
//   _mm_store_si128(reinterpret_cast<__m128i*>(outvec1.data()), resvec1);

//   return {outvec0, outvec1};
// };

inline std::pair<epi32_t<SIMDWidth::AVX>, epi32_t<SIMDWidth::AVX>> aosToSOA4(const Symbol* base, const uint64_t* index) noexcept
{
  __m256i indexVec = _mm256_load_si256(reinterpret_cast<__m256i const*>(index));
  __m256i data = _mm256_i64gather_epi64(reinterpret_cast<const long long int*>(base), indexVec, 8);

  //epi32_t<SIMDWidth::AVX>shuffleMask{0u,2u,4u,6u,1u,3u,5u,7u};
  // __m256i shuffleMaskVec = _mm256_load_si256(reinterpret_cast<__m256i const*>(shuffleMask.data()));

  data = _mm256_shuffle_epi32(data, _MM_SHUFFLE(3, 1, 2, 0));
  //data = _mm256_permutevar8x32_epi32(data,shuffleMaskVec);
  __m128i lo = _mm256_extractf128_si256(data, 0);
  __m128i hi = _mm256_extractf128_si256(data, 1);

  __m128i resvec0 = _mm_unpacklo_epi64(lo, hi);
  __m128i resvec1 = _mm_unpackhi_epi64(lo, hi);

  epi32_t<SIMDWidth::AVX> outvec0;
  epi32_t<SIMDWidth::AVX> outvec1;
  _mm_store_si128(reinterpret_cast<__m128i*>(outvec0.data()), resvec0);
  _mm_store_si128(reinterpret_cast<__m128i*>(outvec1.data()), resvec1);

  return {outvec0, outvec1};
};

inline std::pair<epi32_t<SIMDWidth::SSE>, epi32_t<SIMDWidth::SSE>> aosToSOAavx(const uint32_t* in0, const uint32_t* in1, const uint32_t* in2, const uint32_t* in3) noexcept
{
  __m128i in0vec = _mm_load_si128(reinterpret_cast<__m128i const*>(in0));
  __m128i in1vec = _mm_load_si128(reinterpret_cast<__m128i const*>(in1));
  __m128i in2vec = _mm_load_si128(reinterpret_cast<__m128i const*>(in2));
  __m128i in3vec = _mm_load_si128(reinterpret_cast<__m128i const*>(in3));

  __m128i merged0 = _mm_unpacklo_epi32(in0vec, in1vec);
  __m128i merged1 = _mm_unpacklo_epi32(in2vec, in3vec);
  __m128i resvec0 = _mm_unpacklo_epi64(merged0, merged1);
  __m128i resvec1 = _mm_unpackhi_epi64(merged0, merged1);

  epi32_t<SIMDWidth::SSE> outvec0;
  epi32_t<SIMDWidth::SSE> outvec1;
  _mm_store_si128(reinterpret_cast<__m128i*>(outvec0.data()), resvec0);
  _mm_store_si128(reinterpret_cast<__m128i*>(outvec1.data()), resvec1);

  return {outvec0, outvec1};
};

inline epi64_t<SIMDWidth::AVX> add(const epi32_t<SIMDWidth::SSE>& in0, const epi32_t<SIMDWidth::SSE>& in1, const epi64_t<SIMDWidth::AVX>& in2) noexcept
{
  __m128i in0vec = _mm_load_si128(reinterpret_cast<__m128i const*>(in0.data()));
  __m128i in1vec = _mm_load_si128(reinterpret_cast<__m128i const*>(in1.data()));
  __m256i in2vec = _mm256_load_si256(reinterpret_cast<__m256i const*>(in2.data()));

  __m128i tmp = _mm_add_epi32(in0vec, in1vec);
  __m256i result = _mm256_cvtepi32_epi64(tmp);
  result = _mm256_add_epi64(result, in2vec);

  epi64_t<SIMDWidth::AVX> resultvec;
  _mm256_store_si256(reinterpret_cast<__m256i*>(resultvec.data()), result);
  return resultvec;
};

void accessAOSavx(benchmark::State& s)
{
  using namespace o2::rans;
  using namespace o2::rans::internal;

  using source_T = typename RANSData::source_T;
  using coder_T = uint64_t;

  RANSData data;
  const auto& sourceMessage = data.getSourceMessage();
  const auto frequencyTable = data.buildFrequencyTable();
  std::vector<simd::Symbol> symbolTable;

  uint32_t cumulatedFrequency = 0;
  for (auto frequency : frequencyTable) {
    symbolTable.emplace_back(frequency, cumulatedFrequency);
    cumulatedFrequency += frequency;
  }
  simd::epi64_t<SIMDWidth::AVX> states{0};

  for (auto _ : s) {
#ifdef ENABLE_VTUNE_PROFILER
    __itt_resume();
#endif
    const size_t nParallelElements = simd::elementCount_v<simd::epi64_t<simd::SIMDWidth::AVX>>;
    for (size_t i = 0; i < sourceMessage.size(); i += nParallelElements) {

      std::array<const uint32_t*, 4> symbols;
      for (size_t j = 0; j < nParallelElements; ++j) {
        const auto sourceSymbol = sourceMessage[i + j];
        const auto tableIndex = sourceSymbol - data.getMin();
        symbols[j] = symbolTable[tableIndex].data();
      }

      const auto [freqs, cumuls] = aosToSOAavx(symbols[0], symbols[1], symbols[2], symbols[3]);
      states = add(freqs, cumuls, states);
    }
    benchmark::DoNotOptimize(states.data());
#ifdef ENABLE_VTUNE_PROFILER
    __itt_pause();
#endif
  }

  s.SetItemsProcessed(int64_t(s.iterations()) * sourceMessage.size());
  s.SetBytesProcessed(int64_t(s.iterations()) * sourceMessage.size() * sizeof(source_T));
};

BENCHMARK(accessAOSsse);
BENCHMARK(accessAOSavx);

BENCHMARK_MAIN();
