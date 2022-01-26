// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   test_ransSIMDEncoder.h
/// @author Michael Lettrich
/// @since  2020-04-15
/// @brief  Test rANS SIMD encoder/ decoder

#define BOOST_TEST_MODULE Utility test
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#ifdef __SSE__

#include <vector>
#include <type_traits>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "rANS/internal/backend/simd/types.h"
#include "rANS/internal/backend/simd/kernel.h"
#include "rANS/internal/backend/simd/Symbol.h"

using namespace o2::rans::internal::simd;

// clang-format off
using pd_types = boost::mpl::list<pd_t<SIMDWidth::SSE>
#ifdef __AVX2__
                                      , pd_t<SIMDWidth::AVX>
#endif /* __AVX2__ */
#ifdef __AVX512F__
                                      , pd_t<SIMDWidth::AVX512>
#endif /* __AVX512F__ */
                                      >;

using epi64_types = boost::mpl::list<epi64_t<SIMDWidth::SSE>
#ifdef __AVX2__
                                          , epi64_t<SIMDWidth::AVX>
#endif /* __AVX2__ */
#ifdef __AVX512F__
                                          , epi64_t<SIMDWidth::AVX512>
#endif /* __AVX512F__ */
                                          >;

using epi32_types = boost::mpl::list<epi32_t<SIMDWidth::SSE>
#ifdef __AVX2__
                                          , epi32_t<SIMDWidth::AVX>
#endif /* __AVX2__ */
#ifdef __AVX512F__
                                          , epi32_t<SIMDWidth::AVX512>
#endif /* __AVX512F__ */
                                          >;
// clang-format on

struct ConvertingFixture64 {
  std::vector<uint64_t> uint64Data = {0x0, 0x1, 0xFFFFFFFFFFFFE, 0xFFFFFFFFFFFFF};
  std::vector<double> doubleData;

  ConvertingFixture64()
  {
    for (auto i : uint64Data) {
      doubleData.push_back(static_cast<double>(i));
    }
  };
};

BOOST_FIXTURE_TEST_SUITE(test_SIMDconvert64, ConvertingFixture64)

BOOST_AUTO_TEST_CASE_TEMPLATE(simd_uint64ToDouble, epi64_T, epi64_types)
{
  for (size_t i = 0; i < uint64Data.size(); ++i) {
    const epi64_T src{uint64Data[i]};
    const auto dest = uint64ToDouble(src);

    for (auto elem : dest) {
      BOOST_CHECK_EQUAL(elem, doubleData[i]);
    }
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simd_doubleToUint64, pd_T, pd_types)
{
  for (size_t i = 0; i < doubleData.size(); ++i) {
    const pd_T src{doubleData[i]};
    const auto dest = doubleToUint64(src);

    for (auto elem : dest) {
      BOOST_CHECK_EQUAL(elem, uint64Data[i]);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

struct ConvertingFixture32 {
  std::vector<uint32_t> uint32Data = {0x0, 0x1, 0x7FFFFFFE, 0x7FFFFFFF};
  std::vector<double> doubleData;

  ConvertingFixture32()
  {
    for (auto i : uint32Data) {
      doubleData.push_back(static_cast<double>(i));
    }
  };
};

BOOST_FIXTURE_TEST_SUITE(test_SIMDconvert32, ConvertingFixture32)

BOOST_AUTO_TEST_CASE_TEMPLATE(simd_int32ToDouble, epi32_T, epi32_types)
{
  constexpr SIMDWidth simdWidth_V = simdWidth_v<epi32_T>;
  using simdPD_T = pd_t<simdWidth_V>;

  for (size_t i = 0; i < uint32Data.size(); ++i) {
    const epi32_T src{uint32Data[i]};
    auto dest = int32ToDouble<simdWidth_V>(src);

    for (auto elem : dest) {
      BOOST_CHECK_EQUAL(elem, doubleData[i]);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

struct ModDivFixture {
  std::vector<uint32_t> numerator = {1, 6, 17};
  std::vector<uint32_t> denominator = {1, 3, 4};
  // test 1: mod = 0, div correctly rounded
  // test 2: div = 0, mod correclty rounded
  // test 3: mod, div nonzero and correctly rounded
  std::array<uint32_t, 3> mod;
  std::array<uint32_t, 3> div;

  ModDivFixture()
  {
    for (size_t i = 0; i < numerator.size(); ++i) {
      div[i] = numerator[i] / denominator[i];
      mod[i] = numerator[i] % denominator[i];
    }
  };
};

BOOST_FIXTURE_TEST_SUITE(testModDiv, ModDivFixture)

BOOST_AUTO_TEST_CASE_TEMPLATE(modDiv, pd_T, pd_types)
{
  for (size_t i = 0; i < numerator.size(); ++i) {
    const pd_T numeratorPD{static_cast<double>(numerator[i])};
    const pd_T denominatorPD{static_cast<double>(denominator[i])};

    auto [divPD, modPD] = divMod(numeratorPD, denominatorPD);

    pd_T modResult{static_cast<double>(mod[i])};
    pd_T divResult{static_cast<double>(div[i])};

    BOOST_CHECK_EQUAL_COLLECTIONS(divResult.begin(), divResult.end(), divPD.begin(), divPD.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(modResult.begin(), modResult.end(), modPD.begin(), modPD.end());
  }
}

BOOST_AUTO_TEST_SUITE_END()

struct RANSEncodeFixture {

  uint64_t mState{};
  double mNormalization{};
  std::vector<double> mFrequency{};
  std::vector<double> mCumulative{};
  std::vector<uint64_t> mResultState{};

  RANSEncodeFixture()
  {
    const uint64_t state = 1ul << 21;
    const std::vector<uint32_t> frequency{1, 1, 997, 1234};
    const std::vector<uint32_t> cumulative{0, 321, 1, (1u << 16) - 1234};
    const uint64_t normalization = 1ul << 16;

    //copy and convert to double
    mState = static_cast<double>(state);
    mNormalization = static_cast<double>(normalization);
    std::copy(std::begin(frequency), std::end(frequency), std::back_inserter(mFrequency));
    std::copy(std::begin(cumulative), std::end(cumulative), std::back_inserter(mCumulative));

    // calculate result based on RANS formula
    for (size_t i = 0; i < frequency.size(); ++i) {
      uint64_t resultState = normalization * (state / frequency[i]) + (state % frequency[i]) + cumulative[i];
      mResultState.push_back(resultState);
    }
  };
};

BOOST_FIXTURE_TEST_SUITE(testRANSEncode, RANSEncodeFixture)

BOOST_AUTO_TEST_CASE_TEMPLATE(simd_RansEncode, pd_T, pd_types)
{
  using epi64_T = epi64_t<simdWidth_v<pd_T>>;

  const size_t nTests = mFrequency.size();

  for (size_t i = 0; i < nTests; ++i) {
    const epi64_T state{mState};
    const pd_T frequencyPD{mFrequency[i]};
    const pd_T cumulativePD{mCumulative[i]};
    const pd_T normalizationPD{mNormalization};
    epi64_T result{0};

    result = ransEncode(state, frequencyPD, cumulativePD, normalizationPD);

    epi64_T correctStateVector{mResultState[i]};

    BOOST_CHECK_EQUAL_COLLECTIONS(correctStateVector.begin(), correctStateVector.end(), result.begin(), result.end());
  }
}
BOOST_AUTO_TEST_SUITE_END()

struct AosToSoaFixture {

  std::vector<Symbol> mSource;
  epi32_t<SIMDWidth::AVX512> mFrequencies;
  epi32_t<SIMDWidth::AVX512> mCumulative;

  AosToSoaFixture()
  {
    constexpr size_t nElems = getElementCount<uint32_t>(SIMDWidth::AVX512);
    uint32_t counter = 0;

    for (size_t i = 0; i < nElems; ++i) {
      Symbol symbol{counter++, counter++};
      mFrequencies[i] = symbol.getFrequency();
      mCumulative[i] = symbol.getCumulative();

      mSource.emplace_back(std::move(symbol));
    }
  };
};
using aosToSoa_T = boost::mpl::list<std::integral_constant<size_t, 2>,
                                    std::integral_constant<size_t, 4>,
                                    std::integral_constant<size_t, 8>>;

BOOST_FIXTURE_TEST_SUITE(testAostoSoa, AosToSoaFixture)
BOOST_AUTO_TEST_CASE_TEMPLATE(simd_AosToSOA, sizes_T, aosToSoa_T)
{
  constexpr sizes_T nElements;
  std::array<o2::rans::internal::simd::Symbol, nElements()> aosPtrs;

  for (size_t i = 0; i < nElements(); ++i) {
    aosPtrs[i] = mSource[i];
  }

  auto [frequencies, cumulative] = aosToSoa(aosPtrs);

  for (size_t i = 0; i < nElements(); ++i) {
    BOOST_CHECK_EQUAL(frequencies[i], mFrequencies[i]);
    BOOST_CHECK_EQUAL(cumulative[i], mCumulative[i]);
  };
}
BOOST_AUTO_TEST_SUITE_END()

#endif /* __SSE__ */