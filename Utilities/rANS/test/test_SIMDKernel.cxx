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

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "rANS/internal/backend/simd/kernel.h"

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

using simduint64_types = boost::mpl::list<epi64_t<SIMDWidth::SSE>
#ifdef __AVX2__
                                          , epi64_t<SIMDWidth::AVX>
#endif /* __AVX2__ */
#ifdef __AVX512F__
                                          , epi64_t<SIMDWidth::AVX512>
#endif /* __AVX512F__ */
                                          >;

using simduint32_types = boost::mpl::list<epi32_t<SIMDWidth::SSE>
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

BOOST_AUTO_TEST_CASE_TEMPLATE(simd_uint64ToDouble, simdInt64_T, simduint64_types)
{
  for (size_t i = 0; i < uint64Data.size(); ++i) {
    const simdInt64_T src{uint64Data[i]};
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

BOOST_AUTO_TEST_CASE_TEMPLATE(simd_int32ToDouble, simdInt32_T, simduint32_types)
{
  for (size_t i = 0; i < uint32Data.size(); ++i) {
    const simdInt32_T src{uint32Data[i]};
    const auto dest = int32ToDouble(src);

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
  std::vector<uint32_t> mod = std::vector<uint32_t>(4);
  std::vector<uint32_t> div = std::vector<uint32_t>(4);

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
    pd_T divPD{0};
    pd_T modPD{0};

    if constexpr (std::size(pd_T{}) == 2) {
      //SSE
      __m128d n = _mm_load_pd(numeratorPD.data());
      __m128d d = _mm_load_pd(denominatorPD.data());
      detail::_mm_moddiv_pd(n, d, divPD.data(), modPD.data());
    } else if constexpr (std::size(pd_T{}) == 4) {
      //AVX
      __m256d n = _mm256_load_pd(numeratorPD.data());
      __m256d d = _mm256_load_pd(denominatorPD.data());
      detail::_mm256_moddiv_pd(n, d, divPD.data(), modPD.data());
    } else if constexpr (std::size(pd_T{}) == 8) {
//AVX512
#ifdef __AVX512F__
      __m512d n = _mm512_load_pd(numeratorPD.data());
      __m512d d = _mm512_load_pd(denominatorPD.data());
      detail::_mm512_moddiv_pd(n, d, divPD.data(), modPD.data());
#endif /* __AVX512F__ */
    }

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
  using epi64_T = epi64_t<getSimdWidth_v<pd_T>>;

  const size_t nTests = mFrequency.size();

  for (size_t i = 0; i < nTests; ++i) {
    const epi64_T state{mState};
    const pd_T frequencyPD{mFrequency[i]};
    const pd_T cumulativePD{mCumulative[i]};
    epi64_T result{0};

    result = ransEncode(state, frequencyPD, cumulativePD, mNormalization);

    epi64_T correctStateVector{mResultState[i]};

    BOOST_CHECK_EQUAL_COLLECTIONS(correctStateVector.begin(), correctStateVector.end(), result.begin(), result.end());
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* __SSE__ */