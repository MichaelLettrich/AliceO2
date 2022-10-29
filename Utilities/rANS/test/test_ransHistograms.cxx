// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   test_ransFrequencyTable.cxx
/// @author Michael Lettrich
/// @since  Aug 1, 2020
/// @brief

#define BOOST_TEST_MODULE Utility test
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/mpl/vector.hpp>
#include <gsl/span>

#include "rANS/histogram.h"
#include "rANS/compat.h"

using namespace o2::rans;

using histogram_t = boost::mpl::vector<
  Histogram<char>,
  Histogram<uint8_t>,
  Histogram<int8_t>,
  Histogram<uint16_t>,
  Histogram<int16_t>,
  Histogram<int32_t>>;

namespace std
{
template <typename key_T, typename value_T>
std::ostream& operator<<(std::ostream& os, const std::pair<key_T, value_T>& pair)
{
  os << fmt::format("{}:{}", static_cast<int64_t>(pair.first), static_cast<int64_t>(pair.second));
  return os;
}
} // namespace std

BOOST_AUTO_TEST_CASE_TEMPLATE(test_emptyTables, histogram_T, histogram_t)
{

  using source_type = typename histogram_T::source_type;
  histogram_T histogram{};

  if constexpr (sizeof(source_type) < 4) {
    const size_t tableSize = 1ul << (sizeof(source_type) * 8);

    BOOST_CHECK_EQUAL(histogram.empty(), true);
    BOOST_CHECK_EQUAL(histogram.size(), tableSize);
    BOOST_CHECK(histogram.begin() != histogram.end());
    BOOST_CHECK(histogram.cbegin() != histogram.cend());
  } else {
    using source_type = typename histogram_T::source_type;
    histogram_T histogram{};

    BOOST_CHECK_EQUAL(histogram.empty(), true);

    BOOST_CHECK_EQUAL(histogram.size(), 0);
    BOOST_CHECK(histogram.begin() == histogram.end());
    BOOST_CHECK(histogram.cbegin() == histogram.cend());
  }
};

BOOST_AUTO_TEST_CASE_TEMPLATE(test_addSamples, histogram_T, histogram_t)
{
  using source_type = typename histogram_T::source_type;

  auto computeTableSize = [](const auto& resultsMap) {
    if constexpr (sizeof(source_type) < 4) {
      return 1ul << (sizeof(source_type) * 8);
    } else {
      const auto [minIter, maxIter] = std::minmax_element(std::begin(resultsMap), std::end(resultsMap), [](const auto& a, const auto& b) { return a.first < b.first; });
      return maxIter->first - minIter->first + std::is_signed_v<source_type>;
    }
  };

  const size_t fixedSizeOffset = std::numeric_limits<source_type>::min();

  std::vector<source_type> samples{
    static_cast<source_type>(-5),
    static_cast<source_type>(-2),
    static_cast<source_type>(1),
    static_cast<source_type>(3),
    static_cast<source_type>(5),
    static_cast<source_type>(8),
    static_cast<source_type>(-5),
    static_cast<source_type>(5),
    static_cast<source_type>(1),
    static_cast<source_type>(1),
    static_cast<source_type>(1),
    static_cast<source_type>(14),
    static_cast<source_type>(8),
    static_cast<source_type>(8),
    static_cast<source_type>(8),
    static_cast<source_type>(8),
    static_cast<source_type>(8),
    static_cast<source_type>(8),
    static_cast<source_type>(8),
  };

  std::unordered_map<source_type, uint32_t> results{{static_cast<source_type>(-5), 2},
                                                    {static_cast<source_type>(-2), 1},
                                                    {static_cast<source_type>(-1), 0},
                                                    {static_cast<source_type>(0), 0},
                                                    {static_cast<source_type>(1), 4},
                                                    {static_cast<source_type>(2), 0},
                                                    {static_cast<source_type>(3), 1},
                                                    {static_cast<source_type>(4), 0},
                                                    {static_cast<source_type>(5), 2},
                                                    {static_cast<source_type>(8), 8},
                                                    {static_cast<source_type>(14), 1}};

  size_t tableSize = computeTableSize(results);

  histogram_T histogram{};
  histogram.addSamples(samples.begin(), samples.end());

  histogram_T histogram2{};
  histogram2.addSamples(samples);

  BOOST_CHECK_EQUAL_COLLECTIONS(histogram.begin(), histogram.end(), histogram2.begin(), histogram2.end());

  for (const auto [symbol, value] : results) {
    BOOST_TEST_MESSAGE(fmt::format("testing symbol {}", static_cast<int64_t>(symbol)));
    BOOST_CHECK_EQUAL(histogram[symbol], value);
  }

  BOOST_CHECK_EQUAL(histogram.empty(), false);

  BOOST_CHECK_EQUAL(histogram.size(), tableSize);
  if constexpr (std::is_signed_v<source_type>) {
    BOOST_CHECK_EQUAL(histogram.getOffset(), sizeof(source_type) < 4 ? fixedSizeOffset : -5);
  } else {
    BOOST_CHECK_EQUAL(histogram.getOffset(), sizeof(source_type) < 4 ? fixedSizeOffset : 0);
  }

  BOOST_CHECK(histogram.begin() != histogram.end());
  BOOST_CHECK(histogram.cbegin() != histogram.cend());

  // lets add more frequencies;
  std::vector<source_type> samples2{
    static_cast<source_type>(-10),
    static_cast<source_type>(0),
    static_cast<source_type>(50),
    static_cast<source_type>(-10),
    static_cast<source_type>(0),
    static_cast<source_type>(50),
    static_cast<source_type>(-10),
    static_cast<source_type>(0),
    static_cast<source_type>(50),
    static_cast<source_type>(-10),
    static_cast<source_type>(0),
    static_cast<source_type>(50),
    static_cast<source_type>(-10),
    static_cast<source_type>(0),
    static_cast<source_type>(50),
    static_cast<source_type>(-10),
    static_cast<source_type>(0),
    static_cast<source_type>(50),
  };

  results[static_cast<source_type>(-10)] = 6;
  results[static_cast<source_type>(0)] = 6;
  results[static_cast<source_type>(50)] = 6;

  tableSize = computeTableSize(results);

  histogram.addSamples(samples2.begin(), samples2.end());

  histogram2.addSamples(samples2);

  BOOST_CHECK_EQUAL_COLLECTIONS(histogram.begin(), histogram.end(), histogram2.begin(), histogram2.end());

  for (const auto [symbol, value] : results) {
    BOOST_TEST_MESSAGE(fmt::format("testing symbol {}", static_cast<int64_t>(symbol)));
    BOOST_CHECK_EQUAL(histogram[symbol], value);
  }

  BOOST_CHECK_EQUAL(histogram.empty(), false);
  BOOST_CHECK_EQUAL(histogram.size(), tableSize);
  if constexpr (std::is_signed_v<source_type>) {
    BOOST_CHECK_EQUAL(histogram.getOffset(), sizeof(source_type) < 4 ? fixedSizeOffset : -10);
  } else {
    BOOST_CHECK_EQUAL(histogram.getOffset(), sizeof(source_type) < 4 ? fixedSizeOffset : 0);
  }

  BOOST_CHECK(histogram.begin() != histogram.end());
  BOOST_CHECK(histogram.cbegin() != histogram.cend());

  BOOST_CHECK(histogram.countNUsedAlphabetSymbols() == 10);
  BOOST_CHECK(histogram.getNumSamples() == samples.size() + samples2.size());
};

BOOST_AUTO_TEST_CASE_TEMPLATE(test_addFrequencies, histogram_T, histogram_t)
{
  using source_type = typename histogram_T::source_type;
  using value_type = typename histogram_T::value_type;
  std::vector<value_type> frequencies{0, 1, 2, 3, 4, 5};

  std::unordered_map<source_type, uint32_t> results{
    {static_cast<source_type>(1), 1},
    {static_cast<source_type>(2), 2},
    {static_cast<source_type>(3), 3},
    {static_cast<source_type>(4), 4},
    {static_cast<source_type>(5), 5},
  };

  const size_t fixedtableSize = 1ul << (sizeof(source_type) * 8);
  const size_t fixedSizeOffset = std::numeric_limits<source_type>::min();

  histogram_T histogram{};
  histogram.addFrequencies(frequencies.begin(), frequencies.end(), 0);

  histogram_T histogram2{};
  histogram2.addFrequencies(gsl::make_span(frequencies), 0);

  BOOST_CHECK_EQUAL_COLLECTIONS(histogram.begin(), histogram.end(), histogram2.begin(), histogram2.end());

  for (const auto [symbol, value] : results) {
    BOOST_CHECK_EQUAL(histogram[symbol], value);
  }

  BOOST_CHECK_EQUAL(histogram.empty(), false);
  BOOST_CHECK_EQUAL(histogram.size(), sizeof(source_type) < 4 ? fixedtableSize : 5);
  BOOST_CHECK_EQUAL(histogram.getOffset(), sizeof(source_type) < 4 ? fixedSizeOffset : 1);
  BOOST_CHECK(histogram.begin() != histogram.end());
  BOOST_CHECK(histogram.cbegin() != histogram.cend());
  BOOST_CHECK_EQUAL(histogram.countNUsedAlphabetSymbols(), 5);
  BOOST_CHECK_EQUAL(histogram.getNumSamples(), 15);

  // lets add more frequencies;
  std::vector<value_type> frequencies2{3, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0};

  if constexpr (std::is_signed_v<source_type>) {
    histogram.addFrequencies(frequencies2.begin(), frequencies2.end(), -1);
    histogram2.addFrequencies(gsl::make_span(frequencies2), -1);

    results[static_cast<source_type>(0 + -1)] += 3;
    results[static_cast<source_type>(2 + -1)] += 4;
    results[static_cast<source_type>(11 + -1)] += 5;

    BOOST_CHECK_EQUAL(histogram.size(), sizeof(source_type) < 4 ? fixedtableSize : 12);
    BOOST_CHECK_EQUAL(histogram.getOffset(), sizeof(source_type) < 4 ? fixedSizeOffset : -1);
    BOOST_CHECK_EQUAL(histogram.countNUsedAlphabetSymbols(), 7);
  } else {
    histogram.addFrequencies(frequencies2.begin(), frequencies2.end(), 3);
    histogram2.addFrequencies(gsl::make_span(frequencies2), 3);

    results[static_cast<source_type>(0 + 3)] += 3;
    results[static_cast<source_type>(2 + 3)] += 4;
    results[static_cast<source_type>(11 + 3)] += 5;

    BOOST_CHECK_EQUAL(histogram.size(), sizeof(source_type) < 4 ? fixedtableSize : 14);
    BOOST_CHECK_EQUAL(histogram.getOffset(), sizeof(source_type) < 4 ? fixedSizeOffset : 1);
    BOOST_CHECK_EQUAL(histogram.countNUsedAlphabetSymbols(), 6);
  }
  BOOST_CHECK_EQUAL(histogram.getNumSamples(), 27);
  BOOST_CHECK_EQUAL_COLLECTIONS(histogram.begin(), histogram.end(), histogram2.begin(), histogram2.end());

  for (const auto [symbol, value] : results) {
    BOOST_CHECK_EQUAL(histogram[symbol], value);
  }

  BOOST_CHECK_EQUAL(histogram.empty(), false);
  BOOST_CHECK(histogram.begin() != histogram.end());
  BOOST_CHECK(histogram.cbegin() != histogram.cend());
};

using countingContainer_t = boost::mpl::vector<
  Histogram<uint8_t>,
  Histogram<uint32_t>>;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_renormIncompressible, histogram_T, countingContainer_t)
{
  std::vector<uint32_t> frequencies{1, 1, 2, 2, 2, 2, 6, 8, 4, 10, 8, 14, 10, 19, 26, 30, 31, 35, 41, 45, 51, 44, 47, 39, 58, 52, 42, 53, 50, 34, 50, 30, 32, 24, 30, 20, 17, 12, 16, 6, 8, 5, 6, 4, 4, 2, 2, 2, 1};
  histogram_T histogram{frequencies.begin(), frequencies.end(), static_cast<uint8_t>(0)};

  const size_t scaleBits = 8;

  auto renormedHistogram = renorm(std::move(histogram), scaleBits, true, 1);

  const std::vector<uint32_t> rescaledFrequencies{1, 2, 1, 3, 2, 3, 3, 5, 6, 7, 8, 9, 10, 11, 13, 11, 12, 10, 14, 13, 10, 13, 12, 8, 12, 7, 8, 6, 7, 5, 4, 3, 4, 2, 2, 1, 2, 1, 1};
  BOOST_CHECK_EQUAL(renormedHistogram.isRenormedTo(scaleBits), true);
  BOOST_CHECK_EQUAL(renormedHistogram.getNumSamples(), 1 << scaleBits);
  BOOST_CHECK_EQUAL(renormedHistogram.getIncompressibleSymbolFrequency(), 4);
  BOOST_CHECK_EQUAL_COLLECTIONS(renormedHistogram.begin() + 6, renormedHistogram.begin() + 6 + rescaledFrequencies.size(), rescaledFrequencies.begin(), rescaledFrequencies.end());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_renormLegacy, histogram_T, countingContainer_t)
{
  std::vector<uint32_t> frequencies{1, 1, 2, 2, 2, 2, 6, 8, 4, 10, 8, 14, 10, 19, 26, 30, 31, 35, 41, 45, 51, 44, 47, 39, 58, 52, 42, 53, 50, 34, 50, 30, 32, 24, 30, 20, 17, 12, 16, 6, 8, 5, 6, 4, 4, 2, 2, 2, 1};
  histogram_T histogram{frequencies.begin(), frequencies.end(), static_cast<uint8_t>(0)};

  const size_t scaleBits = 8;

  auto renormedHistogram = compat::renorm(std::move(histogram), scaleBits);
  const std::vector<uint32_t> rescaledFrequencies{1, 1, 1, 2, 1, 2, 1, 2, 2, 2, 2, 3, 3, 4, 6, 7, 7, 9, 9, 11, 12, 10, 11, 9, 13, 12, 10, 13, 11, 8, 12, 7, 7, 6, 7, 4, 4, 3, 4, 1, 2, 1, 2, 2, 2, 1, 2, 1, 1};
  BOOST_CHECK_EQUAL(renormedHistogram.isRenormedTo(scaleBits), true);
  BOOST_CHECK_EQUAL(renormedHistogram.getNumSamples(), 1 << scaleBits);
  BOOST_CHECK_EQUAL(renormedHistogram.getIncompressibleSymbolFrequency(), 2);
  BOOST_CHECK_EQUAL_COLLECTIONS(renormedHistogram.begin(), renormedHistogram.begin() + rescaledFrequencies.size(), rescaledFrequencies.begin(), rescaledFrequencies.end());
}

BOOST_AUTO_TEST_CASE(test_computeMetrics)
{
  std::vector<uint32_t> frequencies{9, 0, 8, 0, 7, 0, 6, 0, 5, 0, 4, 0, 3, 0, 2, 0, 1};
  Histogram<uint32_t> histogram{frequencies.begin(), frequencies.end(), 0};
  const float eps = 1e-4;

  DatasetMetrics<uint32_t> correctMetrics;
  correctMetrics.min = 0;
  correctMetrics.max = 16;
  correctMetrics.alphabetRangeBits = 4;
  correctMetrics.nUsedAlphabetSymbols = 9;
  correctMetrics.entropy = 2.957295041922758;
  correctMetrics.symbolLengthDistribution[2] = 4;
  correctMetrics.symbolLengthDistribution[3] = 3;
  correctMetrics.symbolLengthDistribution[4] = 1;
  correctMetrics.symbolLengthDistribution[5] = 1;
  correctMetrics.weightedSymbolLengthDistribution[2] = 0.66666666;
  correctMetrics.weightedSymbolLengthDistribution[3] = 0.26666666;
  correctMetrics.weightedSymbolLengthDistribution[4] = 0.04444444;
  correctMetrics.weightedSymbolLengthDistribution[5] = 0.02222222;

  const auto metrics = computeDatasetMetrics(histogram);

  BOOST_CHECK_EQUAL(metrics.min, correctMetrics.min);
  BOOST_CHECK_EQUAL(metrics.max, correctMetrics.max);
  BOOST_CHECK_EQUAL(metrics.alphabetRangeBits, correctMetrics.alphabetRangeBits);
  BOOST_CHECK_EQUAL(metrics.nUsedAlphabetSymbols, correctMetrics.nUsedAlphabetSymbols);
  BOOST_CHECK_CLOSE(metrics.entropy, correctMetrics.entropy, eps);

  uint32_t sumUnweighted = 0;
  float_t sumWeighted = 0;
  for (size_t i = 0; i < 32; ++i) {
    BOOST_TEST_MESSAGE(fmt::format("checking length: {}", i));
    BOOST_CHECK_EQUAL(correctMetrics.symbolLengthDistribution[i], metrics.symbolLengthDistribution[i]);
    BOOST_CHECK_CLOSE(correctMetrics.weightedSymbolLengthDistribution[i], metrics.weightedSymbolLengthDistribution[i], eps);

    sumUnweighted += metrics.symbolLengthDistribution[i];
    sumWeighted += metrics.weightedSymbolLengthDistribution[i];
  }

  BOOST_CHECK_EQUAL(preferPacking(correctMetrics), false);
  BOOST_CHECK_EQUAL(computeRenormingPrecision(correctMetrics), internal::MinRenormThreshold);

  BOOST_CHECK_EQUAL(sumUnweighted, correctMetrics.nUsedAlphabetSymbols);
  BOOST_CHECK_CLOSE(sumWeighted, 1.0f, eps);
}

BOOST_AUTO_TEST_CASE(test_RenormingPrecision)
{
  DatasetMetrics<uint32_t> metrics{};
  BOOST_CHECK_EQUAL(preferPacking(metrics), true);

  metrics.nUsedAlphabetSymbols = 7;
  metrics.weightedSymbolLengthDistribution[31] = 1.0f;
  BOOST_CHECK_EQUAL(computeRenormingPrecision(metrics), internal::MaxRenormThreshold);

  metrics.weightedSymbolLengthDistribution[31] = 0.0f;
  metrics.weightedSymbolLengthDistribution[1] = 0.2f;
  metrics.weightedSymbolLengthDistribution[5] = 0.2f;
  metrics.weightedSymbolLengthDistribution[9] = 0.4f;
  metrics.weightedSymbolLengthDistribution[12] = 0.1f;
  metrics.weightedSymbolLengthDistribution[15] = 0.1f;
  BOOST_CHECK_EQUAL(computeRenormingPrecision(metrics), 17);
}

BOOST_AUTO_TEST_CASE(test_PackingDecision)
{
  DatasetMetrics<uint32_t> metrics{};
  BOOST_CHECK_EQUAL(preferPacking(metrics), true);

  metrics.entropy = 1.0f;
  metrics.alphabetRangeBits = 1;
  BOOST_CHECK_EQUAL(preferPacking(metrics), true);

  metrics.entropy = 1.8f;
  metrics.alphabetRangeBits = 2;
  BOOST_CHECK_EQUAL(preferPacking(metrics), false);

  metrics.entropy = 3.14f;
  metrics.alphabetRangeBits = 16;
  BOOST_CHECK_EQUAL(preferPacking(metrics), false);
}