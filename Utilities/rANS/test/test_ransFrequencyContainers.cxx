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

#include "rANS/StaticFrequencyTable.h"
#include "rANS/DynamicFrequencyTable.h"
#include "rANS/HashFrequencyTable.h"
#include "rANS/renorm.h"

using namespace o2::rans;

using staticFrequencyTables_t = boost::mpl::vector<
  StaticFrequencyTable<char>,
  StaticFrequencyTable<uint8_t>,
  StaticFrequencyTable<int8_t>,
  StaticFrequencyTable<uint16_t>,
  StaticFrequencyTable<int16_t>>;

using dynamicFrequencyTables_t = boost::mpl::vector<
  DynamicFrequencyTable<char>,
  DynamicFrequencyTable<uint8_t>,
  DynamicFrequencyTable<int8_t>,
  DynamicFrequencyTable<uint16_t>,
  DynamicFrequencyTable<int16_t>,
  DynamicFrequencyTable<int32_t>>;

using hashFrequencyTables_t = boost::mpl::vector<
  HashFrequencyTable<char>,
  HashFrequencyTable<uint8_t>,
  HashFrequencyTable<int8_t>,
  HashFrequencyTable<uint16_t>,
  HashFrequencyTable<int16_t>,
  HashFrequencyTable<uint32_t>,
  HashFrequencyTable<int32_t>>;

namespace std
{
template <typename key_T, typename value_T>
std::ostream& operator<<(std::ostream& os, const std::pair<key_T, value_T>& pair)
{
  os << fmt::format("{}:{}", static_cast<int64_t>(pair.first), static_cast<int64_t>(pair.second));
  return os;
}
} // namespace std

BOOST_AUTO_TEST_CASE_TEMPLATE(test_emptyStaticTables, frequencyTable_T, staticFrequencyTables_t)
{

  using source_type = typename frequencyTable_T::source_type;
  frequencyTable_T frequencyTable{};

  const size_t tableSize = 1ul << (sizeof(source_type) * 8);

  BOOST_CHECK_EQUAL(frequencyTable.empty(), true);
  BOOST_CHECK_EQUAL(frequencyTable.size(), tableSize);
  BOOST_CHECK(frequencyTable.begin() != frequencyTable.end());
  BOOST_CHECK(frequencyTable.cbegin() != frequencyTable.cend());
};

BOOST_AUTO_TEST_CASE_TEMPLATE(test_emptyHashTables, frequencyTable_T, hashFrequencyTables_t)
{

  using source_type = typename frequencyTable_T::source_type;
  frequencyTable_T frequencyTable{};

  BOOST_CHECK_EQUAL(frequencyTable.empty(), true);

  BOOST_CHECK_EQUAL(frequencyTable.size(), 0);
  BOOST_CHECK(frequencyTable.begin() == frequencyTable.end());
  BOOST_CHECK(frequencyTable.cbegin() == frequencyTable.cend());
};

BOOST_AUTO_TEST_CASE_TEMPLATE(test_emptyDynamicTables, frequencyTable_T, dynamicFrequencyTables_t)
{

  using source_type = typename frequencyTable_T::source_type;
  frequencyTable_T frequencyTable{};

  BOOST_CHECK_EQUAL(frequencyTable.empty(), true);

  BOOST_CHECK_EQUAL(frequencyTable.size(), 0);
  BOOST_CHECK(frequencyTable.begin() == frequencyTable.end());
  BOOST_CHECK(frequencyTable.cbegin() == frequencyTable.cend());
};

BOOST_AUTO_TEST_CASE_TEMPLATE(test_addSamplesStatic, frequencyTable_T, staticFrequencyTables_t)
{
  using source_type = typename frequencyTable_T::source_type;

  const size_t tableSize = 1ul << (sizeof(source_type) * 8);

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

  frequencyTable_T frequencyTable{};
  frequencyTable.addSamples(samples.begin(), samples.end());

  frequencyTable_T frequencyTable2{};
  frequencyTable2.addSamples(samples);

  BOOST_CHECK_EQUAL_COLLECTIONS(frequencyTable.begin(), frequencyTable.end(), frequencyTable2.begin(), frequencyTable2.end());

  for (const auto [symbol, value] : results) {
    BOOST_TEST_MESSAGE(fmt::format("testing symbol {}", symbol));
    BOOST_CHECK_EQUAL(frequencyTable[symbol], value);
  }

  BOOST_CHECK_EQUAL(frequencyTable.empty(), false);

  BOOST_CHECK_EQUAL(frequencyTable.size(), tableSize);

  BOOST_CHECK(frequencyTable.begin() != frequencyTable.end());
  BOOST_CHECK(frequencyTable.cbegin() != frequencyTable.cend());

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

  frequencyTable.addSamples(samples2.begin(), samples2.end());

  frequencyTable2.addSamples(samples2);

  BOOST_CHECK_EQUAL_COLLECTIONS(frequencyTable.begin(), frequencyTable.end(), frequencyTable2.begin(), frequencyTable2.end());

  for (const auto [symbol, value] : results) {
    BOOST_CHECK_EQUAL(frequencyTable[symbol], value);
  }

  BOOST_CHECK_EQUAL(frequencyTable.empty(), false);
  BOOST_CHECK_EQUAL(frequencyTable.size(), tableSize);

  BOOST_CHECK(frequencyTable.begin() != frequencyTable.end());
  BOOST_CHECK(frequencyTable.cbegin() != frequencyTable.cend());

  BOOST_CHECK(frequencyTable.computeNUsedAlphabetSymbols() == 10);
  BOOST_CHECK(frequencyTable.getNumSamples() == samples.size() + samples2.size());
};

BOOST_AUTO_TEST_CASE_TEMPLATE(test_addSamplesDynamic, frequencyTable_T, dynamicFrequencyTables_t)
{
  using source_type = typename frequencyTable_T::source_type;

  auto computeTableSize = [](const auto& resultsMap) {
    const auto [minIter, maxIter] = std::minmax_element(std::begin(resultsMap), std::end(resultsMap), [](const auto& a, const auto& b) { return a.first < b.first; });
    return maxIter->first - minIter->first + std::is_signed_v<source_type>;
  };

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

  frequencyTable_T frequencyTable{};
  frequencyTable.addSamples(samples.begin(), samples.end());

  frequencyTable_T frequencyTable2{};
  frequencyTable2.addSamples(samples);

  BOOST_CHECK_EQUAL_COLLECTIONS(frequencyTable.begin(), frequencyTable.end(), frequencyTable2.begin(), frequencyTable2.end());

  for (const auto [symbol, value] : results) {
    BOOST_TEST_MESSAGE(fmt::format("testing symbol {}", static_cast<int64_t>(symbol)));
    BOOST_CHECK_EQUAL(frequencyTable[symbol], value);
  }

  BOOST_CHECK_EQUAL(frequencyTable.empty(), false);

  BOOST_CHECK_EQUAL(frequencyTable.size(), tableSize);
  if constexpr (std::is_signed_v<source_type>) {
    BOOST_CHECK_EQUAL(frequencyTable.getOffset(), -5);
  } else {
    BOOST_CHECK_EQUAL(frequencyTable.getOffset(), 0);
  }

  BOOST_CHECK(frequencyTable.begin() != frequencyTable.end());
  BOOST_CHECK(frequencyTable.cbegin() != frequencyTable.cend());

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

  frequencyTable.addSamples(samples2.begin(), samples2.end());

  frequencyTable2.addSamples(samples2);

  BOOST_CHECK_EQUAL_COLLECTIONS(frequencyTable.begin(), frequencyTable.end(), frequencyTable2.begin(), frequencyTable2.end());

  for (const auto [symbol, value] : results) {
    BOOST_TEST_MESSAGE(fmt::format("testing symbol {}", static_cast<int64_t>(symbol)));
    BOOST_CHECK_EQUAL(frequencyTable[symbol], value);
  }

  BOOST_CHECK_EQUAL(frequencyTable.empty(), false);
  BOOST_CHECK_EQUAL(frequencyTable.size(), tableSize);
  if constexpr (std::is_signed_v<source_type>) {
    BOOST_CHECK_EQUAL(frequencyTable.getOffset(), -10);
  } else {
    BOOST_CHECK_EQUAL(frequencyTable.getOffset(), 0);
  }

  BOOST_CHECK(frequencyTable.begin() != frequencyTable.end());
  BOOST_CHECK(frequencyTable.cbegin() != frequencyTable.cend());

  BOOST_CHECK(frequencyTable.computeNUsedAlphabetSymbols() == 10);
  BOOST_CHECK(frequencyTable.getNumSamples() == samples.size() + samples2.size());
};

BOOST_AUTO_TEST_CASE_TEMPLATE(test_addSamplesHash, frequencyTable_T, hashFrequencyTables_t)
{
  using source_type = typename frequencyTable_T::source_type;

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

  frequencyTable_T frequencyTable{};
  frequencyTable.addSamples(samples.begin(), samples.end());

  frequencyTable_T frequencyTable2{};
  frequencyTable2.addSamples(samples);

  BOOST_CHECK_EQUAL(frequencyTable.size(), frequencyTable2.size());
  for (const auto [symbol, value] : frequencyTable) {
    BOOST_CHECK_EQUAL(frequencyTable2[symbol], value);
  }

  for (const auto [symbol, value] : results) {
    BOOST_TEST_MESSAGE(fmt::format("testing symbol {}", symbol));
    BOOST_CHECK_EQUAL(frequencyTable[symbol], value);
  }

  BOOST_CHECK_EQUAL(frequencyTable.empty(), false);
  BOOST_CHECK_EQUAL(frequencyTable.size(), 7);
  BOOST_CHECK(frequencyTable.begin() != frequencyTable.end());
  BOOST_CHECK(frequencyTable.cbegin() != frequencyTable.cend());

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

  frequencyTable.addSamples(samples2.begin(), samples2.end());
  frequencyTable2.addSamples(samples2);

  BOOST_CHECK_EQUAL(frequencyTable.size(), frequencyTable2.size());
  for (const auto [symbol, value] : frequencyTable) {
    BOOST_CHECK_EQUAL(frequencyTable2[symbol], value);
  }

  for (const auto [symbol, value] : results) {
    BOOST_CHECK_EQUAL(frequencyTable[symbol], value);
  }

  BOOST_CHECK_EQUAL(frequencyTable.empty(), false);
  BOOST_CHECK_EQUAL(frequencyTable.size(), 10);
  BOOST_CHECK(frequencyTable.begin() != frequencyTable.end());
  BOOST_CHECK(frequencyTable.cbegin() != frequencyTable.cend());
  BOOST_CHECK_EQUAL(frequencyTable.computeNUsedAlphabetSymbols(), 10);
  BOOST_CHECK_EQUAL(frequencyTable.getNumSamples(), samples.size() + samples2.size());
};

BOOST_AUTO_TEST_CASE_TEMPLATE(test_addFrequenciesStatic, frequencyTable_T, staticFrequencyTables_t)
{
  using source_type = typename frequencyTable_T::source_type;
  using value_type = typename frequencyTable_T::value_type;
  const size_t tableSize = 1ul << (sizeof(source_type) * 8);
  std::vector<value_type> frequencies{0, 1, 2, 3, 4, 5};

  std::unordered_map<source_type, uint32_t> results{
    {static_cast<source_type>(0), 0},
    {static_cast<source_type>(1), 1},
    {static_cast<source_type>(2), 2},
    {static_cast<source_type>(3), 3},
    {static_cast<source_type>(4), 4},
    {static_cast<source_type>(5), 5},
  };

  frequencyTable_T frequencyTable{};
  frequencyTable.addFrequencies(frequencies.begin(), frequencies.end(), 0);

  frequencyTable_T frequencyTable2{};
  frequencyTable2.addFrequencies(gsl::make_span(frequencies), 0);

  BOOST_CHECK_EQUAL_COLLECTIONS(frequencyTable.begin(), frequencyTable.end(), frequencyTable2.begin(), frequencyTable2.end());

  for (const auto [symbol, value] : results) {
    BOOST_CHECK_EQUAL(frequencyTable[symbol], value);
  }

  BOOST_CHECK_EQUAL(frequencyTable.empty(), false);
  BOOST_CHECK_EQUAL(frequencyTable.size(), tableSize);
  BOOST_CHECK(frequencyTable.begin() != frequencyTable.end());
  BOOST_CHECK(frequencyTable.cbegin() != frequencyTable.cend());
  BOOST_CHECK_EQUAL(frequencyTable.computeNUsedAlphabetSymbols(), 5);
  BOOST_CHECK_EQUAL(frequencyTable.getNumSamples(), 15);

  // lets add more frequencies;
  std::vector<value_type> frequencies2{3, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0};

  if constexpr (std::is_signed_v<source_type>) {
    frequencyTable.addFrequencies(frequencies2.begin(), frequencies2.end(), -1);
    frequencyTable2.addFrequencies(gsl::make_span(frequencies2), -1);

    results[static_cast<source_type>(0 + -1)] += 3;
    results[static_cast<source_type>(2 + -1)] += 4;
    results[static_cast<source_type>(11 + -1)] += 5;

    BOOST_CHECK_EQUAL(frequencyTable.computeNUsedAlphabetSymbols(), 7);
  } else {
    frequencyTable.addFrequencies(frequencies2.begin(), frequencies2.end(), 3);
    frequencyTable2.addFrequencies(gsl::make_span(frequencies2), 3);

    results[static_cast<source_type>(0 + 3)] += 3;
    results[static_cast<source_type>(2 + 3)] += 4;
    results[static_cast<source_type>(11 + 3)] += 5;

    BOOST_CHECK_EQUAL(frequencyTable.computeNUsedAlphabetSymbols(), 6);
  }
  BOOST_CHECK_EQUAL(frequencyTable.getNumSamples(), 27);
  BOOST_CHECK_EQUAL_COLLECTIONS(frequencyTable.begin(), frequencyTable.end(), frequencyTable2.begin(), frequencyTable2.end());

  for (const auto [symbol, value] : results) {
    BOOST_CHECK_EQUAL(frequencyTable[symbol], value);
  }

  BOOST_CHECK_EQUAL(frequencyTable.empty(), false);
  BOOST_CHECK_EQUAL(frequencyTable.size(), tableSize);
  BOOST_CHECK(frequencyTable.begin() != frequencyTable.end());
  BOOST_CHECK(frequencyTable.cbegin() != frequencyTable.cend());
};

BOOST_AUTO_TEST_CASE_TEMPLATE(test_addFrequenciesDynamic, frequencyTable_T, dynamicFrequencyTables_t)
{
  using source_type = typename frequencyTable_T::source_type;
  using value_type = typename frequencyTable_T::value_type;
  std::vector<value_type> frequencies{0, 1, 2, 3, 4, 5};

  std::unordered_map<source_type, uint32_t> results{
    {static_cast<source_type>(1), 1},
    {static_cast<source_type>(2), 2},
    {static_cast<source_type>(3), 3},
    {static_cast<source_type>(4), 4},
    {static_cast<source_type>(5), 5},
  };

  size_t tableSize = 5;

  frequencyTable_T frequencyTable{};
  frequencyTable.addFrequencies(frequencies.begin(), frequencies.end(), 0);

  frequencyTable_T frequencyTable2{};
  frequencyTable2.addFrequencies(gsl::make_span(frequencies), 0);

  BOOST_CHECK_EQUAL_COLLECTIONS(frequencyTable.begin(), frequencyTable.end(), frequencyTable2.begin(), frequencyTable2.end());

  for (const auto [symbol, value] : results) {
    BOOST_CHECK_EQUAL(frequencyTable[symbol], value);
  }

  BOOST_CHECK_EQUAL(frequencyTable.empty(), false);
  BOOST_CHECK_EQUAL(frequencyTable.size(), tableSize);
  BOOST_CHECK_EQUAL(frequencyTable.getOffset(), 1);
  BOOST_CHECK(frequencyTable.begin() != frequencyTable.end());
  BOOST_CHECK(frequencyTable.cbegin() != frequencyTable.cend());
  BOOST_CHECK_EQUAL(frequencyTable.computeNUsedAlphabetSymbols(), 5);
  BOOST_CHECK_EQUAL(frequencyTable.getNumSamples(), 15);

  // lets add more frequencies;
  std::vector<value_type> frequencies2{3, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0};

  if constexpr (std::is_signed_v<source_type>) {
    frequencyTable.addFrequencies(frequencies2.begin(), frequencies2.end(), -1);
    frequencyTable2.addFrequencies(gsl::make_span(frequencies2), -1);

    results[static_cast<source_type>(0 + -1)] += 3;
    results[static_cast<source_type>(2 + -1)] += 4;
    results[static_cast<source_type>(11 + -1)] += 5;

    BOOST_CHECK_EQUAL(frequencyTable.size(), 12);
    BOOST_CHECK_EQUAL(frequencyTable.getOffset(), -1);
    BOOST_CHECK_EQUAL(frequencyTable.computeNUsedAlphabetSymbols(), 7);
  } else {
    frequencyTable.addFrequencies(frequencies2.begin(), frequencies2.end(), 3);
    frequencyTable2.addFrequencies(gsl::make_span(frequencies2), 3);

    results[static_cast<source_type>(0 + 3)] += 3;
    results[static_cast<source_type>(2 + 3)] += 4;
    results[static_cast<source_type>(11 + 3)] += 5;

    BOOST_CHECK_EQUAL(frequencyTable.size(), 14);
    BOOST_CHECK_EQUAL(frequencyTable.getOffset(), 1);
    BOOST_CHECK_EQUAL(frequencyTable.computeNUsedAlphabetSymbols(), 6);
  }
  BOOST_CHECK_EQUAL(frequencyTable.getNumSamples(), 27);
  BOOST_CHECK_EQUAL_COLLECTIONS(frequencyTable.begin(), frequencyTable.end(), frequencyTable2.begin(), frequencyTable2.end());

  for (const auto [symbol, value] : results) {
    BOOST_CHECK_EQUAL(frequencyTable[symbol], value);
  }

  BOOST_CHECK_EQUAL(frequencyTable.empty(), false);
  BOOST_CHECK(frequencyTable.begin() != frequencyTable.end());
  BOOST_CHECK(frequencyTable.cbegin() != frequencyTable.cend());
};

BOOST_AUTO_TEST_CASE_TEMPLATE(test_addFrequenciesHash, frequencyTable_T, hashFrequencyTables_t)
{
  using source_type = typename frequencyTable_T::source_type;
  using value_type = typename frequencyTable_T::value_type;
  std::vector<value_type> frequencies{0, 1, 2, 3, 4, 5};

  std::unordered_map<source_type, uint32_t> results{
    {static_cast<source_type>(0), 0},
    {static_cast<source_type>(1), 1},
    {static_cast<source_type>(2), 2},
    {static_cast<source_type>(3), 3},
    {static_cast<source_type>(4), 4},
    {static_cast<source_type>(5), 5},
  };

  frequencyTable_T frequencyTable{};
  frequencyTable.addFrequencies(frequencies.begin(), frequencies.end(), 0);

  frequencyTable_T frequencyTable2{};
  frequencyTable2.addFrequencies(gsl::make_span(frequencies), 0);

  BOOST_CHECK_EQUAL(frequencyTable.size(), frequencyTable2.size());
  for (const auto [symbol, value] : frequencyTable) {
    BOOST_CHECK_EQUAL(frequencyTable2[symbol], value);
  }

  for (const auto [symbol, value] : results) {
    BOOST_CHECK_EQUAL(frequencyTable[symbol], value);
  }

  BOOST_CHECK_EQUAL(frequencyTable.empty(), false);
  BOOST_CHECK_EQUAL(frequencyTable.size(), 5);
  BOOST_CHECK(frequencyTable.begin() != frequencyTable.end());
  BOOST_CHECK(frequencyTable.cbegin() != frequencyTable.cend());
  BOOST_CHECK_EQUAL(frequencyTable.computeNUsedAlphabetSymbols(), 5);
  BOOST_CHECK_EQUAL(frequencyTable.getNumSamples(), 15);

  // lets add more frequencies;
  std::vector<value_type> frequencies2{3, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0};

  if constexpr (std::is_signed_v<source_type>) {
    frequencyTable.addFrequencies(frequencies2.begin(), frequencies2.end(), -1);
    frequencyTable2.addFrequencies(gsl::make_span(frequencies2), -1);

    results[static_cast<source_type>(0 + -1)] += 3;
    results[static_cast<source_type>(2 + -1)] += 4;
    results[static_cast<source_type>(11 + -1)] += 5;

    BOOST_CHECK_EQUAL(frequencyTable.computeNUsedAlphabetSymbols(), 7);
    BOOST_CHECK_EQUAL(frequencyTable.size(), 7);
  } else {
    frequencyTable.addFrequencies(frequencies2.begin(), frequencies2.end(), 3);
    frequencyTable2.addFrequencies(gsl::make_span(frequencies2), 3);

    results[static_cast<source_type>(0 + 3)] += 3;
    results[static_cast<source_type>(2 + 3)] += 4;
    results[static_cast<source_type>(11 + 3)] += 5;

    BOOST_CHECK_EQUAL(frequencyTable.computeNUsedAlphabetSymbols(), 6);
    BOOST_CHECK_EQUAL(frequencyTable.size(), 6);
  }
  BOOST_CHECK_EQUAL(frequencyTable.getNumSamples(), 27);

  BOOST_CHECK_EQUAL(frequencyTable.size(), frequencyTable2.size());
  for (const auto [symbol, value] : frequencyTable) {
    BOOST_CHECK_EQUAL(frequencyTable2[symbol], value);
  }

  for (const auto [symbol, value] : results) {
    BOOST_CHECK_EQUAL(frequencyTable[symbol], value);
  }

  BOOST_CHECK_EQUAL(frequencyTable.empty(), false);
  BOOST_CHECK(frequencyTable.begin() != frequencyTable.end());
  BOOST_CHECK(frequencyTable.cbegin() != frequencyTable.cend());
};

using frequencyContainer_t = boost::mpl::vector<
  StaticFrequencyTable<uint8_t>,
  DynamicFrequencyTable<uint8_t>,
  HashFrequencyTable<uint8_t>,
  FrequencyTable>;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_renormIncompressible, frequencyTable_T, frequencyContainer_t)
{
  histogram_t frequencies{1, 1, 2, 2, 2, 2, 6, 8, 4, 10, 8, 14, 10, 19, 26, 30, 31, 35, 41, 45, 51, 44, 47, 39, 58, 52, 42, 53, 50, 34, 50, 30, 32, 24, 30, 20, 17, 12, 16, 6, 8, 5, 6, 4, 4, 2, 2, 2, 1};
  frequencyTable_T frequencyTable{frequencies.begin(), frequencies.end(), static_cast<uint8_t>(0)};

  const size_t scaleBits = 8;

  auto renormedFrequencyTable = renormCutoffIncompressible<frequencyTable_T>(std::move(frequencyTable), scaleBits, 1);

  const o2::rans::histogram_t rescaledFrequencies{1, 2, 1, 3, 2, 3, 3, 5, 6, 7, 8, 9, 10, 11, 13, 11, 12, 10, 14, 13, 10, 13, 12, 8, 12, 7, 8, 6, 7, 5, 4, 3, 4, 2, 2, 1, 2, 1, 1};
  BOOST_CHECK_EQUAL(renormedFrequencyTable.isRenormedTo(scaleBits), true);
  BOOST_CHECK_EQUAL(renormedFrequencyTable.getNumSamples(), 1 << scaleBits);
  BOOST_CHECK_EQUAL(renormedFrequencyTable.getIncompressibleSymbolFrequency(), 4);

  if constexpr (std::is_same_v<frequencyTable_T, HashFrequencyTable<uint8_t>>) {
    for (uint8_t index = 0; index < rescaledFrequencies.size(); ++index) {
      BOOST_CHECK_EQUAL(renormedFrequencyTable[index + 6], rescaledFrequencies[index]);
    }
  } else if constexpr (std::is_same_v<frequencyTable_T, FrequencyTable>) {
    BOOST_CHECK_EQUAL_COLLECTIONS(renormedFrequencyTable.begin(), renormedFrequencyTable.end(), rescaledFrequencies.begin(), rescaledFrequencies.end());
  } else {
    BOOST_CHECK_EQUAL_COLLECTIONS(renormedFrequencyTable.begin() + 6, renormedFrequencyTable.begin() + 6 + rescaledFrequencies.size(), rescaledFrequencies.begin(), rescaledFrequencies.end());
  }
}