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

/// @file   test_ransSymbolStatistics.cxx
/// @author Michael Lettrich
/// @since  2021-06-02
/// @brief

#define BOOST_TEST_MODULE Utility test
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/mpl/vector.hpp>
#include <vector>
#include <numeric>
#include <iterator>

#include "rANS/FrequencyTable.h"

#include "rANS/renorm.h"

#include "rANS/internal/Symbol.h"
#include "rANS/SymbolTable.h"

using namespace o2::rans;

template <typename T>
size_t getNUniqueSymbols(const T& container)
{
  return std::count_if(container.begin(), container.end(), [](uint32_t value) { return value != 0; });
};

using frequencyTable_t = boost::mpl::vector<FrequencyTable<uint8_t>,
                                            FrequencyTable<int8_t>,
                                            FrequencyTable<uint16_t>,
                                            FrequencyTable<int16_t>,
                                            FrequencyTable<uint32_t>,
                                            FrequencyTable<int32_t>>;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_empty, frequencyTable_T, frequencyTable_t)
{
  using namespace o2::rans;

  using source_type = typename frequencyTable_T::source_type;
  std::vector<source_type> samples{};
  frequencyTable_T f{};
  f.addSamples(gsl::make_span(samples));

  const auto symbolTable = SymbolTable<source_type, internal::Symbol>(renormCutoffIncompressible(f));

  BOOST_CHECK_EQUAL(symbolTable.size(), 0);
  BOOST_CHECK_EQUAL(symbolTable.countNUsedAlphabetSymbols(), 0);

  const auto escapeSymbol = symbolTable.getEscapeSymbol();
  BOOST_CHECK_EQUAL(escapeSymbol.getFrequency(), 1);
  BOOST_CHECK_EQUAL(escapeSymbol.getCumulative(), 0);
  BOOST_CHECK_EQUAL(symbolTable.isEscapeSymbol(escapeSymbol), true);

  //test in range
  for (const auto& symbol : symbolTable) {
    BOOST_CHECK_EQUAL(symbolTable.isEscapeSymbol(symbol), true);
  }

  // out of range checks:
  const int outOfRangeSymbols[] = {-100, 100};
  for (auto symbol : outOfRangeSymbols) {
    const auto outOfRange = symbolTable[symbol];
    BOOST_CHECK_EQUAL(symbolTable.isEscapeSymbol(symbol), true);
    BOOST_CHECK_EQUAL(outOfRange.getFrequency(), escapeSymbol.getFrequency());
    BOOST_CHECK_EQUAL(outOfRange.getCumulative(), escapeSymbol.getCumulative());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_symbolTable, frequencyTable_T, frequencyTable_t)
{
  using source_type = typename frequencyTable_T::source_type;

  histogram_t frequencies{1, 1, 2, 2, 2, 2, 6, 8, 4, 10, 8, 14, 10, 19, 26, 30, 31, 35, 41, 45, 51, 44, 47, 39, 58, 52, 42, 53, 50, 34, 50, 30, 32, 24, 30, 20, 17, 12, 16, 6, 8, 5, 6, 4, 4, 2, 2, 2, 1};
  frequencyTable_T frequencyTable{frequencies.begin(), frequencies.end(), static_cast<uint8_t>(0)};
  const size_t scaleBits = 8;
  auto renormedFrequencyTable = renormCutoffIncompressible(std::move(frequencyTable), scaleBits, 1);
  const auto symbolTable = SymbolTable<source_type, internal::Symbol>(renormedFrequencyTable);

  const histogram_t rescaledFrequencies{1, 2, 1, 3, 2, 3, 3, 5, 6, 7, 8, 9, 10, 11, 13, 11, 12, 10, 14, 13, 10, 13, 12, 8, 12, 7, 8, 6, 7, 5, 4, 3, 4, 2, 2, 1, 2, 1, 1};
  histogram_t cumulativeFrequencies;
  std::exclusive_scan(rescaledFrequencies.begin(), rescaledFrequencies.end(), std::back_inserter(cumulativeFrequencies), 0);

  BOOST_CHECK_EQUAL(symbolTable.getPrecision(), scaleBits);

  for (source_type index = 0; index < rescaledFrequencies.size(); ++index) {
    auto symbol = symbolTable[index + 6];
    BOOST_CHECK_EQUAL(symbol.getFrequency(), rescaledFrequencies[index]);
    BOOST_CHECK_EQUAL(symbol.getCumulative(), cumulativeFrequencies[index]);
  }

  const auto escapeSymbol = symbolTable.getEscapeSymbol();
  BOOST_CHECK_EQUAL(escapeSymbol.getFrequency(), 4);
  BOOST_CHECK_EQUAL(escapeSymbol.getCumulative(), cumulativeFrequencies.back() + rescaledFrequencies.back());

  // out of range checks:
  const int outOfRangeSymbols[] = {-100, 100};
  for (auto symbol : outOfRangeSymbols) {
    const auto escapeSymbol = symbolTable.getEscapeSymbol();
    const auto outOfRange = symbolTable[symbol];
    BOOST_CHECK_EQUAL(symbolTable.isEscapeSymbol(symbol), true);
    BOOST_CHECK_EQUAL(outOfRange.getFrequency(), escapeSymbol.getFrequency());
    BOOST_CHECK_EQUAL(outOfRange.getCumulative(), escapeSymbol.getCumulative());
  }
}
