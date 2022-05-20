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

#include "rANS/StaticFrequencyTable.h"
#include "rANS/DynamicFrequencyTable.h"
#include "rANS/HashFrequencyTable.h"

#include "rANS/renorm.h"

#include "rANS/internal/backend/cpp/DecoderSymbol.h"
#include "rANS/StaticSymbolTable.h"
#include "rANS/DynamicSymbolTable.h"
#include "rANS/HashSymbolTable.h"

using namespace o2::rans;

template <typename T>
size_t getNUniqueSymbols(const T& container)
{
  return std::count_if(container.begin(), container.end(), [](uint32_t value) { return value != 0; });
};

using source_type = int8_t;

using frequencyTable_t = boost::mpl::vector<StaticFrequencyTable<source_type>,
                                            DynamicFrequencyTable<source_type>,
                                            HashFrequencyTable<source_type>>;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_empty, frequencyTable_T, frequencyTable_t)
{
  using namespace o2::rans;
  std::vector<source_type> samples{};
  frequencyTable_T f{};
  f.addSamples(gsl::make_span(samples));

  const auto symbolTable = makeSymbolTable<source_type, internal::cpp::DecoderSymbol>(renormCutoffIncompressible(f));

  BOOST_CHECK_EQUAL(symbolTable.size(), f.size());
  BOOST_CHECK_EQUAL(symbolTable.computeNUsedAlphabetSymbols(), 0);

  const auto escapeSymbol = symbolTable.getEscapeSymbol();
  BOOST_CHECK_EQUAL(escapeSymbol.getFrequency(), 1);
  BOOST_CHECK_EQUAL(escapeSymbol.getCumulative(), 0);
  BOOST_CHECK_EQUAL(symbolTable.isEscapeSymbol(escapeSymbol), true);

  //test in range
  for (const auto& symbol : symbolTable) {
    if constexpr (std::is_same_v<frequencyTable_T, HashFrequencyTable<source_type>>) {
      BOOST_CHECK_EQUAL(symbolTable.isEscapeSymbol(symbol.second), true);
    } else {
      BOOST_CHECK_EQUAL(symbolTable.isEscapeSymbol(symbol), true);
    }
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

  histogram_t frequencies{1, 1, 2, 2, 2, 2, 6, 8, 4, 10, 8, 14, 10, 19, 26, 30, 31, 35, 41, 45, 51, 44, 47, 39, 58, 52, 42, 53, 50, 34, 50, 30, 32, 24, 30, 20, 17, 12, 16, 6, 8, 5, 6, 4, 4, 2, 2, 2, 1};
  frequencyTable_T frequencyTable{frequencies.begin(), frequencies.end(), static_cast<uint8_t>(0)};
  const size_t scaleBits = 8;
  auto renormedFrequencyTable = renormCutoffIncompressible<frequencyTable_T>(std::move(frequencyTable), scaleBits, 1);
  const auto symbolTable = makeSymbolTable<source_type, internal::cpp::DecoderSymbol>(renormedFrequencyTable);

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

// BOOST_AUTO_TEST_CASE_TEMPLATE(test_symbolTable, frequencyTable_T, frequencyTable_t)
// {
//   using namespace o2::rans;

//   std::vector<source_type> A{5, 5, 6, 6, 8, 8, 8, 8, 8, -1, -5, 2, 7, 3};
//   std::vector<count_t> histA{1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 2, 2, 1, 5, 1};
//   size_t scaleBits = 17;

//   frequencyTable_T frequencyTable{};
//   frequencyTable.addSamples(gsl::make_span(A));

//   const auto symbolTable = makeSymbolTable<frequencyTable_T>(gsl::make_span(A));

//   BOOST_CHECK_EQUAL(symbolTable.size(), histA.size());
//   BOOST_CHECK_EQUAL(symbolTable.getNUsedAlphabetSymbols(), getNUniqueSymbols(histA));

//   const auto escapeSymbol = symbolTable.getEscapeSymbol();
//   const std::vector<uint32_t> frequencies{8738, 0, 0, 0, 8738, 0, 0, 8738, 8738, 0, 17476, 17477, 8738, 43690, 8739};
//   const std::vector<uint32_t> cumulative{0, 8738, 8738, 8738, 8738, 17476, 17476, 17476, 26214, 34952, 34952, 52428, 69905, 78643, 122333};
//   BOOST_CHECK_EQUAL(symbolTable.size(), frequencies.size());
//   BOOST_CHECK_EQUAL(symbolTable.size(), cumulative.size());

//   // all but last since this is the escape symbol
//   for (size_t i = 0; i < frequencies.size() - 1; ++i) {
//     const uint32_t symbol = min + i;

//     const auto decodeSymbol = symbolTable[symbol];
//     const auto decodeSymbolAt = symbolTable.at(i);
//     BOOST_CHECK_EQUAL(decodeSymbol.getFrequency(), decodeSymbolAt.getFrequency());
//     BOOST_CHECK_EQUAL(decodeSymbol.getCumulative(), decodeSymbolAt.getCumulative());

//     const auto escapeSymbol = symbolTable.getEscapeSymbol();

//     if (frequencies[i] == 0) {
//       BOOST_CHECK_EQUAL(symbolTable.isEscapeSymbol(symbol), true);
//       BOOST_CHECK_EQUAL(decodeSymbol.getFrequency(), escapeSymbol.getFrequency());
//       BOOST_CHECK_EQUAL(decodeSymbol.getCumulative(), escapeSymbol.getCumulative());
//     } else {
//       BOOST_CHECK_EQUAL(symbolTable.isEscapeSymbol(symbol), false);
//       BOOST_CHECK_EQUAL(decodeSymbol.getFrequency(), frequencies[i]);
//       BOOST_CHECK_EQUAL(decodeSymbol.getCumulative(), cumulative[i]);
//     }
//   }
//   // escape symbol:
//   BOOST_CHECK_EQUAL(symbolTable.isEscapeSymbol(0), true);
//   BOOST_CHECK_EQUAL(escapeSymbol.getFrequency(), symbolTable.at(frequencies.size() - 1).getFrequency());
//   BOOST_CHECK_EQUAL(escapeSymbol.getCumulative(), symbolTable.at(frequencies.size() - 1).getCumulative());

//   // out of range checks:
//   const int outOfRangeSymbols[] = {-100, 100};
//   for (auto symbol : outOfRangeSymbols) {
//     const auto escapeSymbol = symbolTable.getEscapeSymbol();
//     const auto outOfRange = symbolTable[symbol];
//     BOOST_CHECK_EQUAL(symbolTable.isEscapeSymbol(symbol), true);
//     BOOST_CHECK_EQUAL(outOfRange.getFrequency(), escapeSymbol.getFrequency());
//     BOOST_CHECK_EQUAL(outOfRange.getCumulative(), escapeSymbol.getCumulative());
//   }
// }
