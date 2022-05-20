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

/// @file   SymbolTable.h
/// @author Michael Lettrich
/// @since  2019-06-21
/// @brief  Container for information needed to encode/decode a symbol of the alphabet

#ifndef RANS_HASHSYMBOLTABLE_H
#define RANS_HASHSYMBOLTABLE_H

#include <vector>
#include <cstdint>
#include <cmath>
#include <fairlogger/Logger.h>

#include "rANS/definitions.h"
#include "rANS/RenormedFrequencyTable.h"
#include "rANS/internal/backend/simd/Symbol.h"

#include "rANS/internal/SymbolTableContainer.h"
#include "rANS/RenormedFrequencies.h"

namespace o2
{
namespace rans
{

template <class source_T, class value_T>
class HashSymbolTable : public internal::SymbolTableContainer<source_T,
                                                              std::make_unsigned_t<source_T>,
                                                              value_T,
                                                              absl::flat_hash_map<source_T, value_T>,
                                                              HashSymbolTable<source_T, value_T>>
{
  using base_type = internal::SymbolTableContainer<source_T,
                                                   std::make_unsigned_t<source_T>,
                                                   value_T,
                                                   std::vector<value_T>,
                                                   HashSymbolTable<source_T, value_T>>;

 public:
  using source_type = typename base_type::source_type;
  using index_type = typename base_type::index_type;
  using value_type = typename base_type::value_type;
  using container_type = typename base_type::container_type;
  using size_type = typename base_type::size_type;
  using difference_type = typename base_type::difference_type;
  using reference = typename base_type::reference;
  using const_reference = typename base_type::const_reference;
  using pointer = typename base_type::pointer;
  using const_pointer = typename base_type::const_pointer;
  using const_iterator = typename base_type::const_iterator;

  HashSymbolTable() = default;

  HashSymbolTable(const RenormedHashFrequencyTable<source_type>& renormedFrequencies);

  [[nodiscard]] inline const_reference operator[](source_type sourceSymbol) const
  {
    auto iter = this->mContainer.find(sourceSymbol);
    if (iter == this->mContainer.end()) {
      return this->mEscapeSymbol;
    } else {
      return iter->second;
    }
  };

  [[nodiscard]] inline size_type size() const noexcept { return this->mContainer.size(); };

  [[nodiscard]] inline size_type computeNUsedAlphabetSymbols() const noexcept { return this->size(); };

  friend void swap(HashSymbolTable& a, HashSymbolTable& b) noexcept
  {
    using std::swap;
    swap(static_cast<typename HashSymbolTable::base_type&>(a),
         static_cast<typename HashSymbolTable::base_type&>(b));
  };

 protected:
}; // namespace rans

template <class source_T, class value_T>
HashSymbolTable<source_T, value_T>::HashSymbolTable(const RenormedHashFrequencyTable<source_type>& frequencyTable)
{
  using renormedTable_type = RenormedHashFrequencyTable<source_type>;
  using count_type = typename value_T::value_type;

  std::vector<typename renormedTable_type::const_iterator> orderedIndices = [&, this]() {
    std::vector<typename renormedTable_type::const_iterator> orderedIndices;
    orderedIndices.reserve(frequencyTable.size());
    for (auto iter = frequencyTable.begin(); iter != frequencyTable.end(); ++iter) {
      orderedIndices.push_back(iter);
    }

    std::stable_sort(orderedIndices.begin(), orderedIndices.end(), [](const auto& a, const auto& b) -> bool { return a->first < b->first; });
    return orderedIndices;
  }();

  this->mContainer.reserve(frequencyTable.size());
  this->mSymbolTablePrecision = frequencyTable.getRenormingBits();
  this->mEscapeSymbol = [&]() -> value_T {
    const count_type symbolFrequency = frequencyTable.getIncompressibleSymbolFrequency();
    count_type cumulatedFrequency = frequencyTable.getNumSamples() - symbolFrequency;
    return {symbolFrequency, cumulatedFrequency, this->getPrecision()};
  }();

  count_type cumulatedFrequency = 0;
  for (auto iter : orderedIndices) {
    const auto [symbol, frequency] = *iter;
    this->mContainer.emplace(symbol, value_type{frequency, cumulatedFrequency, this->getPrecision()});
    cumulatedFrequency += frequency;
  }
};

template <typename source_T, class symbol_T>
inline HashSymbolTable<source_T, symbol_T> makeSymbolTable(RenormedHashFrequencyTable<source_T> table)
{
  return {std::move(table)};
};

} // namespace rans
} // namespace o2

#endif /* RANS_HASHSYMBOLTABLE_H */
