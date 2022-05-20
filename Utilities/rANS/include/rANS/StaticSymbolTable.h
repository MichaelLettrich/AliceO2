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

#ifndef RANS_STATICSYMBOLTABLE_H
#define RANS_STATICSYMBOLTABLE_H

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
class StaticSymbolTable : public internal::SymbolTableContainer<source_T,
                                                                std::make_unsigned_t<source_T>,
                                                                value_T,
                                                                std::vector<value_T>,
                                                                StaticSymbolTable<source_T, value_T>>
{
  using base_type = internal::SymbolTableContainer<source_T,
                                                   std::make_unsigned_t<source_T>,
                                                   value_T,
                                                   std::vector<value_T>,
                                                   StaticSymbolTable<source_T, value_T>>;

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

  StaticSymbolTable() = default;

  StaticSymbolTable(const RenormedStaticFrequencyTable<source_type>& renormedFrequencies);

  static_assert(sizeof(index_type) <= 2, "This datatype requires a <=16Bit datatype for source_T");

  [[nodiscard]] inline const_reference operator[](source_type sourceSymbol) const { return this->mContainer[static_cast<index_type>(sourceSymbol)]; };

  [[nodiscard]] inline const_pointer data() const noexcept { return this->mContainer.data(); };

  [[nodiscard]] inline size_type size() const noexcept { return internal::pow2(internal::toBits(sizeof(index_type))); };

  [[nodiscard]] inline size_type computeNUsedAlphabetSymbols() const noexcept
  {
    return std::count_if(this->begin(), this->end(), [this](const_reference v) { return !this->isEscapeSymbol(v); });
  };

  friend void swap(StaticSymbolTable& a, StaticSymbolTable& b) noexcept
  {
    using std::swap;
    swap(static_cast<typename StaticSymbolTable::base_type&>(a),
         static_cast<typename StaticSymbolTable::base_type&>(b));
  };
}; // namespace rans

template <class source_T, class value_T>
StaticSymbolTable<source_T, value_T>::StaticSymbolTable(const RenormedStaticFrequencyTable<source_type>& frequencyTable) : StaticSymbolTable()
{

  using count_type = typename value_T::value_type;

  this->mContainer.reserve(this->size());
  this->mSymbolTablePrecision = frequencyTable.getRenormingBits();
  this->mEscapeSymbol = [&]() -> value_T {
    const count_type symbolFrequency = frequencyTable.getIncompressibleSymbolFrequency();
    const count_type cumulatedFrequency = static_cast<count_type>(frequencyTable.getNumSamples() - symbolFrequency);
    return {symbolFrequency, cumulatedFrequency, this->getPrecision()};
  }();

  count_type cumulatedFrequency = 0;
  for (const auto symbolFrequency : frequencyTable) {
    if (symbolFrequency) {
      this->mContainer.emplace_back(symbolFrequency, cumulatedFrequency, this->getPrecision());
      cumulatedFrequency += symbolFrequency;
    } else {
      this->mContainer.push_back(this->mEscapeSymbol);
    }
  }
};

template <typename source_T, class symbol_T>
inline StaticSymbolTable<source_T, symbol_T> makeSymbolTable(RenormedStaticFrequencyTable<source_T> table)
{
  return {std::move(table)};
};

} // namespace rans
} // namespace o2

#endif /* RANS_STATICSYMBOLTABLE_H */
