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
/// @brief  Lookup table containing statistical information for each symbol in the alphabet required for encoding/decoding

#ifndef RANS_INTERNAL_CONTAINERS_SYMBOLTABLE_H_
#define RANS_INTERNAL_CONTAINERS_SYMBOLTABLE_H_

#include <vector>
#include <cstdint>
#include <cmath>
#include <fairlogger/Logger.h>

#include "rANS/internal/containers/Container.h"
#include "rANS/internal/containers/RenormedHistogram.h"
#include "rANS/internal/containers/HistogramView.h"

namespace o2::rans
{

template <class source_T, class symbol_T>
class SymbolTable : public internal::Container<source_T, symbol_T, SymbolTable<source_T, symbol_T>>
{
  using base_type = internal::Container<source_T, symbol_T, SymbolTable<source_T, symbol_T>>;
  friend base_type;

 public:
  using source_type = typename base_type::source_type;
  using symbol_type = typename base_type::value_type;
  using container_type = typename base_type::container_type;
  using size_type = typename base_type::size_type;
  using difference_type = typename base_type::difference_type;
  using reference = typename base_type::reference;
  using const_reference = typename base_type::const_reference;
  using pointer = typename base_type::pointer;
  using const_pointer = typename base_type::const_pointer;
  using const_iterator = typename base_type::const_iterator;

  SymbolTable() = default;

  SymbolTable(const RenormedHistogram<source_type>& renormedHistogram);

  [[nodiscard]] inline const_reference operator[](source_type sourceSymbol) const noexcept
  {
    const size_type index = static_cast<size_type>(sourceSymbol - this->getOffset());
    // static cast to unsigned: idx < 0 => (uint)idx > MAX_INT => idx > mIndex.size()
    if (index < this->size()) {
      return this->mContainer[sourceSymbol];
    } else {
      return this->getEscapeSymbol();
    }
  };

  [[nodiscard]] inline const_pointer lookupSafe(source_type sourceSymbol) const
  {
    const size_type index = static_cast<size_type>(sourceSymbol - this->getOffset());
    // static cast to unsigned: idx < 0 => (uint)idx > MAX_INT => idx > mIndex.size()
    if (index < this->size()) {
      return this->mContainer.data() + index;
    } else {
      return nullptr;
    }
  };

  [[nodiscard]] inline const_pointer lookupUnsafe(source_type sourceSymbol) const
  {
    return &this->mContainer[sourceSymbol];
  };

  [[nodiscard]] inline size_type size() const noexcept { return mSize; };

  [[nodiscard]] inline const_reference getEscapeSymbol() const noexcept { return mEscapeSymbol; };

  [[nodiscard]] inline bool isEscapeSymbol(const_reference symbol) const noexcept { return symbol == mEscapeSymbol; };

  [[nodiscard]] inline bool isEscapeSymbol(source_type sourceSymbol) const noexcept { return this->operator[](sourceSymbol) == mEscapeSymbol; };

  [[nodiscard]] inline size_type getPrecision() const noexcept { return mSymbolTablePrecision; };

  friend void swap(SymbolTable& a, SymbolTable& b) noexcept
  {
    using std::swap;
    swap(static_cast<typename SymbolTable::base_type&>(a),
         static_cast<typename SymbolTable::base_type&>(b));
    swap(a.mSize, b.mSize);
    swap(a.mEscapeSymbol, b.mEscapeSymbol);
    swap(a.mSymbolTablePrecision, b.mSymbolTablePrecision);
  }

 protected:
  [[nodiscard]] inline bool isValidSymbol(const symbol_type& value) const noexcept
  {
    return !this->isEscapeSymbol(value);
  };

  size_t mSize = 0;
  symbol_type mEscapeSymbol{};
  size_type mSymbolTablePrecision{};
};

template <class source_T, class value_T>
SymbolTable<source_T, value_T>::SymbolTable(const RenormedHistogram<source_type>& histogram)
{
  using count_type = typename value_T::value_type;

  auto histogramView = internal::trim(internal::HistogramView(histogram.begin(), histogram.end(), histogram.getOffset()));

  this->mContainer.reserve(histogramView.size());
  this->mSymbolTablePrecision = histogram.getRenormingBits();
  this->mEscapeSymbol = [&]() -> value_T {
    const count_type symbolFrequency = histogram.getIncompressibleSymbolFrequency();
    const count_type cumulatedFrequency = histogram.getNumSamples() - symbolFrequency;
    return {symbolFrequency, cumulatedFrequency, this->getPrecision()};
  }();

  count_type cumulatedFrequency = 0;
  for (const auto symbolFrequency : histogramView) {
    if (symbolFrequency) {
      this->mContainer.emplace_back(symbolFrequency, cumulatedFrequency, this->getPrecision());
      cumulatedFrequency += symbolFrequency;
    } else {
      this->mContainer.push_back(this->mEscapeSymbol);
    }
  }
  this->mContainer.setOffset(histogramView.getOffset());
  mSize = this->mContainer.size();
};

} // namespace o2::rans

#endif /* RANS_INTERNAL_CONTAINERS_SYMBOLTABLE_H_ */
