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

#ifndef RANS_SYMBOLTABLECONTAINER_H
#define RANS_SYMBOLTABLECONTAINER_H

#include <vector>
#include <cstdint>
#include <cmath>
#include <fairlogger/Logger.h>

#include "rANS/definitions.h"
#include "rANS/RenormedFrequencyTable.h"
#include "rANS/internal/backend/simd/Symbol.h"

#include "rANS/internal/ContainerInterface.h"

namespace o2
{
namespace rans
{
namespace internal
{

template <class source_T, class index_T, class value_T, class container_T, class derived_T>
class SymbolTableContainer : public ContainerInterface<source_T, index_T, value_T, container_T,
                                                       SymbolTableContainer<source_T, index_T, value_T, container_T, derived_T>>
{
  using base_type = ContainerInterface<source_T, index_T, value_T, container_T,
                                       SymbolTableContainer<source_T, index_T, value_T, container_T, derived_T>>;

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

  [[nodiscard]] inline const_reference operator[](source_type sourceSymbol) const { return static_cast<const derived_T*>(this)->operator[](sourceSymbol); };

  [[nodiscard]] inline size_type size() const noexcept { return static_cast<derived_T*>(this)->size(); };

  [[nodiscard]] inline size_type computeNUsedAlphabetSymbols() const noexcept { return static_cast<const derived_T*>(this)->computeNUsedAlphabetSymbols(); };

  [[nodiscard]] inline const_reference getEscapeSymbol() const noexcept { return mEscapeSymbol; };

  [[nodiscard]] inline bool isEscapeSymbol(const_reference symbol) const noexcept { return symbol == mEscapeSymbol; };

  [[nodiscard]] inline bool isEscapeSymbol(source_type sourceSymbol) const noexcept { return this->operator[](sourceSymbol) == mEscapeSymbol; };

  [[nodiscard]] inline size_type getPrecision() const noexcept { return mSymbolTablePrecision; };

  friend void swap(SymbolTableContainer& a, SymbolTableContainer& b) noexcept
  {
    using std::swap;
    swap(a.mEscapeSymbol, b.mEscapeSymbol);
    swap(a.mSymbolTablePrecision, b.mSymbolTablePrecision);
    swap(static_cast<typename SymbolTableContainer::base_type&>(a),
         static_cast<typename SymbolTableContainer::base_type&>(b));
  };

 protected:
  SymbolTableContainer() = default;

  value_type mEscapeSymbol{};
  size_type mSymbolTablePrecision{};
};

} // namespace internal
} // namespace rans
} // namespace o2

#endif /* RANS_SYMBOLTABLECONTAINER_H */
