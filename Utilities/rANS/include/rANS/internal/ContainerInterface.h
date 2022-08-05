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

/// @file   FrequencyContainer.h
/// @author Michael Lettrich
/// @since  2019-05-08
/// @brief Histogram to depict frequencies of source symbols for rANS compression.

#ifndef INCLUDE_RANS_CONTAINERINTERFACE_H_
#define INCLUDE_RANS_CONTAINERINTERFACE_H_

#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>

namespace o2
{
namespace rans
{

namespace internal
{

template <class source_T, class value_T, class derived_T>
class ContainerInterface
{
 public:
  using source_type = source_T;
  using value_type = value_T;
  using container_type = std::vector<value_type>;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;
  using const_iterator = const value_type*;

  // accessors
  [[nodiscard]] inline const_reference operator[](source_type sourceSymbol) const { return static_cast<const derived_T*>(this)->operator[](sourceSymbol); };

  [[nodiscard]] inline const_pointer data() const noexcept { return this->mContainer.data(); };

  [[nodiscard]] inline const_iterator cbegin() const noexcept { return this->data(); };

  [[nodiscard]] inline const_iterator cend() const noexcept { return this->data() + this->size(); };

  [[nodiscard]] inline const_iterator begin() const noexcept { return cbegin(); };

  [[nodiscard]] inline const_iterator end() const noexcept { return cend(); };

  [[nodiscard]] inline size_type size() const noexcept { return this->mContainer.size(); };

  [[nodiscard]] inline bool empty() const noexcept { return mContainer.empty(); };

  [[nodiscard]] inline source_type getOffset() const noexcept { return mOffset; };

  [[nodiscard]] inline container_type release() && noexcept { return std::move(this->mContainer); };

  friend void swap(ContainerInterface& a, ContainerInterface& b) noexcept
  {
    using std::swap;
    swap(a.mOffset, b.mOffset);
    swap(a.mContainer, b.mContainer);
  };

  [[nodiscard]] inline size_type countNUsedAlphabetSymbols() const noexcept
  {
    return std::count_if(this->begin(), this->end(), [this](const_reference v) { return this->isValidSymbol(v); });
  };

 protected:
  [[nodiscard]] inline bool isValidSymbol(const value_type& value) const noexcept
  {
    return static_cast<const derived_T*>(this)->isValidSymbol(value);
  };

  ContainerInterface() = default;
  ContainerInterface(size_type size, source_type offset) : mOffset(offset), mContainer(size){};

  source_type mOffset{};
  container_type mContainer{};
};
} // namespace internal
} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_CONTAINERINTERFACE_H_ */
