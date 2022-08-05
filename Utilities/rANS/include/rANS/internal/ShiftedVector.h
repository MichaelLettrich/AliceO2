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

#ifndef INCLUDE_RANS_SHIFTEDVECTOR_H_
#define INCLUDE_RANS_SHIFTEDVECTOR_H_

#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cassert>

#include <rANS/internal/helper.h>

#include <fairlogger/Logger.h>

namespace o2
{
namespace rans
{

namespace internal
{

template <class source_T, class value_T>
class ShiftedVector
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
  using iterator = pointer;
  using const_iterator = const_pointer;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  ShiftedVector() = default;

  ShiftedVector(ShiftedVector&&) = default;

  ShiftedVector(const ShiftedVector& vec)
  {
    this->mContainer = vec.mContainer;
    this->setOffset(vec.getOffset());
  }

  ShiftedVector& operator=(ShiftedVector&&) = default;

  ShiftedVector& operator=(const ShiftedVector& vec)
  {
    this->mContainer = vec.mContainer;
    this->setOffset(vec.getOffset());
    return *this;
  }

  ShiftedVector(size_type size, source_type offset = 0) : mContainer(size) { this->setOffset(offset); };

  ShiftedVector(container_type container, source_type offset = 0) : mContainer{container} { this->setOffset(offset); };

  ~ShiftedVector() = default;

  [[nodiscard]] inline const_reference operator[](source_type sourceSymbol) const
  {
    assert(static_cast<size_type>(sourceSymbol) - this->getOffset() < this->size());
    return *this->getAddressAt(sourceSymbol);
  };

  [[nodiscard]] inline reference operator[](source_type sourceSymbol) { return const_cast<reference>(static_cast<const ShiftedVector&>(*this)[sourceSymbol]); };

  [[nodiscard]] inline const_reference operator()(size_type index) const
  {
    assert(index < this->size());
    return mContainer[index];
  };

  [[nodiscard]] inline reference operator()(size_type index) { return const_cast<reference>(static_cast<const ShiftedVector&>(*this)(index)); };

  [[nodiscard]] inline const_pointer data() const noexcept { return mContainer.data(); };

  [[nodiscard]] inline pointer data() noexcept { return mContainer.data(); };

  [[nodiscard]] inline size_type size() const noexcept { return mContainer.size(); };

  [[nodiscard]] inline bool empty() const noexcept { return mContainer.empty(); };

  [[nodiscard]] inline source_type getOffset() const noexcept { return mOffset; };

  inline void setOffset(source_type offset) noexcept
  {
    mOffset = offset;
    mShiftedBegin = reinterpret_cast<intptr_t>(mContainer.data()) - mOffset * sizeof(value_type); // circumvent undefined behavior
    assert(static_cast<difference_type>(this->size() + this->getOffset()) <= static_cast<difference_type>(std::numeric_limits<source_type>::max()) + 1);
    assert(this->getAddressAt(offset) == mContainer.data());
  };

  inline void reserve(size_type newSize)
  {
    mContainer.reserve(newSize);
    this->setOffset(this->getOffset());
  };

  inline void resize(size_type newSize)
  {
    assert(newSize <= pow2(toBits(sizeof(source_type))));
    this->resize(newSize, this->getOffset());
  };

  inline void resize(size_type newSize, source_type offset)
  {
    this->resize(newSize);
    this->setOffset(offset);
  };

  inline void push_back(value_type value)
  {
    mContainer.push_back(std::move(value));
    this->setOffset(this->getOffset()); // update in case of reallocation
  };

  template <class... Args>
  inline void emplace_back(Args&&... args)
  {
    mContainer.emplace_back(std::forward<Args>(args)...);
    this->setOffset(this->getOffset()); // update in case of reallocation
  };

  [[nodiscard]] inline const_iterator cbegin() const noexcept { return this->data(); };

  [[nodiscard]] inline const_iterator cend() const noexcept { return this->data() + this->size(); };

  [[nodiscard]] inline const_iterator begin() const noexcept { return cbegin(); };

  [[nodiscard]] inline const_iterator end() const noexcept { return cend(); };

  [[nodiscard]] inline iterator begin() noexcept { return const_cast<iterator>(static_cast<const ShiftedVector&>(*this).begin()); };

  [[nodiscard]] inline iterator end() noexcept { return const_cast<iterator>(static_cast<const ShiftedVector&>(*this).end()); };

  [[nodiscard]] inline const_reverse_iterator crbegin() const noexcept { return std::reverse_iterator{this->cend()}; };

  [[nodiscard]] inline const_reverse_iterator crend() const noexcept { return std::reverse_iterator{this->cbegin()}; };

  [[nodiscard]] inline const_reverse_iterator rbegin() const noexcept { return crbegin(); };

  [[nodiscard]] inline const_reverse_iterator rend() const noexcept { return crend(); };

  [[nodiscard]] inline reverse_iterator rbegin() noexcept { return std::reverse_iterator{this->end()}; };

  [[nodiscard]] inline reverse_iterator rend() noexcept { return std::reverse_iterator{this->begin()}; };

  [[nodiscard]] inline container_type release() && noexcept { return std::move(this->mContainer); };

  friend void swap(ShiftedVector& a, ShiftedVector& b) noexcept
  {
    using std::swap;
    swap(a.mContainer, b.mContainer);
    swap(a.mOffset, b.mOffset);
    swap(a.mShiftedBegin, b.mShiftedBegin);
  };

 protected:
  [[nodiscard]] inline const_pointer getAddressAt(source_type sourceSymbol) const
  {
    return reinterpret_cast<const value_type*>(mShiftedBegin + sourceSymbol * sizeof(value_type)); // circumvent undefined behavior
  }

  [[nodiscard]] inline pointer getAddressAt(source_type sourceSymbol)
  {
    return const_cast<pointer>(static_cast<const ShiftedVector&>(*this).getAddressAt(sourceSymbol));
  }

  container_type mContainer{};
  source_type mOffset{};
  intptr_t mShiftedBegin{};
};
} // namespace internal
} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_SHIFTEDVECTOR_H_ */
