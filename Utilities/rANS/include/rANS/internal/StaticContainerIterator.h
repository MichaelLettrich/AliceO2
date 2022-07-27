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

/// @file   FrequencyTableImpl.h
/// @author Michael Lettrich
/// @since  2019-05-08
/// @brief Histogram to depict frequencies of source symbols for rANS compression.

#ifndef INCLUDE_RANS_STATICCONTAINERITERATOR_H_
#define INCLUDE_RANS_STATICCONTAINERITERATOR_H_

#include <numeric>
#include <vector>

#include "rANS/definitions.h"
#include "rANS/internal/helper.h"
#include "rANS/internal/FrequencyContainer.h"

namespace o2
{
namespace rans
{
namespace internal
{

template <class staticContainer_T>
class StaticContainerIterator
{
 public:
  using difference_type = typename staticContainer_T::difference_type;
  using value_type = typename staticContainer_T::value_type;
  using pointer = typename staticContainer_T::pointer;
  using reference = typename staticContainer_T::reference;
  using iterator_category = std::random_access_iterator_tag;

  using index_type = int64_t;

  inline constexpr StaticContainerIterator() noexcept = default;
  inline constexpr StaticContainerIterator(staticContainer_T* container, index_type idx) noexcept : mIndex{idx}, mContainer{container} {};
  inline constexpr StaticContainerIterator(const StaticContainerIterator& iter) noexcept = default;
  inline constexpr StaticContainerIterator(StaticContainerIterator&& iter) noexcept = default;
  inline constexpr StaticContainerIterator& operator=(const StaticContainerIterator& other) noexcept = default;
  inline constexpr StaticContainerIterator& operator=(StaticContainerIterator&& other) noexcept = default;
  inline ~StaticContainerIterator() noexcept = default;

  // pointer arithmetics
  inline constexpr StaticContainerIterator& operator++() noexcept
  {
    ++mIndex;
    return *this;
  };

  inline constexpr StaticContainerIterator operator++(int) noexcept
  {
    auto res = *this;
    ++(*this);
    return res;
  };

  inline constexpr StaticContainerIterator& operator--() noexcept
  {
    --mIndex;
    return *this;
  };

  inline constexpr StaticContainerIterator operator--(int) noexcept
  {
    auto res = *this;
    --(*this);
    return res;
  };

  inline constexpr StaticContainerIterator& operator+=(difference_type i) noexcept
  {
    mIndex += i;
    return *this;
  };

  inline constexpr StaticContainerIterator operator+(difference_type i) const noexcept
  {
    auto tmp = *const_cast<StaticContainerIterator*>(this);
    return tmp += i;
  }

  inline constexpr StaticContainerIterator& operator-=(difference_type i) noexcept
  {
    mIndex -= i;
    return *this;
  };

  inline constexpr StaticContainerIterator operator-(difference_type i) const noexcept
  {
    auto tmp = *const_cast<StaticContainerIterator*>(this);
    return tmp -= i;
  };

  inline constexpr difference_type operator-(const StaticContainerIterator& other) const noexcept
  {
    return this->mIter - other.mIter;
  };

  // comparison
  inline constexpr bool operator==(const StaticContainerIterator& other) const noexcept { return this->mIndex == other.mIndex; };
  inline constexpr bool operator!=(const StaticContainerIterator& other) const noexcept { return this->mIndex != other.mIndex; };
  inline constexpr bool operator<(const StaticContainerIterator& other) const noexcept { return this->mIndex < other->mIndex; };
  inline constexpr bool operator>(const StaticContainerIterator& other) const noexcept { return this->mIndex > other->mIndex; };
  inline constexpr bool operator>=(const StaticContainerIterator& other) const noexcept { return this->mIndex >= other->mIndex; };
  inline constexpr bool operator<=(const StaticContainerIterator& other) const noexcept { return this->mIndex <= other->mIndex; };

  // dereference
  inline constexpr decltype(auto) operator*() const { return (*mContainer)[mIndex]; };
  inline constexpr decltype(auto) operator*() { return (*mContainer)[mIndex]; };

  inline constexpr const reference& operator[](difference_type i) const noexcept { return *(*this + i); };

 private:
  index_type mIndex{};
  staticContainer_T* mContainer;
};

} // namespace internal
} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_STATICCONTAINERITERATOR_H_ */