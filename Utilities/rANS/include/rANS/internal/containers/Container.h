// Copyright 2019-2023 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   Container.h
/// @author Michael Lettrich
/// @brief  Abstract container class that defines and  implements basic properties shared by histograms and lookup tables.

#ifndef RANS_INTERNAL_CONTAINERS_CONTAINER_H_
#define RANS_INTERNAL_CONTAINERS_CONTAINER_H_

#include <cstdint>
#include <string>
#include <algorithm>

#include "rANS/internal/containers/ShiftableVector.h"
#include "rANS/internal/containers/SparseVector.h"

namespace o2::rans::internal
{

template <typename container_T>
struct container_traits;

template <typename source_T, typename value_T>
struct container_traits<ShiftableVector<source_T, value_T>> {
  using source_type = typename ShiftableVector<source_T, value_T>::source_type;
  using value_type = typename ShiftableVector<source_T, value_T>::value_type;
};

template <typename source_T, typename value_T>
struct container_traits<SparseVector<source_T, value_T>> {
  using source_type = source_T;
  using value_type = value_T;
};

template <class container_T, class derived_T>
class Container
{
 public:
  using source_type = typename container_traits<container_T>::source_type;
  using value_type = typename container_traits<container_T>::value_type;
  using container_type = container_T;
  using size_type = typename container_type::size_type;
  using difference_type = typename container_type::difference_type;
  using reference = typename container_type::reference;
  using const_reference = typename container_type::const_reference;
  using pointer = typename container_type::pointer;
  using const_pointer = typename container_type::const_pointer;
  using const_iterator = typename container_type::const_iterator;
  using iterator = const_iterator;

  // accessors
  [[nodiscard]] inline const_reference operator[](source_type sourceSymbol) const { return static_cast<const derived_T&>(*this)[sourceSymbol]; };

  [[nodiscard]] inline const_iterator cbegin() const noexcept { return this->mContainer.cbegin(); };

  [[nodiscard]] inline const_iterator cend() const noexcept { return this->mContainer.cend(); };

  [[nodiscard]] inline const_iterator begin() const noexcept { return this->mContainer.begin(); };

  [[nodiscard]] inline const_iterator end() const noexcept { return this->mContainer.end(); };

  [[nodiscard]] inline size_type size() const noexcept { return this->mContainer.size(); };

  [[nodiscard]] inline bool empty() const noexcept { return mNSamples == 0; };

  [[nodiscard]] inline size_type getNumSamples() const noexcept { return mNSamples; };

  [[nodiscard]] inline source_type getOffset() const noexcept { return static_cast<const derived_T*>(this)->getOffset(); };

  [[nodiscard]] inline container_type release() && noexcept { return std::move(this->mContainer); };

 protected:
  template <typename T>
  inline T countSamples(T frequency)
  {
    mNSamples += frequency;
    return frequency;
  };

  Container() = default;
  Container(size_type size, source_type offset) : mContainer{size, offset} {};

  container_type mContainer{};
  size_type mNSamples{};
};

template <typename source_T, typename value_T>
class VectorContainer : public Container<ShiftableVector<source_T, value_T>, VectorContainer<source_T, value_T>>
{
  using base_type = Container<ShiftableVector<source_T, value_T>, VectorContainer<source_T, value_T>>;
  friend base_type;

 public:
  using source_type = typename base_type::source_type;
  using value_type = typename base_type::value_type;
  using container_type = typename base_type::container_type;
  using size_type = typename base_type::size_type;
  using difference_type = typename base_type::difference_type;
  using reference = typename base_type::reference;
  using const_reference = typename base_type::const_reference;
  using pointer = typename base_type::pointer;
  using const_pointer = typename base_type::const_pointer;
  using const_iterator = typename base_type::const_iterator;
  using iterator = typename base_type::iterator;

  [[nodiscard]] inline const_pointer data() const noexcept { return this->mContainer.data(); };

  [[nodiscard]] inline const_reference operator[](source_type sourceSymbol) const { return this->mContainer[sourceSymbol]; };

  [[nodiscard]] inline source_type getOffset() const noexcept { return this->mContainer.getOffset(); };

 protected:
  VectorContainer() = default;
  VectorContainer(size_type size, source_type offset) : base_type{size, offset} {};
};

template <typename source_T, typename value_T>
class SparseVectorContainer : public Container<SparseVector<source_T, value_T>, SparseVectorContainer<source_T, value_T>>
{
  using base_type = Container<SparseVector<source_T, value_T>, SparseVectorContainer<source_T, value_T>>;
  friend base_type;

 public:
  using source_type = typename base_type::source_type;
  using value_type = typename base_type::value_type;
  using container_type = typename base_type::container_type;
  using size_type = typename base_type::size_type;
  using difference_type = typename base_type::difference_type;
  using reference = typename base_type::reference;
  using const_reference = typename base_type::const_reference;
  using pointer = typename base_type::pointer;
  using const_pointer = typename base_type::const_pointer;
  using const_iterator = typename base_type::const_iterator;
  using iterator = typename base_type::iterator;

  [[nodiscard]] inline const_reference operator[](source_type sourceSymbol) const { return this->mContainer[sourceSymbol]; };

  [[nodiscard]] inline const_reference at(source_type sourceSymbol) const { return this->mContainer.at(sourceSymbol); };

  [[nodiscard]] inline source_type getOffset() const noexcept { return this->mContainer.getOffset(); };

 protected:
  SparseVectorContainer() = default;
};

} // namespace o2::rans::internal

#endif /* RANS_INTERNAL_CONTAINERS_CONTAINER_H_ */
