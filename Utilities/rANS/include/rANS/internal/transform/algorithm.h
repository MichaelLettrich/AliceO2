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

/// @file   algorithm.h
/// @author Michael Lettrich
/// @brief  helper functionalities useful for packing operations

#ifndef RANS_INTERNAL_TRANSFORM_ALGORITHM_H_
#define RANS_INTERNAL_TRANSFORM_ALGORITHM_H_

#include <array>
#include <cstdint>
#include <cstring>
#include <type_traits>

#include "rANS/internal/common/utils.h"

namespace o2::rans::internal
{

template <typename container_T>
inline constexpr auto getIndex(const container_T& histogram, typename container_T::const_iterator iter) -> typename container_T::source_type
{
  return histogram.getOffset() + std::distance(histogram.begin(), iter);
};

template <typename IT>
inline auto trim(IT begin, IT end, const typename std::iterator_traits<IT>::value_type& zeroElem = {}) -> std::pair<IT, IT>
{
  using value_type = typename std::iterator_traits<IT>::value_type;

  auto isZero = [&zeroElem](value_type i) { return i == zeroElem; };
  auto nonZeroBegin = std::find_if_not(begin, end, isZero);
  auto nonZeroEnd = nonZeroBegin == end ? end : std::find_if_not(std::make_reverse_iterator(end), std::make_reverse_iterator(begin), isZero).base();

  return {nonZeroBegin, nonZeroEnd};
};

template <typename IT, class F>
inline void forEachValue(IT begin, IT end, F functor)
{
  const auto [trimmedBegin, trimmedEnd] = trim(begin, end);
  std::for_each(begin, end, functor);
};

template <class container_T, class F>
inline void forEachValue(const container_T& container, F functor)
{
  forEachValue(container.cbegin(), container.cend(), functor);
};

template <typename container_T, class F>
inline void forEachIndexValue(const container_T& container, typename container_T::const_iterator begin, typename container_T::const_iterator end, F functor)
{
  const auto [trimmedBegin, trimmedEnd] = trim(begin, end);

  typename container_T::source_type index = container.getOffset() + std::distance(container.cbegin(), trimmedBegin);
  for (auto iter = trimmedBegin; iter != trimmedEnd; ++iter) {
    functor(index++, *iter);
  }
};

template <typename container_T, class F>
inline void forEachIndexValue(const container_T& container, F functor)
{
  return forEachIndexValue(container, container.begin(), container.end(), functor);
};

template <class container_T, typename source_IT>
auto getMinMax(const container_T& container, source_IT begin, source_IT end) -> std::pair<typename container_T::source_type,
                                                                                          typename container_T::source_type>
{
  if (begin != end) {
    const auto min = getIndex(container, begin);
    const auto max = getIndex(container, --end);
    assert(max >= min);
    return {min, max};
  } else {
    return {container.getOffset(), container.getOffset()};
  }
};

template <typename container_T>
auto getMinMax(const container_T& container) -> std::pair<typename container_T::source_type,
                                                          typename container_T::source_type>
{
  auto [trimmedBegin, trimmedEnd] = trim(container.begin(), container.end());
  if (trimmedBegin != trimmedEnd) {
    const auto min = getIndex(container, trimmedBegin);
    const auto max = getIndex(container, --trimmedEnd);
    assert(max >= min);
    return {min, max};
  } else {
    return {container.getOffset(), container.getOffset()};
  }
};

} // namespace o2::rans::internal

#endif /* RANS_INTERNAL_TRANSFORM_ALGORITHM_H_ */