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

#ifndef RANS_INTERNAL_TRANSFORM_HASHALGORITHM_H_
#define RANS_INTERNAL_TRANSFORM_HASHALGORITHM_H_

#include <array>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <cassert>

#include "rANS/internal/common/containertraits.h"
#include "rANS/internal/common/utils.h"
#include "rANS/internal/transform/algorithm.h"

namespace o2::rans::internal
{

template <class IT, std::enable_if_t<isPair_v<typename std::iterator_traits<IT>::value_type>, bool> = true>
inline auto getValue(IT iter) -> typename std::iterator_traits<IT>::value_type::second_type
{
  return iter->second;
}

template <class IT, std::enable_if_t<isPair_v<typename std::iterator_traits<IT>::value_type>, bool> = true>
inline void setValue(IT iter, const typename std::iterator_traits<IT>::value_type::second_type& value)
{
  return iter->second = value;
}

template <typename source_T, typename value_T>
inline auto getValue(const std::pair<source_T, value_T>& pair) -> value_T
{
  return pair.second;
}

template <typename container_T, std::enable_if_t<isHashContainer_v<container_T>, bool> = true>
inline decltype(auto) getValue(const typename container_T::container_type::value_type& pair)
{
  return pair.second;
}

template <class container_T, class IT,
          std::enable_if_t<std::is_same_v<IT, typename container_T::const_iterator> && isPair_v<typename std::iterator_traits<typename container_T::const_iterator>::value_type>, bool> = true>
inline constexpr auto getIndex(const container_T& histogram, IT iter) -> typename container_T::source_type
{
  return iter->first;
};

template <typename container_T, class F, std::enable_if_t<isHashContainer_v<container_T>, bool> = true>
inline void forEachIndexValue(const container_T& container, typename container_T::const_iterator begin, typename container_T::const_iterator end, F functor)
{
  using iterator_type = typename container_T::const_iterator;

  std::vector<iterator_type> orderedIterators{};
  orderedIterators.reserve(container.size());
  for (auto iter = begin; iter != end; ++iter) {
    orderedIterators.push_back(iter);
  };

  std::sort(orderedIterators.begin(), orderedIterators.end(), [](const iterator_type& a, const iterator_type& b) {
    return a->first < b->first;
  });

  for (const auto& iter : orderedIterators) {
    functor(iter->first, iter->second);
  }
};

template <class IT, std::enable_if_t<isPair_v<typename std::iterator_traits<IT>::value_type>, bool> = true>
inline auto trim(IT begin, IT end,
                 const typename std::iterator_traits<IT>::value_type::second_type& zeroElem = {})
  -> std::pair<IT, IT>
{
  return {begin, end};
};

template <typename container_T, std::enable_if_t<isHashContainer_v<container_T>, bool> = true>
decltype(auto) getMinMax(const container_T& container, typename container_T::const_iterator begin, typename container_T::const_iterator end)
{
  using iterator_type = typename container_T::const_iterator;
  using value_type = typename std::iterator_traits<iterator_type>::value_type::second_type;
  using return_type = std::pair<value_type, value_type>;

  bool empty = container.empty();

  if constexpr (isRenormedHistogram_v<container_T>) {
    empty = container.getNumSamples() == container.getIncompressibleSymbolFrequency();
  }

  if (empty) {
    return return_type{container.getOffset(), container.getOffset()};
  } else {
    const auto [minIter, maxIter] = std::minmax_element(begin, end, [](const auto& a, const auto& b) { return a.first < b.first; });
    return return_type{minIter->first, maxIter->first};
  }
};

template <typename container_T, std::enable_if_t<isHashContainer_v<container_T>, bool> = true>
auto getMinMax(const container_T& container) -> std::pair<typename container_T::source_type,
                                                          typename container_T::source_type>
{
  auto [trimmedBegin, trimmedEnd] = trim(container.begin(), container.end());
  return getMinMax(container, container.begin(), container.end());
};

template <typename IT, class F, std::enable_if_t<isPair_v<typename std::iterator_traits<IT>::value_type>, bool> = true>
inline void forEachValue(IT begin, IT end, F functor)
{
  for (auto iter = begin; iter != end; ++iter) {
    functor(iter->second);
  }
};

template <typename container_T, class F, std::enable_if_t<isHashContainer_v<container_T>, bool> = true>
inline void forEachIndexValue(const container_T& container, F functor)
{
  return forEachIndexValue(container, container.begin(), container.end(), functor);
};

} // namespace o2::rans::internal

#endif /* RANS_INTERNAL_TRANSFORM_HASHALGORITHM_H_ */