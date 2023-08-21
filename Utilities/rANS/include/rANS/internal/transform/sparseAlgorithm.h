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

#ifndef RANS_INTERNAL_TRANSFORM_SPARSEALGORITHM_H_
#define RANS_INTERNAL_TRANSFORM_SPARSEALGORITHM_H_

#include <array>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <cassert>

#include "rANS/internal/common/typetraits.h"
#include "rANS/internal/common/utils.h"

namespace o2::rans::internal
{

template <typename container_T, typename C,
          std::enable_if_t<std::is_same_v<typename container_T::container_type, std::remove_cv_t<std::remove_reference_t<C>>>, bool> = true>
inline constexpr auto getIndex(const container_T& histogram, SparseVectorIterator<C> iter) -> typename container_T::source_type
{
  return iter->first;
};

template <class container_T, class F>
inline void forEachValue(SparseVectorIterator<container_T> begin, SparseVectorIterator<container_T> end, F functor)
{
  using iterator_type = SparseVectorIterator<container_T>;
  using lut_iterator = typename iterator_type::lut_iterator;
  using bucket_iterator = typename iterator_type::bucket_iterator;

  // empty
  if (begin == end) {
    return;
  }

  auto lutIter = begin.getLUTIterator();
  if (lutIter != end.getLUTIterator()) {
    //finish first incomplete bucket
    std::for_each(begin.getBucketIterator(), lutIter->end(), functor);

    // go over all remaining buckets
    std::for_each(++lutIter, end.getLUTIterator(), [&functor](const auto& bucket) {
      // and process each element
      for (const auto& elem : bucket) {
        functor(elem);
      };
    });
  }
  // go over the tail, i.e. the last, possibly incomplete or empty bucket
  if (end.getBucketIterator() != bucket_iterator{}) {
    std::for_each(end.getLUTIterator()->begin(), end.getBucketIterator(), functor);
  }
};

template <typename container_T, typename C, class F>
inline void forEachIndexValue(const container_T& container, SparseVectorIterator<C> begin, SparseVectorIterator<C> end, F functor)
{
  using source_type = typename C::source_type;
  using iterator_type = SparseVectorIterator<C>;
  using lut_iterator = typename iterator_type::lut_iterator;
  using bucket_iterator = typename iterator_type::bucket_iterator;

  const C& sparseContainer = begin.getContainer();

  auto makeIndex = [&sparseContainer](lut_iterator lutIter, bucket_iterator bucketIter) -> source_type {
    size_t lut = std::distance(sparseContainer.data(), const_cast<const typename C::bucket_type*>(&(*lutIter)));
    size_t bucket = std::distance(lutIter->begin(), bucketIter);
    return C::joinIndex(lut, bucket);
  };

  // empty
  if (begin == end) {
    return;
  }

  auto lutIter = begin.getLUTIterator();
  if (lutIter != end.getLUTIterator()) {
    //finish first incomplete bucket
    for (auto bucketIter = begin.getBucketIterator(); bucketIter != lutIter->end(); ++bucketIter) {
      functor(makeIndex(lutIter, bucketIter), *bucketIter);
    }

    // go over all remaining buckets
    for (++lutIter; lutIter != end.getLUTIterator(); ++lutIter) {
      // and process each element in each bucket
      for (auto bucketIter = lutIter->begin(); bucketIter != lutIter->end(); ++bucketIter) {
        functor(makeIndex(lutIter, bucketIter), *bucketIter);
      }
    }
  }
  // go over the tail, i.e. the last, possibly incomplete or empty bucket
  if (end.getBucketIterator() != bucket_iterator{}) {
    for (auto bucketIter = lutIter->begin(); bucketIter != end.getBucketIterator(); ++bucketIter) {
      functor(makeIndex(lutIter, bucketIter), *bucketIter);
    }
  }
};

template <class container_T>
inline auto trim(SparseVectorIterator<container_T> begin, SparseVectorIterator<container_T> end,
                 std::add_lvalue_reference_t<std::add_const_t<typename container_T::value_type>> zeroElem = {})
  -> std::pair<SparseVectorIterator<container_T>, SparseVectorIterator<container_T>>
{
  using value_type = typename container_T::value_type;
  using iterator_type = SparseVectorIterator<container_T>;
  using lut_iterator = typename iterator_type::lut_iterator;
  using bucket_iterator = typename iterator_type::bucket_iterator;

  auto isZero = [&zeroElem](const auto& i) { return i == zeroElem; };

  // no range
  if (begin == end) {
    return {end, end};
  }

  iterator_type nonZeroBegin = [&]() -> iterator_type {
    auto lutIter = begin.getLUTIterator();
    if (lutIter != end.getLUTIterator()) {
      //finish first incomplete bucket
      auto nonZeroBegin = std::find_if_not(begin.getBucketIterator(), lutIter->end(), isZero);
      if (nonZeroBegin != lutIter->end()) {
        return {begin.getContainer(), lutIter, nonZeroBegin};
      }

      // go over all remaining buckets
      for (++lutIter; lutIter != end.getLUTIterator(); ++lutIter) {
        // and process each element in each bucket
        auto nonZeroBegin = std::find_if_not(lutIter->begin(), lutIter->end(), isZero);
        if (nonZeroBegin != lutIter->end()) {
          return {begin.getContainer(), lutIter, nonZeroBegin};
        }
      }
    }
    // go over the tail, i.e. the last, possibly incomplete or empty bucket
    if (end.getBucketIterator() != bucket_iterator{}) {
      auto nonZeroBegin = std::find_if_not(lutIter->begin(), end.getBucketIterator(), isZero);
      if (nonZeroBegin != lutIter->end()) {
        return {begin.getContainer(), lutIter, nonZeroBegin};
      }
    }
    return end;
  }();

  //empty
  if (nonZeroBegin == end) {
    return {end, end};
  }

  iterator_type nonZeroEnd = [&]() -> iterator_type {
    auto lutIter = end.getLUTIterator();
    //start at the tail, i.e. the last, possibly incomplete or empty bucket
    if (end.getBucketIterator() != bucket_iterator{}) {
      // if a tail exists, process it
      auto nonZeroEnd = std::find_if_not(std::make_reverse_iterator(end.getBucketIterator()), lutIter->rend(), isZero);
      if (nonZeroEnd != lutIter->rend()) {
        return {begin.getContainer(), lutIter, nonZeroEnd.base()};
      }
    }

    // go over all full buckets, appart from the last one.
    for (; lutIter-- > begin.getLUTIterator();) {
      // and process each element in each bucket
      auto nonZeroEnd = std::find_if_not(lutIter->rbegin(), lutIter->rend(), isZero);
      if (nonZeroEnd != lutIter->rend()) {
        return {begin.getContainer(), lutIter, nonZeroEnd.base()};
      }
    }

    //finish at first ,possibly incomplete bucket
    assert(lutIter == begin.getLUTIterator());
    auto bucketREnd = std::make_reverse_iterator(begin.getBucketIterator());
    auto nonZeroEnd = std::find_if_not(lutIter->rbegin(), bucketREnd, isZero);
    if (nonZeroEnd != bucketREnd) {
      return {begin.getContainer(), lutIter, nonZeroEnd.base()};
    }
    return begin;
  }();

  return {nonZeroBegin, nonZeroEnd};
};

} // namespace o2::rans::internal

#endif /* RANS_INTERNAL_TRANSFORM_SPARSEALGORITHM_H_ */