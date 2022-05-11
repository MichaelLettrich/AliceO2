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

#ifndef INCLUDE_RANS_HASHMAPFREQUENCYCONTAINER_H_
#define INCLUDE_RANS_HASHMAPFREQUENCYCONTAINER_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
#include <type_traits>
#include <vector>

#include <gsl/span>
#include <absl/container/flat_hash_map.h>

#include <fairlogger/Logger.h>

#include "rANS/definitions.h"
#include "rANS/internal/helper.h"
#include "rANS/utils/HistogramView.h"

namespace o2
{
namespace rans
{

template <typename source_T>
class HashMapFrequencyContainer
{
 public:
  using source_type = source_T;
  using index_type = source_type;
  using value_type = count_t;
  using container_type = absl::flat_hash_map<index_type, value_type>;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;
  using const_iterator = typename container_type::const_iterator;

  HashMapFrequencyContainer() = default;
  template <typename freq_IT, std::enable_if_t<internal::isIntegralIter_v<freq_IT>, bool> = true>
  HashMapFrequencyContainer(freq_IT begin, freq_IT end, source_type offset);

  // accessors
  [[nodiscard]] value_type operator[](source_type sourceSymbol) const;

  [[nodiscard]] inline const_pointer data() const noexcept { return mContainer.data(); };

  [[nodiscard]] inline const_iterator cbegin() const noexcept { return mContainer.cbegin(); };

  [[nodiscard]] inline const_iterator cend() const noexcept { return mContainer.cend(); };

  [[nodiscard]] inline const_iterator begin() const noexcept { return cbegin(); };

  [[nodiscard]] inline const_iterator end() const noexcept { return cend(); };

  [[nodiscard]] inline constexpr size_type size() const noexcept { return mContainer.size(); };

  [[nodiscard]] inline bool empty() const noexcept { return mNSamples == 0; };

  [[nodiscard]] inline size_type computeNUsedAlphabetSymbols() const noexcept { return this->size(); };

  [[nodiscard]] const_iterator findMinSymbol() const noexcept;

  [[nodiscard]] const_iterator findMaxSymbol() const noexcept;

  [[nodiscard]] std::pair<const_iterator, const_iterator> findMinMaxSymbol() const noexcept;

  [[nodiscard]] size_type getNumSamples() const noexcept { return mNSamples; };

  // operations
  template <typename source_IT, std::enable_if_t<internal::isCompatibleIter_v<source_type, source_IT>, bool> = true>
  HashMapFrequencyContainer& addSamples(source_IT begin, source_IT end);

  HashMapFrequencyContainer& addSamples(gsl::span<source_type> samples);

  template <typename freq_IT>
  HashMapFrequencyContainer& addFrequencies(freq_IT begin, freq_IT end, source_type offset);

  inline HashMapFrequencyContainer& addFrequencies(gsl::span<value_type> span, source_type offset)
  {
    addFrequencies(span.begin(), span.end(), offset);
  };

  HashMapFrequencyContainer& operator+(HashMapFrequencyContainer& other)
  {
    addFrequencies(other.cbegin(), other.cbegin(), 0);
  };

  container_type release() && noexcept;

 private:
  container_type mContainer{};
  size_type mNSamples{};
};

template <typename source_T>
[[nodiscard]] inline auto HashMapFrequencyContainer<source_T>::operator[](source_type sourceSymbol) const -> value_type
{
  auto iter = mContainer.find(sourceSymbol);
  if (iter == mContainer.end()) {
    return 0;
  } else {
    return iter->second;
  }
};

template <typename source_T>
[[nodiscard]] inline auto HashMapFrequencyContainer<source_T>::findMinSymbol() const noexcept -> const_iterator
{
  if (this->empty()) {
    return this->cend();
  } else {
    return std::min_element(this->cbegin(), this->cend());
  }
}

template <typename source_T>
[[nodiscard]] inline auto HashMapFrequencyContainer<source_T>::findMaxSymbol() const noexcept -> const_iterator
{
  if (this->empty()) {
    return this->cend();
  } else {
    return std::max_element(this->cbegin(), this->cend());
  }
}

template <typename source_T>
[[nodiscard]] auto HashMapFrequencyContainer<source_T>::findMinMaxSymbol() const noexcept -> std::pair<const_iterator, const_iterator>
{
  if (this->empty()) {
    return {this->cend(), this->cend()};
  } else {
    return std::minmax_element(this->begin(), this->end());
  }
}

template <typename source_T>
template <typename source_IT, std::enable_if_t<o2::rans::internal::isCompatibleIter_v<source_T, source_IT>, bool>>
auto HashMapFrequencyContainer<source_T>::addSamples(source_IT begin, source_IT end) -> HashMapFrequencyContainer&
{
  if constexpr (std::is_pointer_v<source_IT>) {
    return addSamples({begin, end});
  } else {
    std::for_each(begin, end, [this](const source_type& symbol) { 
      ++this->mNSamples;
      ++this->mContainer[static_cast<index_type>(symbol)]; });
  }
  return *this;
}

template <typename source_T>
auto HashMapFrequencyContainer<source_T>::addSamples(gsl::span<source_type> samples) -> HashMapFrequencyContainer&
{
  return this->addSamples(samples.begin(), samples.end());
}

template <typename source_T>
template <typename freq_IT>
auto HashMapFrequencyContainer<source_T>::addFrequencies(freq_IT begin, freq_IT end, source_type offset) -> HashMapFrequencyContainer&
{
  source_type sourceSymbol = offset;
  for (auto iter = begin; iter != end; ++iter) {
    auto value = *iter;
    if (value > 0) {
      mContainer[sourceSymbol] += value;
      mNSamples += value;
    }
    ++sourceSymbol;
  }
  return *this;
}

} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_HASHMAPFREQUENCYCONTAINER_H_ */
