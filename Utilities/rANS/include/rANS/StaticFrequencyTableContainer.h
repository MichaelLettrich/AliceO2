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

#ifndef INCLUDE_RANS_STATICFREQUENCYCONTAINER_H_
#define INCLUDE_RANS_STATICFREQUENCYCONTAINER_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
#include <type_traits>
#include <vector>

#include <gsl/span>

#include <fairlogger/Logger.h>

#include "rANS/definitions.h"
#include "rANS/internal/helper.h"
#include "rANS/utils/HistogramView.h"

namespace o2
{
namespace rans
{

namespace internal
{

inline uint64_t load64(const void* __restrict src)
{
  uint64_t ret;
  std::memcpy(&ret, src, 8);
  return ret;
};

template <typename T>
inline size_t itemsPerQWord()
{
  return sizeof(uint64_t) / sizeof(T);
}

} // namespace internal

template <typename source_T>
class StaticFrequencyContainer
{
 public:
  using source_type = source_T;
  using index_type = std::make_unsigned_t<source_type>;
  using value_type = count_t;
  using container_type = std::vector<value_type>;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;
  using const_iterator = typename container_type::const_iterator;

  static_assert(sizeof(index_type) <= 2, "This datatype requires a <=16Bit datatype for source_T");

  StaticFrequencyContainer() : mContainer(this->size(), 0){};
  template <typename freq_IT, std::enable_if_t<internal::isIntegralIter_v<freq_IT>, bool> = true>
  StaticFrequencyContainer(freq_IT begin, freq_IT end, source_type offset);

  // accessors
  [[nodiscard]] inline value_type operator[](source_type sourceSymbol) const { return mContainer[static_cast<index_type>(sourceSymbol)]; };

  [[nodiscard]] inline const_pointer data() const noexcept { return mContainer.data(); };

  [[nodiscard]] inline const_iterator cbegin() const noexcept { return mContainer.begin(); };

  [[nodiscard]] inline const_iterator cend() const noexcept { return mContainer.end(); };

  [[nodiscard]] inline const_iterator begin() const noexcept { return cbegin(); };

  [[nodiscard]] inline const_iterator end() const noexcept { return cend(); };

  [[nodiscard]] inline constexpr size_type size() const noexcept { return internal::pow2(internal::toBits(sizeof(index_type))); };

  [[nodiscard]] inline bool empty() const noexcept { return mNSamples == 0; };

  [[nodiscard]] size_type computeNUsedAlphabetSymbols() const noexcept;

  [[nodiscard]] const_iterator findMinSymbol() const noexcept;

  [[nodiscard]] const_iterator findMaxSymbol() const noexcept;

  [[nodiscard]] std::pair<const_iterator, const_iterator> findMinMaxSymbol() const noexcept;

  [[nodiscard]] size_type getNumSamples() const noexcept { return mNSamples; };

  // operations
  template <typename source_IT, std::enable_if_t<internal::isCompatibleIter_v<source_type, source_IT>, bool> = true>
  StaticFrequencyContainer& addSamples(source_IT begin, source_IT end);

  StaticFrequencyContainer& addSamples(gsl::span<source_type> samples);

  template <typename freq_IT>
  StaticFrequencyContainer& addFrequencies(freq_IT begin, freq_IT end, source_type offset);

  inline StaticFrequencyContainer& addFrequencies(gsl::span<value_type> span, source_type offset)
  {
    return addFrequencies(span.begin(), span.end(), offset);
  };

  StaticFrequencyContainer& operator+(StaticFrequencyContainer& other)
  {
    return addFrequencies(other.cbegin(), other.cbegin(), 0);
  };

  container_type release() && noexcept;

 private:
  container_type mContainer{};
  size_type mNSamples{};
};

template <typename source_T>
[[nodiscard]] inline auto StaticFrequencyContainer<source_T>::computeNUsedAlphabetSymbols() const noexcept -> size_type
{
  return std::count_if(mContainer.begin(), mContainer.end(), [](value_type v) { return v != static_cast<value_type>(0); });
};

template <typename source_T>
[[nodiscard]] inline auto StaticFrequencyContainer<source_T>::findMinSymbol() const noexcept -> const_iterator
{
  if (this->empty()) {
    return this->cend();
  } else {
    return std::find_if(this->cbegin(), this->cend(), [](const_reference frequency) { return frequency != 0; });
  }
}

template <typename source_T>
[[nodiscard]] inline auto StaticFrequencyContainer<source_T>::findMaxSymbol() const noexcept -> const_iterator
{
  if (this->empty()) {
    return this->cend();
  } else {
    return std::find_if(mContainer.crbegin(), mContainer.crend(), [](const_reference frequency) { return frequency != 0; }).base();
  }
}

template <typename source_T>
[[nodiscard]] auto StaticFrequencyContainer<source_T>::findMinMaxSymbol() const noexcept -> std::pair<const_iterator, const_iterator>
{
  if (this->empty()) {
    return {this->cend(), this->cend()};
  } else {
    return {this->findMinSymbol(), this->findMaxSymbol()};
  }
}

template <typename source_T>
template <typename source_IT, std::enable_if_t<o2::rans::internal::isCompatibleIter_v<source_T, source_IT>, bool>>
auto StaticFrequencyContainer<source_T>::addSamples(source_IT begin, source_IT end) -> StaticFrequencyContainer&
{
  if constexpr (std::is_pointer_v<source_IT>) {
    return addSamples({begin, end});
  } else {
    std::for_each(begin, end, [this](const source_type& symbol) { 
      ++this->mNSamples;
      ++this->mContainer[static_cast<const index_type&>(symbol)]; });
  }
  return *this;
}

template <typename source_T>
auto StaticFrequencyContainer<source_T>::addSamples(gsl::span<source_type> samples) -> StaticFrequencyContainer&
{
  const auto begin = samples.data();
  const auto end = begin + samples.size();
  constexpr size_t ElemsPerQWord = sizeof(uint64_t) / sizeof(source_type);
  constexpr size_t nUnroll = 2 * ElemsPerQWord;
  auto iter = begin;

  if constexpr (sizeof(source_type) == 1) {

    alignas(64) value_type histograms[3][256]{0}; //al"ign to cache-line

    auto addQWord = [&, this](uint64_t in64) {
      uint64_t i = in64;
      ++histograms[0][(static_cast<uint8_t>(i))];
      ++histograms[1][(static_cast<uint16_t>(i) >> 8)];
      i >>= 16;
      ++histograms[2][(static_cast<uint8_t>(i))];
      ++this->mContainer[(static_cast<uint16_t>(i) >> 8)];
      i = in64 >>= 32;
      ++histograms[0][(static_cast<uint8_t>(i))];
      ++histograms[1][(static_cast<uint16_t>(i) >> 8)];
      i >>= 16;
      ++histograms[2][(static_cast<uint8_t>(i))];
      ++this->mContainer[(static_cast<uint16_t>(i) >> 8)];
    };

    if (end - nUnroll > begin) {
      for (; iter < end - nUnroll; iter += nUnroll) {
        addQWord(internal::load64(iter));
        addQWord(internal::load64(iter + ElemsPerQWord));
        mNSamples += nUnroll;
        __builtin_prefetch(iter + 512, 0);
      }
    }

    while (iter != end) {
      ++mNSamples;
      ++mContainer[*iter++];
    }

#pragma gcc unroll(3)
    for (size_t j = 0; j < 3; ++j) {
#pragma omp simd
      for (size_t i = 0; i < 256; ++i) {
        mContainer.data()[i] += histograms[j][i];
      }
    }
  } else {
    container_type histogram(this->size(), 0);

    auto addQWord = [&, this](uint64_t in64) {
      uint64_t i = in64;
      ++this->mContainer[static_cast<uint16_t>(i)];
      ++histogram[static_cast<uint32_t>(i) >> 16];
      i = in64 >> 32;
      ++this->mContainer[static_cast<uint16_t>(i)];
      ++histogram[static_cast<uint32_t>(i) >> 16];
    };

    if (end - nUnroll > begin) {
      for (; iter < end - nUnroll; iter += nUnroll) {
        addQWord(internal::load64(iter));
        addQWord(internal::load64(iter + ElemsPerQWord));
        mNSamples += nUnroll;
        __builtin_prefetch(iter + 512, 0);
      }
    }

    while (iter != end) {
      ++mNSamples;
      ++mContainer[*iter++];
    }

#pragma omp simd
    for (size_t i = 0; i < this->size(); ++i) {
      mContainer.data()[i] += histogram.data()[i];
    }
  }

  return *this;
} // namespace rans

template <typename source_T>
template <typename freq_IT>
auto StaticFrequencyContainer<source_T>::addFrequencies(freq_IT begin, freq_IT end, source_type offset) -> StaticFrequencyContainer&
{
  // bounds check
  utils::HistogramView addedHistogram{begin, end, offset};
  utils::HistogramView<container_type::iterator> thisHistogram{};
  if (std::is_signed_v<source_type>) {
    thisHistogram = utils::HistogramView{mContainer.begin(), mContainer.end(), std::numeric_limits<source_type>::min()};
  } else {
    thisHistogram = utils::HistogramView{mContainer.begin(), mContainer.end(), 0};
  }
  const bool invalidBounds = (utils::leftOffset(thisHistogram, addedHistogram) < 0) || (utils::rightOffset(thisHistogram, addedHistogram) > 0);

  if (invalidBounds) {
    throw std::runtime_error(fmt::format("Incompatible Frequency table dimensions: Cannot add [{},{}] to [{}, {}] ",
                                         addedHistogram.getMin(),
                                         addedHistogram.getMax(),
                                         thisHistogram.getMin(),
                                         thisHistogram.getMax()));
  }

  auto iter = addedHistogram.begin();
  for (source_type i = addedHistogram.getOffset(); i < static_cast<source_type>(addedHistogram.getOffset() + addedHistogram.size()); ++i) {
    mNSamples += *iter;
    mContainer[static_cast<index_type>(i)] += *iter;
    ++iter;
  };
  return *this;
}

} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_STATICFREQUENCYCONTAINER_H_ */
