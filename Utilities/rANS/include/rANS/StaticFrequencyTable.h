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

#ifndef INCLUDE_RANS_STATICFREQUENCYTABLE_H_
#define INCLUDE_RANS_STATICFREQUENCYTABLE_H_

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
#include "rANS/internal/FrequencyTableBase.h"
#include "rANS/internal/StaticFrequencyContainer.h"
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
class StaticFrequencyTable : public internal::StaticFrequencyContainer<source_T>,
                             public FrequencyTableBase<source_T,
                                                       typename internal::StaticFrequencyContainer<source_T>::value_type,
                                                       StaticFrequencyTable<source_T>>
{
  using containerBase_type = internal::StaticFrequencyContainer<source_T>;
  using frequencyTableBase_type = FrequencyTableBase<source_T, typename internal::StaticFrequencyContainer<source_T>::value_type, StaticFrequencyTable<source_T>>;

 public:
  using source_type = source_T;
  using index_type = typename containerBase_type::index_type;
  using value_type = typename containerBase_type::value_type;
  using container_type = typename containerBase_type::container_type;
  using size_type = typename containerBase_type::size_type;
  using difference_type = typename containerBase_type::difference_type;
  using reference = typename containerBase_type::reference;
  using const_reference = typename containerBase_type::const_reference;
  using pointer = typename containerBase_type::pointer;
  using const_pointer = typename containerBase_type::const_pointer;
  using const_iterator = typename containerBase_type::const_iterator;

  static_assert(sizeof(index_type) <= 2, "This datatype requires a <=16Bit datatype for source_T");

  StaticFrequencyTable() = default;

  template <typename freq_IT>
  StaticFrequencyTable(freq_IT begin, freq_IT end, source_type offset) : containerBase_type(), frequencyTableBase_type{begin, end, offset} {};

  // operations
  template <typename source_IT>
  StaticFrequencyTable& addSamples(source_IT begin, source_IT end);

  StaticFrequencyTable& addSamples(gsl::span<const source_type> samples);

  template <typename freq_IT>
  StaticFrequencyTable& addFrequencies(freq_IT begin, freq_IT end, source_type offset);

  using frequencyTableBase_type::addFrequencies;

  friend void swap(StaticFrequencyTable& a, StaticFrequencyTable& b) noexcept
  {
    using std::swap;
    swap(static_cast<typename StaticFrequencyTable::containerBase_type&>(a),
         static_cast<typename StaticFrequencyTable::containerBase_type&>(b));
  };
};

template <typename source_T>
template <typename source_IT>
auto StaticFrequencyTable<source_T>::addSamples(source_IT begin, source_IT end) -> StaticFrequencyTable&
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
auto StaticFrequencyTable<source_T>::addSamples(gsl::span<const source_type> samples) -> StaticFrequencyTable&
{

  if (samples.empty()) {
    return *this;
  }

  const auto begin = samples.data();
  const auto end = begin + samples.size();
  constexpr size_t ElemsPerQWord = sizeof(uint64_t) / sizeof(source_type);
  constexpr size_t nUnroll = 2 * ElemsPerQWord;
  auto iter = begin;

  if constexpr (sizeof(source_type) == 1) {

    alignas(64) value_type histograms[3][256]{0}; // al"ign to cache-line

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
        this->mNSamples += nUnroll;
        __builtin_prefetch(iter + 512, 0);
      }
    }

    while (iter != end) {
      ++this->mNSamples;
      ++this->mContainer[*iter++];
    }

#pragma gcc unroll(3)
    for (size_t j = 0; j < 3; ++j) {
#pragma omp simd
      for (size_t i = 0; i < 256; ++i) {
        this->mContainer.data()[i] += histograms[j][i];
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
        this->mNSamples += nUnroll;
        __builtin_prefetch(iter + 512, 0);
      }
    }

    while (iter != end) {
      ++this->mNSamples;
      ++this->mContainer[*iter++];
    }

#pragma omp simd
    for (size_t i = 0; i < this->size(); ++i) {
      this->mContainer.data()[i] += histogram.data()[i];
    }
  }

  return *this;
} // namespace rans

template <typename source_T>
template <typename freq_IT>
auto StaticFrequencyTable<source_T>::addFrequencies(freq_IT begin, freq_IT end, source_type offset) -> StaticFrequencyTable&
{
  // bounds check
  utils::HistogramView addedHistogram{begin, end, offset};
  utils::HistogramView<typename container_type::iterator> thisHistogram{};
  if (std::is_signed_v<source_type>) {
    thisHistogram = utils::HistogramView{this->mContainer.begin(), this->mContainer.end(), std::numeric_limits<source_type>::min()};
  } else {
    thisHistogram = utils::HistogramView{this->mContainer.begin(), this->mContainer.end(), 0};
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
    this->mNSamples += *iter;
    this->mContainer[static_cast<index_type>(i)] += *iter;
    ++iter;
  };
  return *this;
}

} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_STATICFREQUENCYTABLE_H_ */
