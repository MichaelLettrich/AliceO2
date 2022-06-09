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

#ifndef INCLUDE_RANS_DYNAMICFREQUENCYTABLE_H_
#define INCLUDE_RANS_DYNAMICFREQUENCYTABLE_H_

#include <algorithm>
#include <cassert>

#include <gsl/span>

#include <fairlogger/Logger.h>
#include <immintrin.h>
#include <utility>

#include "rANS/internal/helper.h"
#include "rANS/internal/FrequencyTableBase.h"
#include "rANS/internal/DynamicFrequencyContainer.h"
#include "rANS/utils/HistogramView.h"

#include "rANS/internal/backend/simd/kernel.h"

namespace o2
{
namespace rans
{

namespace internal
{

namespace impl
{
template <typename source_T>
inline std::pair<source_T, source_T> minmaxImpl(const source_T* begin, const source_T* end)
{
  const auto [minIter, maxIter] = std::minmax_element(begin, end);
  return {*minIter, *maxIter};
};

template <>
std::pair<uint32_t, uint32_t> minmaxImpl<uint32_t>(const uint32_t* begin, const uint32_t* end)
{
  constexpr size_t ElemsPerLane = 4;
  constexpr size_t nUnroll = 2 * ElemsPerLane;
  auto iter = begin;

  uint32_t min = *iter;
  uint32_t max = *iter;

  if (end - nUnroll > begin) {

    __m128i minVec[2];
    __m128i maxVec[2];
    __m128i tmpVec[2];

    minVec[0] = _mm_loadu_si128(reinterpret_cast<const __m128i_u*>(iter));
    minVec[1] = minVec[0];
    maxVec[0] = minVec[0];
    maxVec[1] = minVec[0];

    for (; iter < end - nUnroll; iter += nUnroll) {
      tmpVec[0] = _mm_loadu_si128(reinterpret_cast<const __m128i_u*>(iter));
      minVec[0] = _mm_min_epi32(minVec[0], tmpVec[0]);
      maxVec[0] = _mm_max_epi32(maxVec[0], tmpVec[0]);

      tmpVec[1] = _mm_loadu_si128(reinterpret_cast<const __m128i_u*>(iter) + 1);
      minVec[1] = _mm_min_epi32(minVec[1], tmpVec[1]);
      maxVec[1] = _mm_max_epi32(maxVec[1], tmpVec[1]);

      __builtin_prefetch(iter + 512, 0);
    }

    minVec[0] = _mm_min_epu32(minVec[0], minVec[1]);
    maxVec[0] = _mm_max_epu32(maxVec[0], maxVec[1]);

    uint32_t tmpMin[ElemsPerLane];
    uint32_t tmpMax[ElemsPerLane];
    _mm_storeu_si128(reinterpret_cast<__m128i*>(tmpMin), minVec[0]);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(tmpMax), maxVec[0]);

    for (size_t i = 0; i < ElemsPerLane; ++i) {
      min = std::min(tmpMin[i], min);
      max = std::max(tmpMax[i], max);
    }
  }

  while (iter != end) {
    min = std::min(*iter, min);
    max = std::max(*iter, max);
    ++iter;
  }

  return {min, max};
};

template <>
std::pair<int32_t, int32_t> minmaxImpl<int32_t>(const int32_t* begin, const int32_t* end)
{
  constexpr size_t ElemsPerLane = 4;
  constexpr size_t nUnroll = 2 * ElemsPerLane;
  auto iter = begin;

  int32_t min = *iter;
  int32_t max = *iter;

  if (end - nUnroll > begin) {

    __m128i minVec[2];
    __m128i maxVec[2];
    __m128i tmpVec[2];

    minVec[0] = _mm_loadu_si128(reinterpret_cast<const __m128i_u*>(iter));
    minVec[1] = minVec[0];
    maxVec[0] = minVec[0];
    maxVec[1] = minVec[0];

    for (; iter < end - nUnroll; iter += nUnroll) {
      tmpVec[0] = _mm_loadu_si128(reinterpret_cast<const __m128i_u*>(iter));
      minVec[0] = _mm_min_epu32(minVec[0], tmpVec[0]);
      maxVec[0] = _mm_max_epu32(maxVec[0], tmpVec[0]);

      tmpVec[1] = _mm_loadu_si128(reinterpret_cast<const __m128i_u*>(iter) + 1);
      minVec[1] = _mm_min_epu32(minVec[1], tmpVec[1]);
      maxVec[1] = _mm_max_epu32(maxVec[1], tmpVec[1]);

      __builtin_prefetch(iter + 512, 0);
    }

    minVec[0] = _mm_min_epu32(minVec[0], minVec[1]);
    maxVec[0] = _mm_max_epu32(maxVec[0], maxVec[1]);

    int32_t tmpMin[ElemsPerLane];
    int32_t tmpMax[ElemsPerLane];
    _mm_storeu_si128(reinterpret_cast<__m128i*>(tmpMin), minVec[0]);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(tmpMax), maxVec[0]);

    for (size_t i = 0; i < ElemsPerLane; ++i) {
      min = std::min(tmpMin[i], min);
      max = std::max(tmpMax[i], max);
    }
  }

  while (iter != end) {
    min = std::min(*iter, min);
    max = std::max(*iter, max);
    ++iter;
  }

  return {min, max};
};
}; // namespace impl

template <typename source_T>
inline std::pair<source_T, source_T> minmax(gsl::span<const source_T> range)
{
  const auto begin = range.data();
  const auto end = begin + range.size();
  return impl::minmaxImpl<source_T>(begin, end);
};

}; // namespace internal

template <typename source_T>
class DynamicFrequencyTable : public internal::DynamicFrequencyContainer<source_T>,
                              public FrequencyTableBase<source_T,
                                                        typename internal::DynamicFrequencyContainer<source_T>::value_type,
                                                        DynamicFrequencyTable<source_T>>
{
  using containerBase_type = internal::DynamicFrequencyContainer<source_T>;
  using frequencyTableBase_type = FrequencyTableBase<source_T, typename internal::DynamicFrequencyContainer<source_T>::value_type, DynamicFrequencyTable<source_T>>;

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

  DynamicFrequencyTable() = default;

  template <typename freq_IT>
  DynamicFrequencyTable(freq_IT begin, freq_IT end, source_type offset) : containerBase_type(), frequencyTableBase_type{begin, end, offset} {};

  DynamicFrequencyTable& addSamples(gsl::span<const source_type> span);

  template <typename source_IT>
  DynamicFrequencyTable& addSamples(source_IT begin, source_IT end);

  template <typename source_IT>
  DynamicFrequencyTable& addSamples(source_IT begin, source_IT end, source_type min, source_type max);

  DynamicFrequencyTable& addSamples(gsl::span<const source_type> span, source_type min, source_type max);

  // operations

  using frequencyTableBase_type::addFrequencies;

  template <typename freq_IT>
  DynamicFrequencyTable& addFrequencies(freq_IT begin, freq_IT end, source_type offset);

  DynamicFrequencyTable& resize(source_type min, source_type max);

  inline DynamicFrequencyTable& resize(size_type newSize)
  {
    return resize(this->getMinSymbol(), this->getMinSymbol() + newSize);
  };

  friend void swap(DynamicFrequencyTable& a, DynamicFrequencyTable& b) noexcept
  {
    using std::swap;
    swap(static_cast<typename DynamicFrequencyTable::containerBase_type&>(a),
         static_cast<typename DynamicFrequencyTable::containerBase_type&>(b));
  };
};
template <typename source_T>
inline auto DynamicFrequencyTable<source_T>::addSamples(gsl::span<const source_type> samples) -> DynamicFrequencyTable&
{
  if (samples.size() > 0) {
    const auto [min, max] = internal::minmax(samples);
    addSamples(samples, min, max);
  } else {
    LOG(warning) << "Passed empty message to " << __func__; // RS this is ok for empty columns
  }
  return *this;
}

template <typename source_T>
template <typename source_IT>
inline auto DynamicFrequencyTable<source_T>::addSamples(source_IT begin, source_IT end) -> DynamicFrequencyTable&
{
  if (begin != end) {
    const auto [minIter, maxIter] = std::minmax_element(begin, end);
    addSamples(begin, end, *minIter, *maxIter);
  } else {
    LOG(warning) << "Passed empty message to " << __func__; // RS this is ok for empty columns
  }
  return *this;
}

template <typename source_T>
inline auto DynamicFrequencyTable<source_T>::addSamples(gsl::span<const source_type> samples, source_type min, source_type max) -> DynamicFrequencyTable&
{
  if (samples.empty()) {
    return *this;
  }

  if constexpr (sizeof(source_T) == 4) {
    this->resize(min, max);

    const auto begin = samples.data();
    const auto end = begin + samples.size();
    constexpr size_t ElemsPerQWord = sizeof(uint64_t) / sizeof(source_type);
    constexpr size_t nUnroll = 4 * ElemsPerQWord;
    auto iter = begin;

    const source_type offset = this->getOffset();

    auto addQWord = [&, this](uint64_t in64) {
      uint64_t i = in64;
      ++this->mContainer[static_cast<source_type>(i) - offset];
      i = in64 >> 32;
      ++this->mContainer[static_cast<source_type>(i) - offset];
    };

    if (end - nUnroll > begin) {
      for (; iter < end - nUnroll; iter += nUnroll) {
        addQWord(internal::load64(iter));
        addQWord(internal::load64(iter + ElemsPerQWord));
        addQWord(internal::load64(iter + 2 * ElemsPerQWord));
        addQWord(internal::load64(iter + 3 * ElemsPerQWord));
        this->mNSamples += nUnroll;
        __builtin_prefetch(iter + 512, 0);
      }
    }

    while (iter != end) {
      ++this->mNSamples;
      ++this->mContainer[(*iter++) - offset];
    }
  } else {
    return addSamples(samples.data(), samples.data() + samples.size(), min, max);
  }
  return *this;
};

template <typename source_T>
template <typename source_IT>
auto DynamicFrequencyTable<source_T>::addSamples(source_IT begin, source_IT end, source_type min, source_type max) -> DynamicFrequencyTable&
{
  if (begin == end) {
    LOG(warning) << "Passed empty message to " << __func__; // RS this is ok for empty columns
  } else {
    this->resize(min, max);
    // add new symbols
    std::for_each(begin, end, [this](source_type symbol) { ++this->getSymbol(symbol);
          ++this->mNSamples; });
  }
  return *this;
}

template <typename source_T>
template <typename freq_IT>
auto DynamicFrequencyTable<source_T>::addFrequencies(freq_IT begin, freq_IT end, source_type offset) -> DynamicFrequencyTable&
{

  auto frequencyCountingDecorator = [this](value_type frequency) {
    this->mNSamples += frequency;
    return frequency;
  };

  auto thisHistogram = utils::HistogramView{this->mContainer.begin(), this->mContainer.end(), this->mOffset};
  auto addedHistogram = utils::trim(utils::HistogramView{begin, end, offset});
  if (addedHistogram.empty()) {
    LOG(warning) << "Passed empty FrequencyTable to " << __func__; // RS this is ok for empty columns
  } else {

    const symbol_t newMin = std::min(thisHistogram.getMin(), addedHistogram.getMin());
    const symbol_t newMax = std::max(thisHistogram.getMax(), addedHistogram.getMax());

    if (thisHistogram.empty()) {
      this->mContainer = histogram_t(addedHistogram.size());
      std::transform(addedHistogram.begin(), addedHistogram.end(), this->mContainer.begin(), [this, frequencyCountingDecorator](count_t frequency) {
        return frequencyCountingDecorator(frequency);
      });
      this->mOffset = addedHistogram.getOffset();
    } else {
      const symbol_t newSize = newMax - newMin + 1;
      histogram_t newFreequencyTable(newSize, 0);
      auto newHistogram = utils::HistogramView{newFreequencyTable.begin(), newFreequencyTable.end(), newMin};
      auto histogramOverlap = utils::intersection(newHistogram, thisHistogram);
      assert(!histogramOverlap.empty());
      assert(histogramOverlap.size() == thisHistogram.size());
      std::copy(thisHistogram.begin(), thisHistogram.end(), histogramOverlap.begin());

      histogramOverlap = utils::intersection(newHistogram, addedHistogram);
      assert(!histogramOverlap.empty());
      assert(histogramOverlap.size() == addedHistogram.size());
      std::transform(addedHistogram.begin(), addedHistogram.end(),
                     histogramOverlap.begin(), histogramOverlap.begin(),
                     [this, frequencyCountingDecorator](const count_t& a, const count_t& b) { return frequencyCountingDecorator(a) + b; });

      this->mContainer = std::move(newFreequencyTable);
      this->mOffset = newHistogram.getOffset();
    }
  }

  return *this;
}

template <typename source_T>
auto DynamicFrequencyTable<source_T>::resize(source_type min, source_type max) -> DynamicFrequencyTable&
{

  auto getMaxSymbol = [this]() {
    return static_cast<source_type>(this->getOffset() + std::max(0l, static_cast<int32_t>(this->size()) - 1l));
  };

  min = std::min(min, this->getOffset());
  max = std::max(max, getMaxSymbol());

  if (min > max) {
    throw std::runtime_error(fmt::format("{} failed: min {} > max {} ", __func__, min, max));
  }

  const size_type newSize = max - min + 1;
  const source_type oldOffset = this->mOffset;
  this->mOffset = min;
  this->mNSamples = 0;

  if (this->mContainer.empty()) {
    this->mContainer.resize(newSize, 0);
    return *this;
  } else {
    container_type oldFrequencyTable = std::move(this->mContainer);
    auto oldHistogram = utils::HistogramView{oldFrequencyTable.begin(), oldFrequencyTable.end(), oldOffset};
    this->mContainer = histogram_t(newSize, 0);
    return this->addFrequencies(oldHistogram.begin(), oldHistogram.end(), oldHistogram.getMin());
  }
}

} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_DYNAMICFREQUENCYTABLE_H_ */
