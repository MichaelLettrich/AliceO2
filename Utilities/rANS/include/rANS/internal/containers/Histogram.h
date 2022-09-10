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

/// @file   Histogram.h
/// @author Michael Lettrich
/// @brief  Histogram for source symbols used to estimate symbol probabilities for entropy coding

#ifndef RANS_INTERNAL_CONTAINERS_HISTOGRAM_H_
#define RANS_INTERNAL_CONTAINERS_HISTOGRAM_H_

#include <algorithm>
#include <cassert>

#include <gsl/span>

#include <fairlogger/Logger.h>
#include <immintrin.h>
#include <utility>

#include "rANS/internal/common/utils.h"
#include "rANS/internal/containers/HistogramInterface.h"
#include "rANS/internal/containers/CountingContainer.h"
#include "rANS/internal/containers/HistogramView.h"

#include "rANS/internal/common/simdops.h"
#include "rANS/internal/common/simdtypes.h"

namespace o2::rans
{

namespace internal
{
namespace histogramImpl
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
}; // namespace histogramImpl

template <typename source_T>
inline std::pair<source_T, source_T> minmax(gsl::span<const source_T> range)
{
  const auto begin = range.data();
  const auto end = begin + range.size();
  return histogramImpl::minmaxImpl<source_T>(begin, end);
};

} // namespace internal

template <typename source_T, typename = void>
class Histogram;

template <typename source_T>
class Histogram<source_T, std::enable_if_t<sizeof(source_T) == 4>> : public internal::CountingContainer<source_T>,
                                                                     public internal::HistogramInterface<source_T,
                                                                                                         typename internal::CountingContainer<source_T>::value_type,
                                                                                                         Histogram<source_T>>
{
  using containerBase_type = internal::CountingContainer<source_T>;
  using HistogramInterface_type = internal::HistogramInterface<source_T, typename internal::CountingContainer<source_T>::value_type, Histogram<source_T>>;

 public:
  using source_type = source_T;
  using value_type = typename containerBase_type::value_type;
  using container_type = typename containerBase_type::container_type;
  using size_type = typename containerBase_type::size_type;
  using difference_type = typename containerBase_type::difference_type;
  using reference = typename containerBase_type::reference;
  using const_reference = typename containerBase_type::const_reference;
  using pointer = typename containerBase_type::pointer;
  using const_pointer = typename containerBase_type::const_pointer;
  using const_iterator = typename containerBase_type::const_iterator;
  using iterator = typename containerBase_type::iterator;
  using const_reverse_iterator = typename containerBase_type::const_reverse_iterator;
  using reverse_iterator = typename containerBase_type::reverse_iterator;

  Histogram() = default;

  template <typename freq_IT>
  Histogram(freq_IT begin, freq_IT end, source_type offset) : containerBase_type(), HistogramInterface_type{begin, end, offset} {};

  Histogram& addSamples(gsl::span<const source_type> span);

  template <typename source_IT>
  Histogram& addSamples(source_IT begin, source_IT end);

  template <typename source_IT>
  Histogram& addSamples(source_IT begin, source_IT end, source_type min, source_type max);

  Histogram& addSamples(gsl::span<const source_type> span, source_type min, source_type max);

  // operations

  using HistogramInterface_type::addFrequencies;

  template <typename freq_IT>
  Histogram& addFrequencies(freq_IT begin, freq_IT end, source_type offset);

  Histogram& resize(source_type min, source_type max);

  inline Histogram& resize(size_type newSize)
  {
    return resize(this->getMinSymbol(), this->getMinSymbol() + newSize);
  };

  friend void swap(Histogram& a, Histogram& b) noexcept
  {
    using std::swap;
    swap(static_cast<typename Histogram::containerBase_type&>(a),
         static_cast<typename Histogram::containerBase_type&>(b));
  };
};
template <typename source_T>
inline auto Histogram<source_T, std::enable_if_t<sizeof(source_T) == 4>>::addSamples(gsl::span<const source_type> samples) -> Histogram&
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
inline auto Histogram<source_T, std::enable_if_t<sizeof(source_T) == 4>>::addSamples(source_IT begin, source_IT end) -> Histogram&
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
inline auto Histogram<source_T, std::enable_if_t<sizeof(source_T) == 4>>::addSamples(gsl::span<const source_type> samples, source_type min, source_type max) -> Histogram&
{
  using namespace internal;

  if (samples.empty()) {
    return *this;
  }
  this->resize(min, max);

  const auto begin = samples.data();
  const auto end = begin + samples.size();
  constexpr size_t ElemsPerQWord = sizeof(uint64_t) / sizeof(source_type);
  constexpr size_t nUnroll = 4 * ElemsPerQWord;
  auto iter = begin;

  const source_type offset = this->getOffset();

  auto addQWord = [&, this](uint64_t in64) {
    uint64_t i = in64;
    ++this->mContainer[static_cast<source_type>(i)];
    i = in64 >> 32;
    ++this->mContainer[static_cast<source_type>(i)];
  };

  if (end - nUnroll > begin) {
    for (; iter < end - nUnroll; iter += nUnroll) {
      addQWord(load64(iter));
      addQWord(load64(iter + ElemsPerQWord));
      addQWord(load64(iter + 2 * ElemsPerQWord));
      addQWord(load64(iter + 3 * ElemsPerQWord));
      this->mNSamples += nUnroll;
      __builtin_prefetch(iter + 512, 0);
    }
  }

  while (iter != end) {
    ++this->mNSamples;
    ++this->mContainer[*iter++];
  }
  return *this;
};

template <typename source_T>
template <typename source_IT>
auto Histogram<source_T, std::enable_if_t<sizeof(source_T) == 4>>::addSamples(source_IT begin, source_IT end, source_type min, source_type max) -> Histogram&
{
  if (begin == end) {
    LOG(warning) << "Passed empty message to " << __func__; // RS this is ok for empty columns
  } else {
    this->resize(min, max);
    // add new symbols
    std::for_each(begin, end, [this](source_type symbol) { 
      ++this->mContainer[symbol];
          ++this->mNSamples; });
  }
  return *this;
}

template <typename source_T>
template <typename freq_IT>
auto Histogram<source_T, std::enable_if_t<sizeof(source_T) == 4>>::addFrequencies(freq_IT begin, freq_IT end, source_type offset) -> Histogram&
{

  using namespace internal;

  auto frequencyCountingDecorator = [this](value_type frequency) {
    this->mNSamples += frequency;
    return frequency;
  };

  auto thisHistogramView = HistogramView{this->mContainer.begin(), this->mContainer.end(), this->mContainer.getOffset()};
  auto addedHistogramView = trim(HistogramView{begin, end, offset});
  if (addedHistogramView.empty()) {
    LOG(warning) << "Passed empty Histogram to " << __func__; // RS this is ok for empty columns
  } else {

    const source_type newMin = std::min(thisHistogramView.getMin(), addedHistogramView.getMin());
    const source_type newMax = std::max(thisHistogramView.getMax(), addedHistogramView.getMax());

    if (thisHistogramView.empty()) {
      this->mContainer = container_type(addedHistogramView.size(), addedHistogramView.getOffset());
      std::transform(addedHistogramView.begin(), addedHistogramView.end(), this->mContainer.begin(), [this, frequencyCountingDecorator](count_t frequency) {
        return frequencyCountingDecorator(frequency);
      });
    } else {
      const difference_type newSize = newMax - newMin + 1;
      typename container_type::container_type newHistogram(newSize, 0);
      auto newHistogramView = HistogramView{newHistogram.begin(), newHistogram.end(), newMin};
      auto histogramOverlap = getIntersection(newHistogramView, thisHistogramView);
      assert(!histogramOverlap.empty());
      assert(histogramOverlap.size() == thisHistogramView.size());
      std::copy(thisHistogramView.begin(), thisHistogramView.end(), histogramOverlap.begin());

      histogramOverlap = getIntersection(newHistogramView, addedHistogramView);
      assert(!histogramOverlap.empty());
      assert(histogramOverlap.size() == addedHistogramView.size());
      std::transform(addedHistogramView.begin(), addedHistogramView.end(),
                     histogramOverlap.begin(), histogramOverlap.begin(),
                     [this, frequencyCountingDecorator](const count_t& a, const count_t& b) { return frequencyCountingDecorator(a) + b; });

      this->mContainer = container_type{std::move(newHistogram), static_cast<source_type>(newHistogramView.getOffset())};
    }
  }

  return *this;
}

template <typename source_T>
auto Histogram<source_T, std::enable_if_t<sizeof(source_T) == 4>>::resize(source_type min, source_type max) -> Histogram&
{
  using namespace internal;

  auto getMaxSymbol = [this]() {
    return static_cast<source_type>(this->getOffset() + std::max(0l, static_cast<int32_t>(this->size()) - 1l));
  };

  min = std::min(min, this->getOffset());
  max = std::max(max, getMaxSymbol());

  if (min > max) {
    throw std::runtime_error(fmt::format("{} failed: min {} > max {} ", __func__, min, max));
  }

  const size_type newSize = max - min + 1;
  const source_type oldOffset = this->getOffset();
  this->mNSamples = 0;

  if (this->mContainer.empty()) {
    this->mContainer = container_type{newSize, min};
    return *this;
  } else {
    container_type oldHistogram = std::move(this->mContainer);
    auto oldHistogramView = HistogramView{oldHistogram.begin(), oldHistogram.end(), oldOffset};
    this->mContainer = container_type{newSize, min};
    return this->addFrequencies(oldHistogramView.begin(), oldHistogramView.end(), oldHistogramView.getMin());
  }
}

template <typename source_T>
class Histogram<source_T, std::enable_if_t<sizeof(source_T) <= 2>> : public internal::CountingContainer<source_T>,
                                                                     public internal::HistogramInterface<source_T,
                                                                                                         typename internal::CountingContainer<source_T>::value_type,
                                                                                                         Histogram<source_T>>
{
  using containerBase_type = internal::CountingContainer<source_T>;
  using HistogramInterface_type = internal::HistogramInterface<source_T, typename internal::CountingContainer<source_T>::value_type, Histogram<source_T>>;

 public:
  using source_type = source_T;
  using value_type = typename containerBase_type::value_type;
  using container_type = typename containerBase_type::container_type;
  using size_type = typename containerBase_type::size_type;
  using difference_type = typename containerBase_type::difference_type;
  using reference = typename containerBase_type::reference;
  using const_reference = typename containerBase_type::const_reference;
  using pointer = typename containerBase_type::pointer;
  using const_pointer = typename containerBase_type::const_pointer;
  using const_iterator = typename containerBase_type::const_iterator;
  using iterator = typename containerBase_type::iterator;
  using const_reverse_iterator = typename containerBase_type::const_reverse_iterator;
  using reverse_iterator = typename containerBase_type::reverse_iterator;

  Histogram() = default;

  template <typename freq_IT>
  Histogram(freq_IT begin, freq_IT end, source_type offset) : containerBase_type(), HistogramInterface_type{begin, end, offset} {};

  // operations
  template <typename source_IT>
  Histogram& addSamples(source_IT begin, source_IT end);

  Histogram& addSamples(gsl::span<const source_type> samples);

  template <typename freq_IT>
  Histogram& addFrequencies(freq_IT begin, freq_IT end, source_type offset);

  using HistogramInterface_type::addFrequencies;

  friend void swap(Histogram& a, Histogram& b) noexcept
  {
    using std::swap;
    swap(static_cast<typename Histogram::containerBase_type&>(a),
         static_cast<typename Histogram::containerBase_type&>(b));
  };
};

template <typename source_T>
template <typename source_IT>
auto Histogram<source_T, std::enable_if_t<sizeof(source_T) <= 2>>::addSamples(source_IT begin, source_IT end) -> Histogram&
{
  if constexpr (std::is_pointer_v<source_IT>) {
    return addSamples({begin, end});
  } else {
    std::for_each(begin, end, [this](const source_type& symbol) { 
      ++this->mNSamples;
       ++this->mContainer[symbol]; });
  }
  return *this;
}

template <typename source_T>
auto Histogram<source_T, std::enable_if_t<sizeof(source_T) <= 2>>::addSamples(gsl::span<const source_type> samples) -> Histogram&
{

  using namespace internal;

  if (samples.empty()) {
    return *this;
  }

  const auto begin = samples.data();
  const auto end = begin + samples.size();
  constexpr size_t ElemsPerQWord = sizeof(uint64_t) / sizeof(source_type);
  constexpr size_t nUnroll = 2 * ElemsPerQWord;
  auto iter = begin;

  if constexpr (sizeof(source_type) == 1) {

    std::array<ShiftableVector<source_type, value_type>, 3> histograms{
      {{this->mContainer.size(), this->mContainer.getOffset()},
       {this->mContainer.size(), this->mContainer.getOffset()},
       {this->mContainer.size(), this->mContainer.getOffset()}}};

    auto addQWord = [&, this](uint64_t in64) {
      uint64_t i = in64;
      ++histograms[0][static_cast<source_type>(i)];
      ++histograms[1][static_cast<source_type>(static_cast<uint16_t>(i) >> 8)];
      i >>= 16;
      ++histograms[2][static_cast<source_type>(i)];
      ++this->mContainer[static_cast<source_type>(static_cast<uint16_t>(i) >> 8)];
      i = in64 >>= 32;
      ++histograms[0][static_cast<source_type>(i)];
      ++histograms[1][static_cast<source_type>(static_cast<uint16_t>(i) >> 8)];
      i >>= 16;
      ++histograms[2][static_cast<source_type>(i)];
      ++this->mContainer[static_cast<source_type>(static_cast<uint16_t>(i) >> 8)];
    };

    if (end - nUnroll > begin) {
      for (; iter < end - nUnroll; iter += nUnroll) {
        addQWord(load64(iter));
        addQWord(load64(iter + ElemsPerQWord));
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
        this->mContainer(i) += histograms[j](i);
      }
    }
  } else {
    container_type histogram{this->mContainer.size(), this->mContainer.getOffset()};

    auto addQWord = [&, this](uint64_t in64) {
      uint64_t i = in64;
      ++histogram[static_cast<source_type>(i)];
      ++this->mContainer[static_cast<source_type>(static_cast<uint32_t>(i) >> 16)];
      i = in64 >> 32;
      ++histogram[static_cast<source_type>(i)];
      ++this->mContainer[static_cast<source_type>(static_cast<uint32_t>(i) >> 16)];
    };

    if (end - nUnroll > begin) {
      for (; iter < end - nUnroll; iter += nUnroll) {
        addQWord(load64(iter));
        addQWord(load64(iter + ElemsPerQWord));
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
}

template <typename source_T>
template <typename freq_IT>
auto Histogram<source_T, std::enable_if_t<sizeof(source_T) <= 2>>::addFrequencies(freq_IT begin, freq_IT end, source_type offset) -> Histogram&
{
  using namespace internal;

  // bounds check
  HistogramView addedHistogramView{begin, end, offset};
  HistogramView<typename container_type::iterator> thisHistogramView{};
  thisHistogramView = HistogramView{this->mContainer.begin(), this->mContainer.end(), this->mContainer.getOffset()};
  const bool invalidBounds = (getLeftOffset(thisHistogramView, addedHistogramView) < 0) || (getRightOffset(thisHistogramView, addedHistogramView) > 0);

  if (invalidBounds) {
    throw std::runtime_error(fmt::format("Incompatible Frequency table dimensions: Cannot add [{},{}] to [{}, {}] ",
                                         addedHistogramView.getMin(),
                                         addedHistogramView.getMax(),
                                         thisHistogramView.getMin(),
                                         thisHistogramView.getMax()));
  }

  auto iter = addedHistogramView.begin();
  for (source_type i = addedHistogramView.getOffset(); i < static_cast<source_type>(addedHistogramView.getOffset() + addedHistogramView.size()); ++i) {
    this->mNSamples += *iter;
    this->mContainer[i] += *iter;
    ++iter;
  };
  return *this;
}

} // namespace o2::rans

#endif /* RANS_INTERNAL_CONTAINERS_HISTOGRAM_H_ */
