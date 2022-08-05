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

#ifndef INCLUDE_RANS_FREQUENCYTABLE_H_
#define INCLUDE_RANS_FREQUENCYTABLE_H_

#include <algorithm>
#include <cassert>

#include <gsl/span>

#include <fairlogger/Logger.h>
#include <immintrin.h>
#include <utility>

#include "rANS/internal/helper.h"
#include "rANS/internal/FrequencyTableBase.h"
#include "rANS/internal/FrequencyContainer.h"
#include "rANS/utils/HistogramView.h"

#include "rANS/internal/backend/simd/kernel.h"

namespace o2
{
namespace rans
{

template <typename source_T, typename = void>
class FrequencyTable;

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
class FrequencyTable<source_T, std::enable_if_t<sizeof(source_T) == 4>> : public internal::FrequencyContainer<source_T>,
                                                                          public FrequencyTableBase<source_T,
                                                                                                    typename internal::FrequencyContainer<source_T>::value_type,
                                                                                                    FrequencyTable<source_T>>
{
  using containerBase_type = internal::FrequencyContainer<source_T>;
  using frequencyTableBase_type = FrequencyTableBase<source_T, typename internal::FrequencyContainer<source_T>::value_type, FrequencyTable<source_T>>;

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

  FrequencyTable() = default;

  template <typename freq_IT>
  FrequencyTable(freq_IT begin, freq_IT end, source_type offset) : containerBase_type(), frequencyTableBase_type{begin, end, offset} {};

  FrequencyTable& addSamples(gsl::span<const source_type> span);

  template <typename source_IT>
  FrequencyTable& addSamples(source_IT begin, source_IT end);

  template <typename source_IT>
  FrequencyTable& addSamples(source_IT begin, source_IT end, source_type min, source_type max);

  FrequencyTable& addSamples(gsl::span<const source_type> span, source_type min, source_type max);

  // operations

  using frequencyTableBase_type::addFrequencies;

  template <typename freq_IT>
  FrequencyTable& addFrequencies(freq_IT begin, freq_IT end, source_type offset);

  FrequencyTable& resize(source_type min, source_type max);

  inline FrequencyTable& resize(size_type newSize)
  {
    return resize(this->getMinSymbol(), this->getMinSymbol() + newSize);
  };

  friend void swap(FrequencyTable& a, FrequencyTable& b) noexcept
  {
    using std::swap;
    swap(static_cast<typename FrequencyTable::containerBase_type&>(a),
         static_cast<typename FrequencyTable::containerBase_type&>(b));
  };
};
template <typename source_T>
inline auto FrequencyTable<source_T, std::enable_if_t<sizeof(source_T) == 4>>::addSamples(gsl::span<const source_type> samples) -> FrequencyTable&
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
inline auto FrequencyTable<source_T, std::enable_if_t<sizeof(source_T) == 4>>::addSamples(source_IT begin, source_IT end) -> FrequencyTable&
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
inline auto FrequencyTable<source_T, std::enable_if_t<sizeof(source_T) == 4>>::addSamples(gsl::span<const source_type> samples, source_type min, source_type max) -> FrequencyTable&
{
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
    ++this->mContainer[*iter++];
  }
  return *this;
};

template <typename source_T>
template <typename source_IT>
auto FrequencyTable<source_T, std::enable_if_t<sizeof(source_T) == 4>>::addSamples(source_IT begin, source_IT end, source_type min, source_type max) -> FrequencyTable&
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
auto FrequencyTable<source_T, std::enable_if_t<sizeof(source_T) == 4>>::addFrequencies(freq_IT begin, freq_IT end, source_type offset) -> FrequencyTable&
{

  auto frequencyCountingDecorator = [this](value_type frequency) {
    this->mNSamples += frequency;
    return frequency;
  };

  auto thisHistogram = utils::HistogramView{this->mContainer.begin(), this->mContainer.end(), this->mContainer.getOffset()};
  auto addedHistogram = utils::trim(utils::HistogramView{begin, end, offset});
  if (addedHistogram.empty()) {
    LOG(warning) << "Passed empty FrequencyTable to " << __func__; // RS this is ok for empty columns
  } else {

    const symbol_t newMin = std::min(thisHistogram.getMin(), addedHistogram.getMin());
    const symbol_t newMax = std::max(thisHistogram.getMax(), addedHistogram.getMax());

    if (thisHistogram.empty()) {
      this->mContainer = container_type(addedHistogram.size(), addedHistogram.getOffset());
      std::transform(addedHistogram.begin(), addedHistogram.end(), this->mContainer.begin(), [this, frequencyCountingDecorator](count_t frequency) {
        return frequencyCountingDecorator(frequency);
      });
    } else {
      const symbol_t newSize = newMax - newMin + 1;
      typename container_type::container_type newFreequencyTable(newSize, 0);
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

      this->mContainer = container_type{std::move(newFreequencyTable), static_cast<source_type>(newHistogram.getOffset())};
    }
  }

  return *this;
}

template <typename source_T>
auto FrequencyTable<source_T, std::enable_if_t<sizeof(source_T) == 4>>::resize(source_type min, source_type max) -> FrequencyTable&
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
  const source_type oldOffset = this->getOffset();
  this->mNSamples = 0;

  if (this->mContainer.empty()) {
    this->mContainer = container_type{newSize, min};
    return *this;
  } else {
    container_type oldFrequencyTable = std::move(this->mContainer);
    auto oldHistogram = utils::HistogramView{oldFrequencyTable.begin(), oldFrequencyTable.end(), oldOffset};
    this->mContainer = container_type{newSize, min};
    return this->addFrequencies(oldHistogram.begin(), oldHistogram.end(), oldHistogram.getMin());
  }
}

template <typename source_T>
class FrequencyTable<source_T, std::enable_if_t<sizeof(source_T) <= 2>> : public internal::FrequencyContainer<source_T>,
                                                                          public FrequencyTableBase<source_T,
                                                                                                    typename internal::FrequencyContainer<source_T>::value_type,
                                                                                                    FrequencyTable<source_T>>
{
  using containerBase_type = internal::FrequencyContainer<source_T>;
  using frequencyTableBase_type = FrequencyTableBase<source_T, typename internal::FrequencyContainer<source_T>::value_type, FrequencyTable<source_T>>;

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

  FrequencyTable() = default;

  template <typename freq_IT>
  FrequencyTable(freq_IT begin, freq_IT end, source_type offset) : containerBase_type(), frequencyTableBase_type{begin, end, offset} {};

  // operations
  template <typename source_IT>
  FrequencyTable& addSamples(source_IT begin, source_IT end);

  FrequencyTable& addSamples(gsl::span<const source_type> samples);

  template <typename freq_IT>
  FrequencyTable& addFrequencies(freq_IT begin, freq_IT end, source_type offset);

  using frequencyTableBase_type::addFrequencies;

  friend void swap(FrequencyTable& a, FrequencyTable& b) noexcept
  {
    using std::swap;
    swap(static_cast<typename FrequencyTable::containerBase_type&>(a),
         static_cast<typename FrequencyTable::containerBase_type&>(b));
  };
};

template <typename source_T>
template <typename source_IT>
auto FrequencyTable<source_T, std::enable_if_t<sizeof(source_T) <= 2>>::addSamples(source_IT begin, source_IT end) -> FrequencyTable&
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
auto FrequencyTable<source_T, std::enable_if_t<sizeof(source_T) <= 2>>::addSamples(gsl::span<const source_type> samples) -> FrequencyTable&
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

    std::array<internal::ShiftedVector<source_type, value_type>, 3> histograms{
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
auto FrequencyTable<source_T, std::enable_if_t<sizeof(source_T) <= 2>>::addFrequencies(freq_IT begin, freq_IT end, source_type offset) -> FrequencyTable&
{
  // bounds check
  utils::HistogramView addedHistogram{begin, end, offset};
  utils::HistogramView<typename container_type::iterator> thisHistogram{};
  thisHistogram = utils::HistogramView{this->mContainer.begin(), this->mContainer.end(), this->mContainer.getOffset()};
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
    this->mContainer[i] += *iter;
    ++iter;
  };
  return *this;
}

} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_FREQUENCYTABLE_H_ */
