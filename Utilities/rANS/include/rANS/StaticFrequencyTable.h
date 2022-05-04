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

/// @file   StaticFrequencyTable.h
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

#include <fairlogger/Logger.h>

#include "rANS/definitions.h"
#include "rANS/internal/helper.h"
#include "rANS/utils/HistogramView.h"

namespace o2
{
namespace rans
{

namespace impl
{

inline uint64_t load64(const void* __restrict src)
{
  uint64_t ret;
  std::memcpy(&ret, src, 8);
  return ret;
};

} // namespace impl

template <typename source_T>
class StaticFrequencyTable
{
 public:
  using symbol_t = source_T;
  using iterator_t = count_t*;
  using constIterator_t = const count_t*;

  // Constructors

  // TODO(milettri): fix once ROOT cling respects the standard http://wg21.link/p1286r2
  StaticFrequencyTable() noexcept : mFrequencyTable(1ull << (sizeof(source_T) * 8), 0){}; // NOLINT

  template <typename Freq_IT, std::enable_if_t<internal::isIntegralIter_v<Freq_IT>, bool> = true>
  StaticFrequencyTable(Freq_IT begin, Freq_IT end, symbol_t min);

  // accessors

  inline count_t operator[](symbol_t symbol) const
  {
    return getSymbol(symbol);
  };

  count_t at(size_t index) const;

  inline const count_t* data() const noexcept { return mFrequencyTable.data(); };

  inline constIterator_t cbegin() const noexcept { return data(); };

  inline constIterator_t cend() const noexcept { return data() + size(); };

  inline constIterator_t begin() const noexcept { return cbegin(); };

  inline constIterator_t end() const noexcept { return cend(); };

  inline iterator_t begin() noexcept { return const_cast<iterator_t>(cbegin()); };

  inline iterator_t end() noexcept { return const_cast<iterator_t>(cend()); };

  inline size_t size() const noexcept { return mFrequencyTable.size(); };

  size_t getNUsedAlphabetSymbols() const;

  inline size_t getAlphabetRangeBits() const noexcept { return internal::numBitsForNSymbols(this->size()); };

  inline symbol_t getMinSymbol() const noexcept { return std::numeric_limits<source_T>::min(); };
  inline symbol_t getMaxSymbol() const noexcept { return std::numeric_limits<source_T>::max(); };

  inline size_t getNumSamples() const noexcept { return mNumSamples; };

  // operations
  template <typename Source_IT, std::enable_if_t<internal::isCompatibleIter_v<source_T, Source_IT>, bool> = true>
  StaticFrequencyTable& addSamples(Source_IT begin, Source_IT end);

  template <typename Freq_IT, std::enable_if_t<internal::isIntegralIter_v<Freq_IT>, bool> = true>
  StaticFrequencyTable& addFrequencies(Freq_IT begin, Freq_IT end, symbol_t min);

  StaticFrequencyTable& operator+(StaticFrequencyTable& other);

  histogram_t release() && noexcept;

 private:
  const count_t& getSymbol(symbol_t symbol) const;
  count_t& getSymbol(symbol_t symbol);

  count_t frequencyCountingDecorator(count_t frequency) noexcept;

  histogram_t mFrequencyTable{};
  size_t mNumSamples{};
}; // namespace rans

} // namespace rans
} // namespace o2

// IMPL
///////////////////////////////////////////////////////////////////////////////////////////

namespace o2
{
namespace rans
{

template <typename source_T>
template <typename Freq_IT, std::enable_if_t<internal::isIntegralIter_v<Freq_IT>, bool>>
StaticFrequencyTable<source_T>::StaticFrequencyTable(Freq_IT begin, Freq_IT end, symbol_t min) : StaticFrequencyTable()
{

  auto histogram = utils::HistogramView{begin, end, min};
  this->addFrequencies(histogram.begin(), histogram.end(), histogram.getMin(), true);
}

template <typename source_T>
template <typename Source_IT, std::enable_if_t<internal::isCompatibleIter_v<source_T, Source_IT>, bool>>
inline StaticFrequencyTable<source_T>& StaticFrequencyTable<source_T>::addSamples(Source_IT begin, Source_IT end)
{
  if (begin == end) {
    LOG(warning) << "Passed empty message to " << __func__;
  } else {
    const auto inputSize = std::distance(begin, end);
    auto iter = begin;

    constexpr size_t iterSize = sizeof(uint64_t) / sizeof(source_T);

    count_t hists[3][256]{0};

    if constexpr (iterSize == 8) {
      while (iter != begin + (inputSize & ~(2 * iterSize - 1))) {
        mNumSamples += 2 * iterSize;
        uint64_t in0 = impl::load64(&(*iter));
        uint64_t in1 = impl::load64(&(*iter) + iterSize);

        uint64_t i0 = in0;
        ++hists[0][(static_cast<uint8_t>(i0))];
        ++hists[1][(static_cast<uint16_t>(i0) >> 8)];
        i0 >>= 16;
        ++hists[2][(static_cast<uint8_t>(i0))];
        ++mFrequencyTable[(static_cast<uint16_t>(i0) >> 8)];
        i0 = in0 >>= 32;
        ++hists[0][(static_cast<uint8_t>(i0))];
        ++hists[1][(static_cast<uint16_t>(i0) >> 8)];
        i0 >>= 16;
        ++hists[2][(static_cast<uint8_t>(i0))];
        ++mFrequencyTable[(static_cast<uint16_t>(i0) >> 8)];

        uint64_t i1 = in1;
        ++hists[0][(static_cast<uint8_t>(i1))];
        ++hists[1][(static_cast<uint16_t>(i1) >> 8)];
        i1 >>= 16;
        ++hists[2][(static_cast<uint8_t>(i1))];
        ++mFrequencyTable[(static_cast<uint16_t>(i1) >> 8)];
        i1 = in1 >>= 32;
        ++hists[0][(static_cast<uint8_t>(i1))];
        ++hists[1][(static_cast<uint16_t>(i1) >> 8)];
        i1 >>= 16;
        ++hists[2][(static_cast<uint8_t>(i1))];
        ++mFrequencyTable[(static_cast<uint16_t>(i1) >> 8)];

        iter += 2 * iterSize;

        __builtin_prefetch(&(*iter) + 512, 0);
      }
      while (iter != end) {
        ++mNumSamples;
        ++this->getSymbol(*iter++);
      }

      for (size_t i = 0; i < 256; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          mFrequencyTable[i] += hists[j][i];
        }
      }
    } else {
      auto secondFreq = mFrequencyTable;

      while (iter != begin + (inputSize & ~(2 * iterSize - 1))) {
        mNumSamples += 2 * iterSize;
        uint64_t in0 = impl::load64(&(*iter));
        uint64_t in1 = impl::load64(&(*iter) + iterSize);

        // 32 bit
        if constexpr (iterSize == 2) {
          ++this->getSymbol(static_cast<source_T>(in0));
          ++this->getSymbol(static_cast<source_T>(in0) >> 32);
          ++this->getSymbol(static_cast<source_T>(in1));
          ++this->getSymbol(static_cast<source_T>(in1) >> 32);

        } else if constexpr (iterSize == 4) {
          uint64_t i0 = in0;
          ++mFrequencyTable[static_cast<uint16_t>(i0)];
          ++secondFreq[static_cast<uint32_t>(i0) >> 16];
          i0 = in0 >> 32;
          ++mFrequencyTable[static_cast<uint16_t>(i0)];
          ++secondFreq[static_cast<uint32_t>(i0) >> 16];

          uint64_t i1 = in1;
          ++mFrequencyTable[static_cast<uint16_t>(i1)];
          ++secondFreq[static_cast<uint32_t>(i1) >> 16];
          i1 = in1 >> 32;
          ++mFrequencyTable[static_cast<uint16_t>(i1)];
          ++secondFreq[static_cast<uint32_t>(i1) >> 16];
        }
        if constexpr (iterSize == 8) {
          uint64_t i0 = in0;
          ++mFrequencyTable[(static_cast<uint8_t>(i0))];
          ++secondFreq[(static_cast<uint16_t>(i0) >> 8)];
          i0 >>= 16;
          ++mFrequencyTable[(static_cast<uint8_t>(i0))];
          ++secondFreq[(static_cast<uint16_t>(i0) >> 8)];
          i0 = in0 >>= 32;
          ++mFrequencyTable[(static_cast<uint8_t>(i0))];
          ++secondFreq[(static_cast<uint16_t>(i0) >> 8)];
          i0 >>= 16;
          ++mFrequencyTable[(static_cast<uint8_t>(i0))];
          ++secondFreq[(static_cast<uint16_t>(i0) >> 8)];

          uint64_t i1 = in1;
          ++mFrequencyTable[(static_cast<uint8_t>(i1))];
          ++secondFreq[(static_cast<uint16_t>(i1) >> 8)];
          i1 >>= 16;
          ++mFrequencyTable[(static_cast<uint8_t>(i1))];
          ++secondFreq[(static_cast<uint16_t>(i1) >> 8)];
          i1 = in1 >>= 32;
          ++mFrequencyTable[(static_cast<uint8_t>(i1))];
          ++secondFreq[(static_cast<uint16_t>(i1) >> 8)];
          i1 >>= 16;
          ++mFrequencyTable[(static_cast<uint8_t>(i1))];
          ++secondFreq[(static_cast<uint16_t>(i1) >> 8)];
        }

        iter += 2 * iterSize;

        // ++this->getSymbol(*iter++);
        // ++this->getSymbol(*iter++);
        // ++this->getSymbol(*iter++);
        // ++this->getSymbol(*iter++);
        // ++this->getSymbol(*iter++);
        // ++this->getSymbol(*iter++);
        // ++this->getSymbol(*iter++);
        // ++this->getSymbol(*iter++);

        __builtin_prefetch(&(*iter) + 512, 0);
      }
      while (iter != end) {
        ++mNumSamples;
        ++this->getSymbol(*iter++);
      }

      for (size_t i = 0; i < mFrequencyTable.size(); ++i) {
        mFrequencyTable[i] += secondFreq[i];
      }
    }
  }
  return *this;
}

template <typename source_T>
template <typename Freq_IT, std::enable_if_t<internal::isIntegralIter_v<Freq_IT>, bool>>
StaticFrequencyTable<source_T>& StaticFrequencyTable<source_T>::addFrequencies(Freq_IT begin, Freq_IT end, symbol_t min)
{
  // bounds check
  utils::HistogramView addedHistogram{begin, end, min};
  utils::HistogramView thisHistogram{mFrequencyTable.begin(), mFrequencyTable.end(), this->getMinSymbol()};
  const bool invalidBounds = (utils::leftOffset(thisHistogram, addedHistogram) < 0) || (utils::rightOffset(thisHistogram, addedHistogram) > 0);

  if (invalidBounds) {
    throw std::runtime_error(fmt::format("Incompatible Frequency table dimensions: Cannot add [{},{}] to [{}, {}] ",
                                         addedHistogram.getMin(),
                                         addedHistogram.getMax(),
                                         thisHistogram.getMin(),
                                         thisHistogram.getMax()));
  }

  // intersection
  auto overlapAdded = utils::intersection(addedHistogram, thisHistogram);
  auto overlapThis = utils::intersection(thisHistogram, addedHistogram);
  if (!overlapAdded.empty()) {
    assert(overlapAdded.getMin() == overlapThis.getMin());
    assert(overlapAdded.size() == overlapThis.size());
    std::transform(overlapAdded.begin(), overlapAdded.end(), overlapThis.begin(), overlapThis.begin(), [this](const count_t& a, const count_t& b) { return this->frequencyCountingDecorator(a) + b; });
  }
  return &this;
}

template <typename source_T>
inline auto StaticFrequencyTable<source_T>::at(size_t index) const -> count_t
{
  assert(index < size());
  return mFrequencyTable[index];
};

template <typename source_T>
inline StaticFrequencyTable<source_T>& StaticFrequencyTable<source_T>::operator+(StaticFrequencyTable<source_T>& other)
{
  addFrequencies(other.cbegin(), other.cend(), other.getMinSymbol());
  return *this;
}

template <typename source_T>
inline auto StaticFrequencyTable<source_T>::release() && noexcept -> histogram_t
{
  auto frequencies = std::move(mFrequencyTable);

  return frequencies;
};

template <typename source_T>
inline size_t StaticFrequencyTable<source_T>::getNUsedAlphabetSymbols() const
{
  return std::count_if(mFrequencyTable.begin(), mFrequencyTable.end(), [](count_t count) { return count > 0; }) + static_cast<count_t>(this->hasIncompressibleSymbols());
};

template <typename source_T>
inline auto StaticFrequencyTable<source_T>::getSymbol(symbol_t symbol) const -> const count_t&
{
  if constexpr (std::is_signed_v<source_T>) {
    symbol -= getMinSymbol();
  };
  return mFrequencyTable[symbol];
}
template <typename source_T>
inline auto StaticFrequencyTable<source_T>::getSymbol(symbol_t symbol) -> count_t&
{
  return const_cast<count_t&>(static_cast<const StaticFrequencyTable&>(*this).getSymbol(symbol));
}

template <typename source_T>
inline count_t StaticFrequencyTable<source_T>::frequencyCountingDecorator(count_t frequency) noexcept
{
  mNumSamples += frequency;
  return frequency;
};

} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_STATICFREQUENCYTABLE_H_ */
