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

#include "rANS/internal/FrequencyTableBase.h"
#include "rANS/internal/DynamicFrequencyContainer.h"
#include "rANS/utils/HistogramView.h"

namespace o2
{
namespace rans
{

template <typename source_T>
class DynamicFrequencyTable : public DynamicFrequencyContainer<source_T>,
                              public FrequencyTableBase<source_T,
                                                        typename DynamicFrequencyContainer<source_T>::value_type,
                                                        DynamicFrequencyTable<source_T>>
{
  using containerBase_type = DynamicFrequencyContainer<source_T>;
  using frequencyTableBase_type = FrequencyTableBase<source_T, typename DynamicFrequencyContainer<source_T>::value_type, DynamicFrequencyTable<source_T>>;

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

  using frequencyTableBase_type::addSamples;

  template <typename source_IT>
  DynamicFrequencyTable& addSamples(source_IT begin, source_IT end);

  template <typename source_IT>
  DynamicFrequencyTable& addSamples(source_IT begin, source_IT end, source_type min, source_type max);

  template <typename source_IT>
  DynamicFrequencyTable& addSamples(gsl::span<const source_type> span, source_type min, source_type max);

  // operations

  using frequencyTableBase_type::addFrequencies;

  template <typename freq_IT>
  DynamicFrequencyTable& addFrequencies(freq_IT begin, freq_IT end, source_type offset);

  DynamicFrequencyTable& resize(source_type min, source_type max);

  inline DynamicFrequencyTable& resize(size_type newSize) { return resize(this->getMinSymbol(), this->getMinSymbol() + newSize); };
};

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
template <typename source_IT>
inline auto DynamicFrequencyTable<source_T>::addSamples(gsl::span<const source_type> span, source_type min, source_type max) -> DynamicFrequencyTable&
{
  return addSamples(span.begin(), span.end(), min, max);
}

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
