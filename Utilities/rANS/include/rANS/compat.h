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

/// @file   compat.h
/// @author Michael Lettrich
/// @brief  functionality to maintain compatibility with previous version of this library

#ifndef RANS_COMPAT_H_
#define RANS_COMPAT_H_

#include <numeric>

#include "rANS/histogram.h"

namespace o2::rans::compat
{

namespace compatImpl
{
inline constexpr uint32_t MinRenormThreshold = 10;
inline constexpr uint32_t MaxRenormThreshold = 20;

size_t computeRenormingPrecision(size_t nUsedAlphabetSymbols)
{
  const uint32_t minBits = o2::rans::internal::log2UInt(nUsedAlphabetSymbols);
  const uint32_t estimate = minBits * 3u / 2u;
  const uint32_t maxThreshold = std::max(minBits, MaxRenormThreshold);
  const uint32_t minThreshold = std::max(estimate, MinRenormThreshold);

  return std::min(minThreshold, maxThreshold);
};
} // namespace compatImpl

template <typename source_T>
RenormedHistogram<source_T> renorm(Histogram<source_T> histogram, size_t newPrecision = 0)
{
  using namespace o2::rans::internal;

  constexpr size_t IncompressibleSymbolFrequency = 1;

  if (histogram.empty()) {
    LOG(warning) << "rescaling empty histogram";
  }

  size_t nUsedAlphabetSymbols = histogram.countNUsedAlphabetSymbols();

  if (newPrecision == 0) {
    newPrecision = compatImpl::computeRenormingPrecision(nUsedAlphabetSymbols);
  }

  std::vector<uint64_t> cumulativeFrequencies(histogram.size() + 2);
  const size_t alphabetSize = histogram.size() + 1;
  cumulativeFrequencies[0] = 0;
  std::inclusive_scan(histogram.begin(), histogram.end(), ++cumulativeFrequencies.begin(), std::plus<>(), 0ull);
  cumulativeFrequencies.back() = histogram.getNumSamples() + IncompressibleSymbolFrequency;

  auto getFrequency = [&cumulativeFrequencies](count_t i) {
    assert(cumulativeFrequencies[i + 1] >= cumulativeFrequencies[i]);
    return cumulativeFrequencies[i + 1] - cumulativeFrequencies[i];
  };

  // we will memorize only those entries which can be used
  const auto sortIdx = [&]() {
    std::vector<size_t> indices;
    indices.reserve(nUsedAlphabetSymbols + 1);

    for (size_t i = 0; i < alphabetSize; ++i) {
      if (getFrequency(i) != 0) {
        indices.push_back(i);
      }
    }

    std::sort(indices.begin(), indices.end(), [&](count_t i, count_t j) { return getFrequency(i) < getFrequency(j); });
    return indices;
  }();

  // resample distribution based on cumulative frequencies
  const count_t newCumulatedFrequency = pow2(newPrecision);
  size_t nSamples = histogram.getNumSamples() + 1;
  assert(newCumulatedFrequency >= nUsedAlphabetSymbols);
  size_t needsShift = 0;
  for (size_t i = 0; i < sortIdx.size(); i++) {
    if (static_cast<count_t>(getFrequency(sortIdx[i])) * (newCumulatedFrequency - needsShift) / nSamples >= 1) {
      break;
    }
    needsShift++;
  }

  size_t shift = 0;
  auto beforeUpdate = cumulativeFrequencies[0];
  for (size_t i = 0; i < alphabetSize; i++) {
    auto& nextCumulative = cumulativeFrequencies[i + 1];
    uint64_t oldFrequeny = nextCumulative - beforeUpdate;
    if (oldFrequeny && oldFrequeny * (newCumulatedFrequency - needsShift) / nSamples < 1) {
      shift++;
    }
    beforeUpdate = cumulativeFrequencies[i + 1];
    nextCumulative = (static_cast<uint64_t>(newCumulatedFrequency - needsShift) * nextCumulative) / nSamples + shift;
  }
  assert(shift == needsShift);

  // verify
#if !defined(NDEBUG)
  assert(cumulativeFrequencies.front() == 0);
  assert(cumulativeFrequencies.back() == newCumulatedFrequency);
  size_t i = 0;
  for (auto frequency : histogram) {
    if (frequency == 0) {
      assert(cumulativeFrequencies[i + 1] == cumulativeFrequencies[i]);
    } else {
      assert(cumulativeFrequencies[i + 1] > cumulativeFrequencies[i]);
    }
    ++i;
  }
#endif

  typename RenormedHistogram<source_T>::container_type rescaledFrequencies(histogram.size(), histogram.getOffset());

  // calculate updated frequencies
  for (size_t i = 0; i < histogram.size(); i++) {
    rescaledFrequencies[i] = getFrequency(i);
  }
  const typename RenormedHistogram<source_T>::value_type incompressibleSymbolFrequency = getFrequency(histogram.size());

  return RenormedHistogram<source_T>{std::move(rescaledFrequencies), newPrecision, incompressibleSymbolFrequency};
};
} // namespace o2::rans::compat

#endif /* RANS_COMPAT_H_ */