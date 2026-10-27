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

/// @file   utils.h
/// @author Michael Lettrich
/// @brief  utility functions related to calculating different dataset metrics

#ifndef RANS_INTERNAL_METRICS_UTILS_H_
#define RANS_INTERNAL_METRICS_UTILS_H_

#include <cstddef>
#include <cmath>
#include <cstdint>
#include <numeric>

#include "rANS/internal/common/utils.h"
#include "rANS/internal/containers/HistogramView.h"

namespace o2::rans
{

template <class container_T, typename renormed_container_T>
double_t computeExpectedCodewordLength(const container_T& histogram, const renormed_container_T& rescaledHistogram)
{
  assert(histogram.getNumSamples() > 0);
  assert(rescaledHistogram.getNumSamples() > 0);

  using namespace internal;
  using namespace utils;

  using source_type = typename container_T::source_type;
  using value_type = typename container_T::value_type;

  const auto histogramView = makeHistogramView(histogram);
  const auto renormedView = makeHistogramView(rescaledHistogram);

  auto getRescaledFrequency = [&renormedView](source_type sourceSymbol) -> value_type {
    if (sourceSymbol >= renormedView.getMin() && sourceSymbol <= renormedView.getMax()) {
      return renormedView[sourceSymbol];
    } else {
      return static_cast<value_type>(0);
    }
  };

  double_t expectedCodewordLength{0.0};
  value_type trueIncompressibleFrequency{0};

  assert(countNUsedAlphabetSymbols(histogram) >= countNUsedAlphabetSymbols(rescaledHistogram));

  const double_t reciprocalNumSamples = 1.0 / histogram.getNumSamples();
  const double_t reciprocalNumSamplesRescaled = 1.0 / rescaledHistogram.getNumSamples();
  source_type incompressibleMin = std::numeric_limits<source_type>::max();
  source_type incompressibleMax = std::numeric_limits<source_type>::min();

  // all "normal symbols"
  for (int64_t sourceSymbol = histogramView.getMin(); sourceSymbol < histogramView.getMax() + 1; ++sourceSymbol) {

    const value_type frequency = histogramView[sourceSymbol];
    if (frequency) {
      const value_type rescaledFrequency = getRescaledFrequency(sourceSymbol);
      const double_t trueProbability = static_cast<double_t>(frequency) * reciprocalNumSamples;

      if (rescaledFrequency) {
        const double_t rescaledProbability = static_cast<double_t>(rescaledFrequency) * reciprocalNumSamplesRescaled;
        expectedCodewordLength -= trueProbability * fastlog2(rescaledProbability);
      } else {
        trueIncompressibleFrequency += frequency;
        incompressibleMin = std::min<source_type>(incompressibleMin, sourceSymbol);
        incompressibleMax = std::max<source_type>(incompressibleMax, sourceSymbol);
      }
    }
  }

  // incompressibleSymbol:
  const double_t trueIncompressibleProbability = static_cast<double_t>(trueIncompressibleFrequency) * reciprocalNumSamples;
  if (trueIncompressibleProbability) {
    const double_t rescaledProbability = static_cast<double_t>(rescaledHistogram.getIncompressibleSymbolFrequency()) * reciprocalNumSamplesRescaled;
    expectedCodewordLength -= trueIncompressibleProbability * fastlog2(rescaledProbability);
    expectedCodewordLength += trueIncompressibleProbability * fastlog2(static_cast<double_t>(incompressibleMax - incompressibleMin) + 1.0);
  }

  return expectedCodewordLength;
};

} // namespace o2::rans

#endif /* RANS_INTERNAL_METRICS_UTILS_H_ */