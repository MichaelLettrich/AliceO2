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

/// @file   renorm.h
/// @author Michael Lettrich
/// @brief  Renorm histogram to sum of frequencies = 2^P for use in fast rans coding. Includes estimation of P.

#ifndef RANS_INTERNAL_TRANSFORM_RENORM_H_
#define RANS_INTERNAL_TRANSFORM_RENORM_H_

#include <fairlogger/Logger.h>

#include "rANS/internal/containers/RenormedHistogram.h"
#include "rANS/internal/containers/Histogram.h"
#include "rANS/internal/containers/HistogramView.h"
#include "rANS/internal/common/utils.h"

namespace o2::rans
{

namespace internal
{
inline constexpr size_t MinRenormThreshold = 10;
inline constexpr size_t MaxRenormThreshold = 20;
} // namespace internal

template <typename source_T>
struct DatasetMetrics {
  source_T min{};
  source_T max{};
  uint32_t alphabetRangeBits{};
  uint32_t nUsedAlphabetSymbols{};
  float_t entropy{};
  std::array<uint32_t, 32> symbolLengthDistribution{{}};
  std::array<float_t, 32> weightedSymbolLengthDistribution{{}};
};

template <typename source_T, template <class> class Histogram_T>
DatasetMetrics<source_T> computeDatasetMetrics(const Histogram_T<source_T>& frequencyTable)
{
  using namespace internal;
  using source_type = source_T;
  DatasetMetrics<source_type> metrics{};

  const auto trimmedFrequencyView = trim(HistogramView{frequencyTable.begin(), frequencyTable.end()});
  metrics.min = trimmedFrequencyView.getMin();
  metrics.max = trimmedFrequencyView.getMax();
  assert(metrics.max >= metrics.min);
  metrics.alphabetRangeBits = log2UInt(static_cast<uint32_t>(metrics.max - metrics.min));

  for (auto iter = trimmedFrequencyView.begin(); iter != trimmedFrequencyView.end(); ++iter) {
    auto frequency = *iter;
    if (frequency) {
      ++metrics.nUsedAlphabetSymbols;

      const float_t probability = static_cast<float_t>(frequency) / static_cast<float_t>(frequencyTable.getNumSamples());
      const float_t length = std::log2(probability);

      metrics.entropy -= probability * length;
      ++metrics.symbolLengthDistribution[static_cast<uint32_t>(-length)];
      metrics.weightedSymbolLengthDistribution[static_cast<uint32_t>(-length)] += probability;
    }
  }
  return metrics;
}

template <typename source_T>
double_t computeExpectedCodewordLength(const Histogram<source_T>& histogram, const RenormedHistogram<source_T>& rescaledHistogram)
{
  using namespace internal;
  using value_type = typename Histogram<source_T>::value_type;

  HistogramView histogramView{histogram.begin(), histogram.end(), histogram.getOffset()};
  HistogramView renormedView{rescaledHistogram.begin(), rescaledHistogram.end(), rescaledHistogram.getOffset()};

  auto getRescaledFrequency = [&renormedView](source_T sourceSymbol) -> value_type {
    if (sourceSymbol >= renormedView.getMin() && sourceSymbol <= renormedView.getMax()) {
      return renormedView[sourceSymbol];
    } else {
      return static_cast<value_type>(0);
    }
  };

  double_t expectedCodewordLength = 0;
  value_type trueIncompressibleFrequency = 0;

  assert(histogram.countNUsedAlphabetSymbols() >= rescaledHistogram.countNUsedAlphabetSymbols());

  // all "normal symbols"
  for (value_type sourceSymbol = histogramView.getMin(); sourceSymbol <= histogramView.getMax(); ++sourceSymbol) {

    const value_type frequency = histogramView[sourceSymbol];
    if (frequency) {
      const value_type rescaledFrequency = getRescaledFrequency(sourceSymbol);

      const double_t trueProbability = static_cast<double_t>(frequency) / histogram.getNumSamples();

      if (rescaledFrequency) {
        const double_t rescaledProbability = static_cast<double_t>(rescaledFrequency) / rescaledHistogram.getNumSamples();
        expectedCodewordLength -= trueProbability * std::log2(rescaledProbability);
      } else {
        trueIncompressibleFrequency += frequency;
      }
    }
  }
  // incompressibleSymbol:
  const double_t trueProbability = static_cast<double_t>(trueIncompressibleFrequency) / histogram.getNumSamples();
  const double_t rescaledProbability = static_cast<double_t>(rescaledHistogram.getIncompressibleSymbolFrequency()) / rescaledHistogram.getNumSamples();

  expectedCodewordLength -= trueProbability * std::log2(rescaledProbability);
  expectedCodewordLength += trueProbability * std::log2(toBits(sizeof(source_T)));

  return expectedCodewordLength;
};

template <typename source_T>
inline constexpr bool preferPacking(const DatasetMetrics<source_T>& metrics, float_t threshold = 0.1) noexcept
{
  if (metrics.entropy == 0) {
    return true;
  } else {
    const float_t ratio = static_cast<float_t>(metrics.alphabetRangeBits) / static_cast<float_t>(metrics.entropy);
    return ratio < (1.0f + threshold);
  }
};

template <typename source_T>
inline constexpr size_t computeRenormingPrecision(const DatasetMetrics<source_T>& metrics, float_t cutoffPrecision = 0.999) noexcept
{
  size_t computedRenormingBits = [&]() -> size_t {
    float_t cumulatedPrecision = 0;
    size_t renormingBits = 0;
    for (; renormingBits < metrics.weightedSymbolLengthDistribution.size() && cumulatedPrecision < cutoffPrecision; ++renormingBits) {
      cumulatedPrecision += metrics.weightedSymbolLengthDistribution[renormingBits];
    }

    if (cumulatedPrecision == 0) {
      return 0;
    } else {
      return renormingBits;
    }
  }();

  // select renorming from [MinThreshold, MaxThreshold]
  return std::max(internal::MinRenormThreshold, std::min(internal::MaxRenormThreshold, computedRenormingBits + 1));
};

template <typename source_T>
size_t estimateEncodeBufferBytes(const DatasetMetrics<source_T>& metrics, const Histogram<source_T>& frequencyTable, float_t safetyMargin = 1.2f)
{
  float_t estimatedSize = metrics.entropy * frequencyTable.getNumSamples(); //min Size in Bits, based on entropy
  estimatedSize *= safetyMargin;                                            // add safety margin;
  return internal::toBytes(std::ceil(estimatedSize));
};

template <typename source_T>
size_t estimateDictBytes(const DatasetMetrics<source_T>& metrics, float_t safetyMargin = 1.2f)
{

  auto computeDeltaCodedSize = [](uint64_t x) -> uint64_t {
    return internal::log2UInt(x) + 2ul * internal::log2UInt(internal::log2UInt(x) + 1ul) + 1ul;
  };

  constexpr size_t maxBitsFrequency = computeDeltaCodedSize(internal::MaxRenormThreshold);
  constexpr size_t maxBitsDelta = computeDeltaCodedSize(internal::toBits(sizeof(source_T)));
  constexpr size_t generalOverhead = 0;

  float_t estimatedSize = (maxBitsFrequency + maxBitsDelta) * metrics.nUsedAlphabetSymbols;
  estimatedSize += generalOverhead;
  estimatedSize *= safetyMargin;
  return internal::toBytes(std::ceil(estimatedSize));
};

template <typename source_T>
inline constexpr size_t computeRenormingPrecision(size_t nUsedAlphabetSymbols)
{
  if (nUsedAlphabetSymbols == 0) {
    return 0;
  }

  if constexpr (sizeof(source_T) == 1) {
    return 14;
  } else {
    const size_t minBits = internal::log2UInt(nUsedAlphabetSymbols);
    const size_t estimate = minBits * 3u / 2u;
    const size_t maxThreshold = std::max(minBits, internal::MaxRenormThreshold);
    const size_t minThreshold = std::max(estimate, internal::MinRenormThreshold);

    return std::min(minThreshold, maxThreshold);
  };
}

namespace renormImpl
{

template <typename source_T>
inline size_t getNUsedAlphabetSymbols(const Histogram<source_T>& f)
{
  if constexpr (sizeof(source_T) <= 2) {
    const size_t nUsedAlphabetSymbols = f.empty() ? 0 : f.size();
    return nUsedAlphabetSymbols;
  } else {
    return f.countNUsedAlphabetSymbols();
  }
}

} // namespace renormImpl

template <typename source_T>
RenormedHistogram<source_T> renormCutoffIncompressible(Histogram<source_T> histogram, uint8_t newPrecision = 0, uint8_t lowProbabilityCutoffBits = 3)
{
  if (histogram.empty()) {
    LOG(warning) << "rescaling Frequency Table for empty message";
  }

  using source_type = source_T;
  using count_type = typename Histogram<source_T>::value_type;
  using difference_type = typename Histogram<source_T>::difference_type;
  using container_type = typename Histogram<source_T>::container_type;
  using iterator_type = typename container_type::iterator;

  const source_type offset = histogram.getOffset();
  const double_t nSamples = histogram.getNumSamples();
  const count_type nUsedAlphabetSymbols = renormImpl::getNUsedAlphabetSymbols(histogram);

  if (newPrecision == 0) {
    newPrecision = computeRenormingPrecision<source_type>(nUsedAlphabetSymbols);
  }

  const count_type nSamplesRescaled = 1 << newPrecision;
  const double_t probabilityCutOffThreshold = 1 / static_cast<double_t>(1ul << (newPrecision + lowProbabilityCutoffBits));

  // scaling
  double_t incompressibleSymbolProbability = 0;
  count_type nSamplesRescaledUncorrected = 0;
  std::vector<iterator_type> correctableIndices;
  correctableIndices.reserve(nUsedAlphabetSymbols);

  auto scaleFrequency = [nSamplesRescaled](double_t symbolProbability) -> double_t { return symbolProbability * nSamplesRescaled; };
  auto roundDownFrequency = [](double_t x) -> count_type { return static_cast<count_type>(x); };
  auto roundUpFrequency = [roundDownFrequency](double_t x) -> count_type { return roundDownFrequency(x) + 1; };
  auto roundFrequency = [roundDownFrequency, roundUpFrequency](double_t rescaledFrequency) -> count_type {
    if (rescaledFrequency * rescaledFrequency <= (roundDownFrequency(rescaledFrequency) * roundUpFrequency(rescaledFrequency))) {
      return roundDownFrequency(rescaledFrequency);
    } else {
      return roundUpFrequency(rescaledFrequency);
    }
  };

  container_type rescaledHistogram = std::move(histogram).release();

  for (auto frequencyIter = rescaledHistogram.begin(); frequencyIter != rescaledHistogram.end(); ++frequencyIter) {
    const count_type frequency = *frequencyIter;
    // LOGP(info, "freq: {}", frequency);
    if (frequency > 0) {
      const double_t symbolProbability = static_cast<double_t>(frequency) / nSamples;
      if (symbolProbability < probabilityCutOffThreshold) {
        incompressibleSymbolProbability += symbolProbability;
        *frequencyIter = 0;
      } else {
        const double_t scaledFrequencyD = scaleFrequency(symbolProbability);
        count_type rescaledFrequency = roundFrequency(scaledFrequencyD);
        assert(rescaledFrequency > 0);
        *frequencyIter = rescaledFrequency;
        nSamplesRescaledUncorrected += rescaledFrequency;
        if (rescaledFrequency > 1) {
          correctableIndices.push_back(frequencyIter);
        }
      }
    }
  }

  // treat incompressible symbol
  const count_type incompressibleSymbolFrequency = std::max(static_cast<count_type>(1), static_cast<count_type>(incompressibleSymbolProbability * nSamplesRescaled));
  nSamplesRescaledUncorrected += incompressibleSymbolFrequency;

  // correction
  std::stable_sort(correctableIndices.begin(), correctableIndices.end(), [&rescaledHistogram](const iterator_type& a, const iterator_type& b) { return *a < *b; });

  difference_type nCorrections = static_cast<difference_type>(nSamplesRescaled) - static_cast<difference_type>(nSamplesRescaledUncorrected);
  const double_t rescalingFactor = static_cast<double_t>(nSamplesRescaled) / static_cast<double_t>(nSamplesRescaledUncorrected);

  for (auto iter : correctableIndices) {
    if (std::abs(nCorrections) > 0) {
      const difference_type uncorrectedFrequency = *iter;
      difference_type correction = uncorrectedFrequency - roundFrequency(uncorrectedFrequency * rescalingFactor);

      if (nCorrections < 0) {
        // overshoot - correct downwards by subtracting correction in [1,|nCorrections|]
        correction = std::max(1l, std::min(correction, std::abs(nCorrections)));
      } else {
        // correct upwards by subtracting correction in [-1, -nCorrections]
        correction = std::min(-1l, std::max(correction, -nCorrections));
      }

      // the corrected frequency must be at least 1 though
      const count_type correctedFrequency = std::max(1l, uncorrectedFrequency - correction);
      nCorrections += uncorrectedFrequency - correctedFrequency;
      *iter = correctedFrequency;
    } else {
      break;
    }
  }

  if (std::abs(nCorrections) > 0) {
    throw std::runtime_error(fmt::format("rANS rescaling incomplete: {} corrections Remaining", nCorrections));
  }

  return RenormedHistogram<source_type>(std::move(rescaledHistogram), newPrecision, incompressibleSymbolFrequency);
};

} // namespace o2::rans

#endif /* RANS_INTERNAL_TRANSFORM_RENORM_H_ */
