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

/// @file   RenormedFrequencyTable.h
/// @author Michael Lettrich
/// @since  2019-05-08
/// @brief Histogram to depict frequencies of source symbols for rANS compression.

#ifndef INCLUDE_RANS_RENORM_H_
#define INCLUDE_RANS_RENORM_H_

#include <fairlogger/Logger.h>

#include "rANS/definitions.h"

#include "rANS/RenormedFrequencies.h"
#include "rANS/StaticFrequencyTable.h"
#include "rANS/DynamicFrequencyTable.h"
#include "rANS/internal/helper.h"

namespace o2
{
namespace rans
{

template <typename source_T>
struct DatasetMetrics {
  source_T min{};
  source_T max{};
  uint32_t nUsedAlphabetSymbols{};
  float_t entropy{};
  std::array<uint32_t, 32> symbolLengthDistribution;
  std::array<float_t, 32> weightedSymbolLengthDistribution;
};

template <typename frequencyTable_T>
DatasetMetrics<typename frequencyTable_T::source_T> computeDatasetMetrics(const frequencyTable_T& container)
{
  using source_type = typename frequencyTable_T::source_T;
  DatasetMetrics<source_type> metrics{};

  for (auto iter = container.begin(); iter != container.end(); ++iter) {
    auto frequency = *iter;
    if (frequency) {
      ++metrics.nUsedAlphabetSymbols;

      const float_t probability = static_cast<float_t>(frequency) / static_cast<float_t>(container.getNumSamples());
      const float_t length = std::log2(probability);

      metrics.entropy -= probability * length;
      ++metrics.symbolLengthDistribution[static_cast<uint32_t>(length)];
      metrics.weightedSymbolLengthDistribution[static_cast<uint32_t>(length)] += static_cast<float_t>(frequency) * probability;
    }
  }

  // normalize weightedSymbolLengthDistribution
  for (float& elem : metrics.weightedSymbolLengthDistribution) {
    elem /= static_cast<float_t>(container.getNumSamples());
  }

  return metrics;
}

template <typename source_T>
inline constexpr bool shouldBePacked(const DatasetMetrics<source_T>& metrics, float_t threshold = 0.1) noexcept
{
  constexpr float_t alphabetRangeBits = internal::numBitsForNSymbols(metrics.max - metrics.min);

  return alphabetRangeBits / metrics.entropy > (1.0f + threshold);
};

template <typename source_T>
inline constexpr size_t computeRenormingPrecision(const DatasetMetrics<source_T>& metrics, float_t precision = 0.999) noexcept
{
  if constexpr (sizeof(source_T) == 1) {
    return MinRenormThreshold;
  } else {
    size_t computedRenormingBits = [&]() {
      size_t computedRenormingBits = 0;
      for (size_t i = 0; i < metrics.symbolLengthDistribution.size(); ++i) {
      }
    }();
  }
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
    const uint8_t minBits = internal::log2UInt(nUsedAlphabetSymbols);
    const uint8_t estimate = minBits * 3u / 2u;
    const uint8_t maxThreshold = std::max(minBits, MaxRenormThreshold);
    const uint8_t minThreshold = std::max(estimate, MinRenormThreshold);

    return std::min(minThreshold, maxThreshold);
  };
}

namespace renorm_impl
{

template <typename frequencyTable_T>
inline size_t getIncompressibleSymbolFrequency(const frequencyTable_T&)
{
  return 0;
}

template <typename source_T>
inline size_t getNUsedAlphabetSymbols(const DynamicFrequencyTable<source_T>& f)
{
  return f.computeNUsedAlphabetSymbols();
}

template <typename frequencyTable_T>
inline size_t getNUsedAlphabetSymbols(const frequencyTable_T& f)
{
  if (f.empty()) {
    return 0;
  } else {
    return f.size();
  }
}

template <typename frequencyTable_T>
inline uint32_t getOffset(const frequencyTable_T& f)
{
  return f.getOffset();
}

} // namespace renorm_impl

template <class frequencyTable_T>
decltype(auto) renormCutoffIncompressible(frequencyTable_T unrenomredTable, uint8_t newPrecision = 0, uint8_t lowProbabilityCutoffBits = 3)
{
  if (unrenomredTable.empty()) {
    LOG(warning) << "rescaling Frequency Table for empty message";
  }

  using source_type = typename frequencyTable_T::source_type;
  using index_type = typename frequencyTable_T::index_type;
  using count_type = typename frequencyTable_T::value_type;
  using difference_type = typename frequencyTable_T::difference_type;
  using container_type = typename frequencyTable_T::container_type;
  using iterator_type = typename container_type::iterator;

  const source_type offset = renorm_impl::getOffset(unrenomredTable);
  const double_t nSamples = unrenomredTable.getNumSamples();
  const double_t nIncompressibleSymbols = renorm_impl::getIncompressibleSymbolFrequency(unrenomredTable);
  const count_type nUsedAlphabetSymbols = renorm_impl::getNUsedAlphabetSymbols(unrenomredTable);

  if (newPrecision == 0) {
    newPrecision = computeRenormingPrecision<source_type>(nUsedAlphabetSymbols);
  }

  const count_type nSamplesRescaled = 1 << newPrecision;
  const double_t probabilityCutOffThreshold = 1 / static_cast<double_t>(1ul << (newPrecision + lowProbabilityCutoffBits));

  // scaling
  double_t incompressibleSymbolProbability = nIncompressibleSymbols / nSamplesRescaled;
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

  container_type rescaledFrequencies = std::move(unrenomredTable).release();

  for (auto frequencyIter = rescaledFrequencies.begin(); frequencyIter != rescaledFrequencies.end(); ++frequencyIter) {
    const count_type frequency = *frequencyIter;
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
  std::stable_sort(correctableIndices.begin(), correctableIndices.end(), [&rescaledFrequencies](const iterator_type& a, const iterator_type& b) { return *a < *b; });

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

  if constexpr (std::is_same_v<frequencyTable_T, StaticFrequencyTable<source_type>>) {
    return RenormedStaticFrequencyTable<source_type>(std::move(rescaledFrequencies), offset, newPrecision, incompressibleSymbolFrequency);
  } else if constexpr (std::is_same_v<frequencyTable_T, DynamicFrequencyTable<source_type>>) {
    return RenormedDynamicFrequencyTable<source_type>(std::move(rescaledFrequencies), offset, newPrecision, incompressibleSymbolFrequency);
  }
}

} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_RENORM_H_ */
