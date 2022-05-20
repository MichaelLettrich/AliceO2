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
#include "rANS/HashFrequencyTable.h"
#include "rANS/FrequencyTable.h"
#include "rANS/RenormedFrequencyTable.h"

namespace o2
{
namespace rans
{

template <typename source_T>
inline constexpr size_t computeRenormingPrecision(size_t nUsedAlphabetSymbols)
{
  if (nUsedAlphabetSymbols == 0) {
    return 0;
  }

  if constexpr (sizeof(source_T) == 1) {
    return 14;
  }
  if constexpr (sizeof(source_T) == 2) {
    return 18;
  } else {
    const uint8_t estimate = internal::log2UInt(nUsedAlphabetSymbols) * 3u / 2u;
    const uint8_t upperBound = std::min(estimate, MaxRenormThreshold);
    const uint8_t lowerBound = std::max(estimate, MinRenormThreshold);
    return std::min(lowerBound, upperBound);
  };
}

namespace renorm_impl
{

template <typename frequencyTable_T>
inline size_t getIncompressibleSymbolFrequency(const frequencyTable_T&)
{
  return 0;
}

inline size_t getIncompressibleSymbolFrequency(const FrequencyTable& f)
{
  return f.getIncompressibleSymbolFrequency();
}

inline size_t getNUsedAlphabetSymbols(const FrequencyTable& f)
{
  return f.getNUsedAlphabetSymbols();
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

inline auto getOffset(const FrequencyTable& f)
{
  return f.getMinSymbol();
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

  auto expandSymbol = [&](auto iter) -> std::pair<index_type, count_type> {
    if constexpr (std::is_same_v<container_type, typename HashFrequencyTable<source_type>::container_type>) {
      return *iter;
    } else {
      return {static_cast<index_type>(rescaledFrequencies.begin, iter), *iter};
    }
  };

  auto getFrequency = [&](auto iter) -> count_type {
    if constexpr (std::is_same_v<container_type, typename HashFrequencyTable<source_type>::container_type>) {
      return iter->second;
    } else {
      return *iter;
    }
  };

  auto setFrequency = [&](auto iter, count_type value) -> void {
    if constexpr (std::is_same_v<container_type, typename HashFrequencyTable<source_type>::container_type>) {
      iter->second = value;
    } else {
      *iter = value;
    }
  };

  auto setToZero = [&](auto iter) -> void {
    if constexpr (std::is_same_v<container_type, typename HashFrequencyTable<source_type>::container_type>) {
      rescaledFrequencies.erase(iter);
    } else {
      *iter = 0;
    }
  };

  for (auto frequencyIter = rescaledFrequencies.begin(); frequencyIter != rescaledFrequencies.end(); ++frequencyIter) {
    const count_type frequency = getFrequency(frequencyIter);
    if (frequency > 0) {
      const double_t symbolProbability = static_cast<double_t>(frequency) / nSamples;
      if (symbolProbability < probabilityCutOffThreshold) {
        incompressibleSymbolProbability += symbolProbability;
        setToZero(frequencyIter);
      } else {
        const double_t scaledFrequencyD = scaleFrequency(symbolProbability);
        count_type rescaledFrequency = roundFrequency(scaledFrequencyD);
        assert(rescaledFrequency > 0);
        setFrequency(frequencyIter, rescaledFrequency);
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
  std::stable_sort(correctableIndices.begin(), correctableIndices.end(), [&rescaledFrequencies](const iterator_type& a, const iterator_type& b) {
    if constexpr (std::is_same_v<container_type, typename HashFrequencyTable<source_type>::container_type>) {
      // iterator over unodered_map is by storage location and not key. For same (stable) results of array and hash map approach, we need to enforce the same order.
      if (a->second == b->second) {
        return a->first < b->first;
      } else {
        return a->second < b->second;
      }
    } else {
      return *a < *b;
    }
  });

  difference_type nCorrections = static_cast<difference_type>(nSamplesRescaled) - static_cast<difference_type>(nSamplesRescaledUncorrected);
  const double_t rescalingFactor = static_cast<double_t>(nSamplesRescaled) / static_cast<double_t>(nSamplesRescaledUncorrected);

  for (auto iter : correctableIndices) {
    if (std::abs(nCorrections) > 0) {
      const difference_type uncorrectedFrequency = getFrequency(iter);
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
      setFrequency(iter, correctedFrequency);
    } else {
      break;
    }
  }

  if (std::abs(nCorrections) > 0) {
    throw std::runtime_error(fmt::format("rANS rescaling incomplete: {} corrections Remaining", nCorrections));
  }

  if constexpr (std::is_same_v<frequencyTable_T, FrequencyTable>) {
    frequencyTable_T newFrequencyTable{rescaledFrequencies.begin(), rescaledFrequencies.end(), offset, incompressibleSymbolFrequency};
    return RenormedFrequencyTable{std::move(newFrequencyTable), newPrecision};
  } else if constexpr (std::is_same_v<frequencyTable_T, StaticFrequencyTable<source_type>>) {
    return RenormedStaticFrequencyTable<source_type>(std::move(rescaledFrequencies), offset, newPrecision, incompressibleSymbolFrequency);
  } else if constexpr (std::is_same_v<frequencyTable_T, DynamicFrequencyTable<source_type>>) {
    return RenormedDynamicFrequencyTable<source_type>(std::move(rescaledFrequencies), offset, newPrecision, incompressibleSymbolFrequency);
  } else if constexpr (std::is_same_v<frequencyTable_T, HashFrequencyTable<source_type>>) {
    return RenormedHashFrequencyTable<source_type>(std::move(rescaledFrequencies), offset, newPrecision, incompressibleSymbolFrequency);
  }
}

} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_RENORM_H_ */
