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

/// @file   ReverseSymbolLookupTable.h
/// @author Michael Lettrich
/// @since  2020-04-06
/// @brief  Maps CDF back to source symbol - needed for the decoder

#ifndef RANS_INTERNAL_REVERSESYMBOLLOOKUPTABLE_H
#define RANS_INTERNAL_REVERSESYMBOLLOOKUPTABLE_H

#include <vector>
#include <type_traits>
#include <fairlogger/Logger.h>

#include "rANS/definitions.h"
#include "rANS/internal/helper.h"
#include "rANS/RenormedFrequencies.h"

namespace o2
{
namespace rans
{
namespace internal
{

template <typename source_T>
class RLUT
{
 public:
  using source_type = source_T;
  using index_type = source_type;
  using count_type = uint32_t;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using container_type = std::vector<source_type>;
  using iterator_type = source_type*;

  // TODO(milettri): fix once ROOT cling respects the standard http://wg21.link/p1286r2
  inline RLUT() noexcept {}; // NOLINT

  template <typename renormedFrequencyTable_T>
  explicit RLUT(const renormedFrequencyTable_T& frequencyTable)
  {
    if (frequencyTable.empty()) {
      LOG(warning) << "SymbolStatistics of empty message passed to " << __func__;
    }

    mLut.reserve(frequencyTable.getNumSamples());

    index_type symbol = frequencyTable.getOffset();
    for (count_type symbolFrequency : frequencyTable) {
      mLut.insert(mLut.end(), symbolFrequency, symbol++);
    }
  };

  inline size_type size() const noexcept { return mLut.size(); };

  inline bool isIncompressible(count_type cumul) const noexcept
  {
    return cumul >= this->size();
  };

  inline source_type operator[](count_type cumul) const noexcept
  {
    assert(cumul < this->size());
    return mLut[cumul];
  };

  inline const iterator_type begin() const noexcept { return mLut.data(); };
  inline const iterator_type end() const noexcept { return mLut.data() + size(); };

  container_type mLut{};
  count_type mIncompressibleFrequency{};
};

} // namespace internal
} // namespace rans
} // namespace o2

#endif /* RANS_INTERNAL_REVERSESYMBOLLOOKUPTABLE_H */
