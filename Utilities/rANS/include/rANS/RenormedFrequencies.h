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

#ifndef INCLUDE_RANS_RENORMEDFREQUENCYTABLEIMPL_H_
#define INCLUDE_RANS_RENORMEDFREQUENCYTABLEIMPL_H_

#include <fairlogger/Logger.h>

#include "rANS/internal/StaticFrequencyContainer.h"
#include "rANS/internal/DynamicFrequencyContainer.h"
#include "rANS/internal/HashFrequencyContainer.h"

namespace o2
{
namespace rans
{

template <class frequencyContainer_T>
class RenormedFrequencyTable_Impl : public frequencyContainer_T
{
  using base_type = frequencyContainer_T;

 public:
  using source_type = typename base_type::source_type;
  using index_type = typename base_type::index_type;
  using value_type = typename base_type::value_type;
  using container_type = typename base_type::container_type;
  using size_type = typename base_type::size_type;
  using difference_type = typename base_type::difference_type;
  using reference = typename base_type::reference;
  using const_reference = typename base_type::const_reference;
  using pointer = typename base_type::pointer;
  using const_pointer = typename base_type::const_pointer;
  using const_iterator = typename base_type::const_iterator;

  RenormedFrequencyTable_Impl() : base_type(){};

  inline RenormedFrequencyTable_Impl(container_type frequencies, source_type offset, size_t renormingBits, value_type nIncompressible) : mNIncompressible(nIncompressible)
  {
    this->mContainer = std::move(frequencies);
    this->mOffset = offset;
    this->mNSamples = internal::pow2(renormingBits);

    // TODO(milettri): do some checks when nDebug is active;
  };

  [[nodiscard]] inline size_t getRenormingBits() const noexcept { return internal::log2UInt(this->mNSamples); };

  [[nodiscard]] inline bool isRenormedTo(size_t nBits) const noexcept { return nBits == this->getRenormingBits(); };

  [[nodiscard]] inline value_type getIncompressibleSymbolFrequency() const noexcept { return mNIncompressible; };

  [[nodiscard]] inline bool hasIncompressibleSymbol() const noexcept { return mNIncompressible != 0; };

 private:
  value_type mNIncompressible{};
};

template <typename source_T>
using RenormedStaticFrequencyTable = RenormedFrequencyTable_Impl<StaticFrequencyContainer<source_T>>;
template <typename source_T>
using RenormedDynamicFrequencyTable = RenormedFrequencyTable_Impl<DynamicFrequencyContainer<source_T>>;
template <typename source_T>
using RenormedHashFrequencyTable = RenormedFrequencyTable_Impl<HashFrequencyContainer<source_T>>;

} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_RENORMEDFREQUENCYTABLEIMPL_H_ */
