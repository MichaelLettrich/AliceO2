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

#include "rANS/internal/FrequencyContainer.h"

namespace o2
{
namespace rans
{

template <typename source_T>
class RenormedFrequencyTable : public internal::FrequencyContainer<source_T>
{
  using base_type = internal::FrequencyContainer<source_T>;

 public:
  using source_type = typename base_type::source_type;
  using value_type = typename base_type::value_type;
  using container_type = typename base_type::container_type;
  using size_type = typename base_type::size_type;
  using difference_type = typename base_type::difference_type;
  using reference = typename base_type::reference;
  using const_reference = typename base_type::const_reference;
  using pointer = typename base_type::pointer;
  using const_pointer = typename base_type::const_pointer;
  using const_iterator = typename base_type::const_iterator;

  RenormedFrequencyTable() : base_type(){};

  inline RenormedFrequencyTable(container_type frequencies, source_type offset, size_t renormingBits, value_type nIncompressible) : mNIncompressible(nIncompressible)
  {
    this->mContainer = std::move(frequencies);
    this->setOffset(offset);
    this->mNSamples = internal::pow2(renormingBits);

    // TODO(milettri): do some checks when nDebug is active;
  };

  [[nodiscard]] inline size_t getRenormingBits() const noexcept { return internal::log2UInt(this->mNSamples); };

  [[nodiscard]] inline bool isRenormedTo(size_t nBits) const noexcept { return nBits == this->getRenormingBits(); };

  [[nodiscard]] inline value_type getIncompressibleSymbolFrequency() const noexcept { return mNIncompressible; };

  [[nodiscard]] inline bool hasIncompressibleSymbol() const noexcept { return mNIncompressible != 0; };

  friend void swap(RenormedFrequencyTable& a, RenormedFrequencyTable& b) noexcept
  {
    using std::swap;
    swap(static_cast<typename RenormedFrequencyTable::base_type&>(a),
         static_cast<typename RenormedFrequencyTable::base_type&>(b));
  };

 private:
  value_type mNIncompressible{};
};

} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_RENORMEDFREQUENCYTABLEIMPL_H_ */
