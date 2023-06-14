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

/// @file   DecoderSymbol.h
/// @author Michael Lettrich
/// @since  2019-05-21
/// @brief  Structure containing all relevant information for decoding a rANS encoded symbol

#ifndef RANSLEGACY_INTERNAL_CPP_DECODERSYMBOL_H
#define RANSLEGACY_INTERNAL_CPP_DECODERSYMBOL_H

#include <cstdint>
#include <cstring>
#include <cassert>

#include "rANSLegacy/definitions.h"

namespace o2
{
namespace ranslegacy
{
namespace internal
{
namespace cpp
{

// Decoder symbols are straightforward.
class DecoderSymbol
{
 public:
  using value_type = count_t;

  //TODO(milettri): fix once ROOT cling respects the standard http://wg21.link/p1286r2
  constexpr DecoderSymbol() noexcept {}; //NOLINT
  // Initialize a decoder symbol to start "start" and frequency "freq"
  constexpr DecoderSymbol(value_type frequency, value_type cumulative, size_t symbolTablePrecision)
    : mCumulative(cumulative), mFrequency(frequency)
  {
    (void)symbolTablePrecision; // silence compiler warnings if assert not compiled.
    assert(mCumulative <= pow2(symbolTablePrecision));
    assert(mFrequency <= pow2(symbolTablePrecision) - mCumulative);
  };
  inline constexpr value_type getFrequency() const noexcept { return mFrequency; };
  inline constexpr value_type getCumulative() const noexcept { return mCumulative; };

  inline bool operator==(const DecoderSymbol& other) const { return this->mCumulative == other.mCumulative; };

 private:
  value_type mCumulative{}; // Start of range.
  value_type mFrequency{};  // Symbol frequency.
};

} // namespace cpp
} // namespace internal
} // namespace ranslegacy
} // namespace o2

#endif /* RANSLEGACY_INTERNAL_CPP_DECODERSYMBOL_H */
