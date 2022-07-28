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

/// @file   DecoderBase.h
/// @author michael.lettrich@cern.ch
/// @since  Feb 8, 2021
/// @brief

#ifndef INCLUDE_RANSLEGACY_INTERNAL_DECODERBASE_H_
#define INCLUDE_RANSLEGACY_INTERNAL_DECODERBASE_H_

#include <cstddef>
#include <type_traits>
#include <iostream>
#include <memory>

#include <fairlogger/Logger.h>

#include "rANSLegacy/RenormedFrequencyTable.h"
#include "rANSLegacy/internal/backend/cpp/DecoderSymbol.h"
#include "rANSLegacy/internal/ReverseSymbolLookupTable.h"
#include "rANSLegacy/internal/SymbolTable.h"
#include "rANSLegacy/internal/backend/cpp/Decoder.h"
#include "rANSLegacy/internal/helper.h"

namespace o2
{
namespace ranslegacy
{
namespace internal
{

template <typename coder_T, typename stream_T, typename source_T>
class DecoderBase
{
 protected:
  using decoderSymbolTable_t = SymbolTable<cpp::DecoderSymbol>;
  using reverseSymbolLookupTable_t = ReverseSymbolLookupTable;
  using ransDecoder_t = cpp::Decoder<coder_T, stream_T>;

 public:
  using coder_t = coder_T;
  using stream_t = stream_T;
  using source_t = source_T;

  // TODO(milettri): fix once ROOT cling respects the standard http://wg21.link/p1286r2
  DecoderBase() noexcept {}; // NOLINT
  explicit DecoderBase(const RenormedFrequencyTable& frequencyTable) : mSymbolTable{frequencyTable}, mReverseLUT{frequencyTable} {};

  inline size_t getAlphabetRangeBits() const noexcept { return mSymbolTable.getAlphabetRangeBits(); }
  inline size_t getSymbolTablePrecision() const noexcept { return mSymbolTable.getPrecision(); }
  inline int getMinSymbol() const noexcept { return mSymbolTable.getMinSymbol(); }
  inline int getMaxSymbol() const noexcept { return mSymbolTable.getMaxSymbol(); }

 protected:
  decoderSymbolTable_t mSymbolTable{};
  reverseSymbolLookupTable_t mReverseLUT{};
};
} // namespace internal
} // namespace ranslegacy
} // namespace o2

#endif /* INCLUDE_RANSLEGACY_INTERNAL_DECODERBASE_H_ */
