// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   Encoder.h
/// @author Michael Lettrich
/// @since  2021-03-18
/// @brief  class for encoding symbols using rANS

#ifndef RANS_INTERNAL_ENCODECOMMANDINTERFACE_H
#define RANS_INTERNAL_ENCODECOMMANDINTERFACE_H

#include <vector>
#include <cstdint>
#include <cassert>
#include <type_traits>
#include <tuple>

#include <fairlogger/Logger.h>

#include "rANS/internal/helper.h"

namespace o2
{
namespace rans
{
namespace internal
{

template <typename symbol_T, typename derived_T>
class EncodeCommandInterface
{
 public:
  using stream_type = uint32_t;
  using state_type = uint64_t;
  using symbol_type = symbol_T;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;

  [[nodiscard]] inline static constexpr size_type getNstreams() noexcept
  {
    return derived_T::getNstreams();
  };

  // Flushes the rANS encoder.
  template <typename Stream_IT>
  [[nodiscard]] inline Stream_IT flush(Stream_IT outputIter)
  {
    return static_cast<derived_T*>(this)->flush(outputIter);
  };

  template <typename Stream_IT>
  [[nodiscard]] inline Stream_IT putSymbols(Stream_IT outputIter, const symbol_type& encodeSymbols)
  {
    return static_cast<derived_T*>(this)->putSymbols(outputIter, encodeSymbols);
  };

  template <typename Stream_IT>
  [[nodiscard]] inline Stream_IT putSymbols(Stream_IT outputIter, const symbol_type& encodeSymbols, size_type nActiveStreams)
  {
    return static_cast<derived_T*>(this)->putSymbols(outputIter, encodeSymbols, nActiveStreams);
  };

  [[nodiscard]] inline static constexpr state_type getStreamingLowerBound() noexcept
  {
    return derived_T::getStreamingLowerBound();
  };

 protected:
  [[nodiscard]] inline static constexpr state_type getStreamOutTypeBits() noexcept
  {
    return toBits(sizeof(stream_type));
  };

  EncodeCommandInterface() = default;
  explicit EncodeCommandInterface(size_t symbolTablePrecision) noexcept : mSymbolTablePrecision{symbolTablePrecision} {};

  size_type mSymbolTablePrecision{};
};

} // namespace internal
} // namespace rans
} // namespace o2

#endif /* RANS_INTERNAL_ENCODECOMMANDINTERFACE_H */
