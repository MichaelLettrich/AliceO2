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
/// @since  2019-05-10
/// @brief  class for encoding symbols using rANS

#ifndef RANS_INTERNAL_CPP_SIMPLEENCODER_H
#define RANS_INTERNAL_CPP_SIMPLEENCODER_H

#include <vector>
#include <cstdint>
#include <cassert>
#include <type_traits>
#include <tuple>

#include "rANS/internal/backend/cpp/EncoderSymbol.h"
#include "rANS/internal/helper.h"

namespace o2
{
namespace rans
{
namespace internal
{
namespace cpp
{

template <typename state_T, typename stream_T, uint8_t lowerBound_V>
class SimpleEncoder
{

  static_assert(sizeof(state_T) > sizeof(stream_T), "We cannot stream out more than the size of our state");
  __extension__ using uint128_t = unsigned __int128;

 public:
  explicit SimpleEncoder(size_t symbolTablePrecission) noexcept;

  template <typename stream_IT>
  stream_IT flush(stream_IT outputIter);

  // Encodes a given symbol.
  template <typename stream_IT>
  stream_IT putSymbol(stream_IT outputIter, const EncoderSymbol<state_T>& symbol);

 private:
  state_T mState{LowerBound};
  size_t mSymbolTablePrecission{};

  // Renormalize the encoder.
  template <typename stream_IT>
  std::tuple<state_T, stream_IT> renorm(state_T state, stream_IT outputIter, uint32_t frequency);

  // L ('l' in the paper) is the lower bound of our normalization interval.
  // Between this and our byte-aligned emission, we use 31 (not 32!) bits.
  // This is done intentionally because exact reciprocals for 31-bit uints
  // fit in 32-bit uints: this permits some optimizations during encoding.
  inline static constexpr state_T LowerBound = pow2(lowerBound_V);

  inline static constexpr state_T StreamBits = toBits(sizeof(stream_T));
};

template <typename state_T, typename stream_T, uint8_t lowerBound_V>
SimpleEncoder<state_T, stream_T, lowerBound_V>::SimpleEncoder(size_t symbolTablePrecission) noexcept : mSymbolTablePrecission(symbolTablePrecission){};

template <typename state_T, typename stream_T, uint8_t lowerBound_V>
template <typename stream_IT>
stream_IT SimpleEncoder<state_T, stream_T, lowerBound_V>::flush(stream_IT outputIter)
{
  constexpr size_t StateBits = toBits(sizeof(state_T));

  stream_IT streamPosition = outputIter;
  for (size_t shift = StateBits - StreamBits; shift > 0; shift -= StreamBits) {
    *(++streamPosition) = static_cast<stream_T>(mState >> shift);
  }
  *(++streamPosition) = static_cast<stream_T>(mState >> 0);

  mState = 0;
  return streamPosition;
};

template <typename state_T, typename stream_T, uint8_t lowerBound_V>
template <typename stream_IT>
stream_IT SimpleEncoder<state_T, stream_T, lowerBound_V>::putSymbol(stream_IT outputIter, const EncoderSymbol<state_T>& symbol)
{

  assert(symbol.getFrequency() != 0); // can't encode symbol with freq=0

  // renormalize
  const auto [newState, streamPosition] = renorm(mState, outputIter, symbol.getFrequency());

  // x = C(s,x)
  state_T quotient = 0;

  if constexpr (needs64Bit<state_T>()) {
    // This code needs support for 64-bit long multiplies with 128-bit result
    // (or more precisely, the top 64 bits of a 128-bit result).
    quotient = static_cast<state_T>((static_cast<uint128_t>(newState) * symbol.getReciprocalFrequency()) >> 64);
  } else {
    quotient = static_cast<state_T>((static_cast<uint64_t>(newState) * symbol.getReciprocalFrequency()) >> 32);
  }
  quotient = quotient >> symbol.getReciprocalShift();

  mState = newState + symbol.getBias() + quotient * symbol.getFrequencyComplement();
  return streamPosition;
};

template <typename state_T, typename stream_T, uint8_t lowerBound_V>
template <typename stream_IT>
inline std::tuple<state_T, stream_IT> SimpleEncoder<state_T, stream_T, lowerBound_V>::renorm(state_T state, stream_IT outputIter, uint32_t frequency)
{
  state_T maxState = ((LowerBound >> mSymbolTablePrecission) << StreamBits) * frequency; // this turns into a shift.
  if (state >= maxState) {
    *(++outputIter) = static_cast<stream_T>(state);
    state >>= StreamBits;
    assert(state < maxState);
  }

  return std::make_tuple(state, outputIter);
}; // namespace cpp

} // namespace cpp
} // namespace internal
} // namespace rans
} // namespace o2

#endif /* RANS_INTERNAL_CPP_SIMPLEENCODER_H */

/*
// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   SimpleEncoder.h
/// @author Michael Lettrich
/// @since  2019-05-10
/// @brief  class for encoding symbols using rANS

#ifndef RANS_INTERNAL_CPP_SIMPLEENCODER_H
#define RANS_INTERNAL_CPP_SIMPLEENCODER_H

#include <vector>
#include <cstdint>
#include <cassert>
#include <type_traits>
#include <tuple>

#include "rANS/internal/backend/cpp/EncoderSymbol.h"
#include "rANS/internal/helper.h"

namespace o2
{
namespace rans
{
namespace internal
{
namespace cpp
{

template <typename state_T, typename stream_T, uint8_t lowerBound_V>
class SimpleEncoder
{
  static_assert(sizeof(state_T) > sizeof(stream_T), "We cannot stream out more than the size of our state");

 public:
  explicit SimpleEncoder(size_t symbolTablePrecission) noexcept;

  template <typename stream_IT>
  stream_IT flush(stream_IT outputIter);

  // Encodes a given symbol.
  template <typename stream_IT>
  stream_IT putSymbol(stream_IT outputIter, const EncoderSymbol<state_T>& symbol);

 private:
  state_T mState{pow2(lowerBound_V)};
  size_t mSymbolTablePrecision{};

  // Renormalize the encoder.
  template <typename stream_IT>
  std::tuple<state_T, stream_IT> renorm(state_T state, stream_IT outputIter, uint32_t frequency);

  // L ('l' in the paper) is the lower bound of our normalization interval.
  // Between this and our byte-aligned emission, we use 31 (not 32!) bits.
  // This is done intentionally because exact reciprocals for 31-bit uints
  // fit in 32-bit uints: this permits some optimizations during encoding.
  //inline static constexpr state_T LOWER_BOUND = static_cast<state_T>(1ull << lowerBound_V); // lower bound of our normalization interval
};

template <typename state_T, typename stream_T, uint8_t lowerBound_V>
SimpleEncoder<state_T, stream_T, lowerBound_V>::SimpleEncoder(size_t symbolTablePrecission) noexcept : mSymbolTablePrecision(symbolTablePrecission){};

template <typename state_T, typename stream_T, uint8_t lowerBound_V>
template <typename stream_IT>
stream_IT SimpleEncoder<state_T, stream_T, lowerBound_V>::flush(stream_IT outputIter)
{
  constexpr size_t StateBits = toBits(sizeof(state_T));
  constexpr size_t StreamBits = toBits(sizeof(stream_T));

  stream_IT streamPosition = outputIter;
  for (size_t shift = StateBits - StreamBits; shift > 0; shift -= StreamBits) {
    *(++streamPosition) = static_cast<stream_T>(mState >> shift);
  }
  *(++streamPosition) = static_cast<stream_T>(mState >> 0);

  mState = 0;
  return streamPosition;
};

template <typename state_T, typename stream_T, uint8_t lowerBound_V>
template <typename stream_IT>
stream_IT SimpleEncoder<state_T, stream_T, lowerBound_V>::putSymbol(stream_IT outputIter, const EncoderSymbol<state_T>& symbol)
{
  assert(symbol.getFrequency() != 0); // can't encode symbol with freq=0
  // renormalize
  const auto [newState, streamPosition] = renorm(mState, outputIter, symbol.getFrequency());
  // x = C(s,x)
  mState = ((newState / symbol.getFrequency()) << mSymbolTablePrecision) + symbol.getBias() + (newState % symbol.getFrequency());
  return streamPosition;
};

template <typename state_T, typename stream_T, uint8_t lowerBound_V>
template <typename stream_IT>
inline std::tuple<state_T, stream_IT> SimpleEncoder<state_T, stream_T, lowerBound_V>::renorm(state_T state, stream_IT outputIter, uint32_t frequency)
{
  constexpr size_t StreamBits = toBits(sizeof(stream_T));
  constexpr size_t LowerBound = pow2(lowerBound_V);

  state_T maxState = ((LowerBound >> mSymbolTablePrecision) << StreamBits) * frequency; // this turns into a shift.
  if (state >= maxState) {
    ++outputIter;
    *outputIter = static_cast<stream_T>(state);
    state >>= StreamBits;
    assert(state < maxState);
  }

  return std::make_tuple(state, outputIter);
}; // namespace cpp

} // namespace cpp
} // namespace internal
} // namespace rans
} // namespace o2

#endif /* RANS_INTERNAL_CPP_SIMPLEENCODER_H */