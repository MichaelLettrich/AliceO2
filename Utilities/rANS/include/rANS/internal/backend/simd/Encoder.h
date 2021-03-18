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

#ifndef RANS_INTERNAL_SIMD__ENCODER_H
#define RANS_INTERNAL_SIMD__ENCODER_H

#include <vector>
#include <cstdint>
#include <cassert>
#include <type_traits>
#include <tuple>

#include <fairlogger/Logger.h>

#include "rANS/internal/backend/simd/EncoderSymbol.h"
#include "rANS/internal/backend/simd/utils.h"
#include "rANS/internal/helper.h"

namespace o2
{
namespace rans
{
namespace internal
{
namespace simd
{

template <typename State_T, typename Stream_T, size_t nHardwareStreams_V>
class Encoder
{
  // the Coder works either with a 64Bit state and 32 bit streaming or
  //a 32 Bit state and 8 Bit streaming We need to make sure it gets initialized with
  //the right template arguments at compile time.
  static_assert((sizeof(State_T) == sizeof(uint64_t) && sizeof(Stream_T) == sizeof(uint32_t)),
                "Coder can either be 32Bit with 8 Bit stream type or 64 Bit Type with 32 Bit stream type");

 public:
  Encoder();

  // Flushes the rANS encoder.
  template <typename Stream_IT>
  Stream_IT flush(Stream_IT iter);

  template <typename Stream_IT>
  Stream_IT putSymbols(Stream_IT iter, const std::array<EncoderSymbol, nHardwareStreams_V>& sym, uint32_t scale_bits);

  template <typename Stream_IT>
  Stream_IT putSymbols(Stream_IT iter, const std::array<EncoderSymbol, nHardwareStreams_V>& sym, uint32_t scale_bits, size_t mask);

 private:
  std::array<State_T, nHardwareStreams_V> mStates;

  template <typename Stream_IT>
  Stream_IT putSymbol(Stream_IT iter, const EncoderSymbol& sym, uint32_t scale_bits, State_T& state);

  template <typename Stream_IT>
  Stream_IT flushState(State_T& state, Stream_IT iter);

  // Renormalize the encoder.
  template <typename Stream_IT>
  std::tuple<State_T, Stream_IT> renorm(State_T x, Stream_IT iter, uint32_t freq, uint32_t scale_bits);

  // L ('l' in the paper) is the lower bound of our normalization interval.
  // Between this and our byte-aligned emission, we use 31 (not 32!) bits.
  // This is done intentionally because exact reciprocals for 31-bit uints
  // fit in 32-bit uints: this permits some optimizations during encoding.
  inline static constexpr State_T LOWER_BOUND = (1u << 31); // lower bound of our normalization interval

  inline static constexpr State_T STREAM_BITS = sizeof(Stream_T) * 8; // lower bound of our normalization interval
};

template <typename State_T, typename Stream_T, size_t nHardwareStreams_V>
Encoder<State_T, Stream_T, nHardwareStreams_V>::Encoder()
{
  LOG(INFO) << "lower bound" << LOWER_BOUND;
  for (auto& state : mStates) {
    state = LOWER_BOUND;
  }
};

template <typename State_T, typename Stream_T, size_t nHardwareStreams_V>
template <typename Stream_IT>
Stream_IT Encoder<State_T, Stream_T, nHardwareStreams_V>::flush(Stream_IT iter)
{
  Stream_IT streamPos = iter;
  for (auto stateIter = mStates.rbegin(); stateIter != mStates.rend(); ++stateIter) {
    streamPos = flushState(*stateIter, streamPos);
  }
  return streamPos;
};

template <typename State_T, typename Stream_T, size_t nHardwareStreams_V>
template <typename Stream_IT>
Stream_IT Encoder<State_T, Stream_T, nHardwareStreams_V>::putSymbols(Stream_IT iter, const std::array<EncoderSymbol, nHardwareStreams_V>& symbols, uint32_t scale_bits)
{

  // can't encode symbol with freq=0
  for (const auto& symbol : symbols) {
    assert(symbol.getFrequency() != 0);
  }

  // normalize in reverse direction for decoding to work propperly
  Stream_IT streamPosition = [this, &symbols, iter, scale_bits]() {
    auto streamIter = iter;
    for (size_t i = nHardwareStreams_V; i-- > 0;) {
      auto [tmpState, tmpStreamIter] = renorm(mStates[i], streamIter, symbols[i].getFrequency(), scale_bits);
      mStates[i] = tmpState;
      streamIter = tmpStreamIter;
    }
    return streamIter;
  }();

  //calculate div and mod
  auto [div, mod] = [this, &symbols]() {
    std::array<uint64_t, nHardwareStreams_V> div;
    std::array<uint64_t, nHardwareStreams_V> mod;

    for (size_t i = 0; i < nHardwareStreams_V; ++i) {
      div[i] = mStates[i] / symbols[i].getFrequency();
      mod[i] = mStates[i] - div[i] * symbols[i].getFrequency();
    }
    return std::make_tuple(div, mod);
  }();

  // encode
  // x = C(s,x)
  [this, &symbols, scale_bits](auto& div, auto& mod) {
    for (size_t i = 0; i < nHardwareStreams_V; ++i) {
      mStates[i] = (div[i] << scale_bits) + mod[i] + symbols[i].getCumulative();
    }
  }(div, mod);
  return streamPosition;
} // namespace fp64

template <typename State_T, typename Stream_T, size_t nHardwareStreams_V>
template <typename Stream_IT>
Stream_IT Encoder<State_T, Stream_T, nHardwareStreams_V>::putSymbols(Stream_IT iter, const std::array<EncoderSymbol, nHardwareStreams_V>& symbols, uint32_t scale_bits, size_t nActiveStreams)
{
  Stream_IT streamPos = iter;

  for (size_t i = nActiveStreams; i-- > 0;) {
    streamPos = putSymbol(streamPos, symbols[i], scale_bits, mStates[i]);
  }

  return streamPos;
};

template <typename State_T, typename Stream_T, size_t nHardwareStreams_V>
template <typename Stream_IT>
Stream_IT Encoder<State_T, Stream_T, nHardwareStreams_V>::putSymbol(Stream_IT iter, const EncoderSymbol& sym, uint32_t scale_bits, State_T& state)
{
  assert(sym.getFrequency() != 0); // can't encode symbol with freq=0
  // renormalize
  const auto [x, streamPos] = renorm(state, iter, sym.getFrequency(), scale_bits);

  // x = C(s,x)
  state = ((x / sym.getFrequency()) << scale_bits) + (x % sym.getFrequency()) + sym.getCumulative();
  return streamPos;
}

template <typename State_T, typename Stream_T, size_t nHardwareStreams_V>
template <typename Stream_IT>
Stream_IT Encoder<State_T, Stream_T, nHardwareStreams_V>::flushState(State_T& state, Stream_IT iter)
{
  Stream_IT streamPos = iter;

  ++streamPos;
  *streamPos = static_cast<Stream_T>(state >> 32);
  ++streamPos;
  *streamPos = static_cast<Stream_T>(state >> 0);

  state = 0;
  return streamPos;
}

template <typename State_T, typename Stream_T, size_t nHardwareStreams_V>
template <typename Stream_IT>
inline std::tuple<State_T, Stream_IT> Encoder<State_T, Stream_T, nHardwareStreams_V>::renorm(State_T x, Stream_IT iter, uint32_t freq, uint32_t scale_bits)
{
  Stream_IT streamPos = iter;

  State_T x_max = ((LOWER_BOUND >> scale_bits) << STREAM_BITS) * freq; // this turns into a shift.
  if (x >= x_max) {
    ++streamPos;
    *streamPos = static_cast<Stream_T>(x);
    x >>= STREAM_BITS;
    assert(x < x_max);
  }
  return std::make_tuple(x, streamPos);
};

} // namespace simd
} // namespace internal
} // namespace rans
} // namespace o2

#endif /* RANS_INTERNAL_SIMD__ENCODER_H */
