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

#ifndef RANS_INTERNAL_SIMD_ENCODER_H
#define RANS_INTERNAL_SIMD_ENCODER_H

#include <vector>
#include <cstdint>
#include <cassert>
#include <type_traits>
#include <tuple>

#include <fairlogger/Logger.h>

#ifdef ENABLE_VTUNE_PROFILER
#include <ittnotify.h>
#endif

#include "rANS/internal/backend/simd/Symbol.h"
#include "rANS/internal/backend/simd/kernel.h"
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

template <SIMDWidth simdWidth_V>
class Encoder
{
 public:
  using stream_t = uint32_t;
  using state_t = uint64_t;
  inline static constexpr size_t nHardwareStreams = getElementCount<state_t>(simdWidth_V);

  Encoder(size_t symbolTablePrecision);
  Encoder() : Encoder{0} {};

  // Flushes the rANS encoder.
  template <typename Stream_IT>
  Stream_IT flush(Stream_IT outputIter);

  template <typename Stream_IT>
  Stream_IT putSymbols(Stream_IT outputIter, ArrayView<const Symbol*, nHardwareStreams> encodeSymbols);

  template <typename Stream_IT>
  Stream_IT putSymbols(Stream_IT outputIter, ArrayView<const Symbol*, nHardwareStreams> encodeSymbols, size_t mask);

 private:
  epi64_t<simdWidth_V> mStates;
  size_t mSymbolTablePrecision{};
  pd_t<simdWidth_V> mNSamples{};

  template <typename Stream_IT>
  Stream_IT putSymbol(Stream_IT outputIter, const Symbol& symbol, state_t& state);

  template <typename Stream_IT>
  Stream_IT flushState(state_t& state, Stream_IT outputIter);

  // Renormalize the encoder.
  template <typename Stream_IT>
  std::tuple<state_t, Stream_IT> renorm(state_t state, Stream_IT outputIter, uint32_t frequency);

  // L ('l' in the paper) is the lower bound of our normalization interval.
  // Between this and our byte-aligned emission, we use 31 (not 32!) bits.
  // This is done intentionally because exact reciprocals for 31-bit uints
  // fit in 32-bit uints: this permits some optimizations during encoding.
  inline static constexpr state_t LOWER_BOUND = pow2(20); // lower bound of our normalization interval

  inline static constexpr state_t STREAM_BITS = toBits(sizeof(stream_t)); // lower bound of our normalization interval
};

template <SIMDWidth simdWidth_V>
Encoder<simdWidth_V>::Encoder(size_t symbolTablePrecision) : mStates{LOWER_BOUND}, mSymbolTablePrecision{symbolTablePrecision}, mNSamples{static_cast<double>(pow2(mSymbolTablePrecision))}
{
  if (mSymbolTablePrecision > LOWER_BOUND) {
    throw std::runtime_error(fmt::format("[{}]: SymbolTable Precision of {} Bits is larger than allowed by the rANS Encoder (max {} Bits)", __PRETTY_FUNCTION__, mSymbolTablePrecision, LOWER_BOUND));
  }
};

template <SIMDWidth simdWidth_V>
template <typename Stream_IT>
Stream_IT Encoder<simdWidth_V>::flush(Stream_IT iter)
{
  Stream_IT streamPos = iter;
  for (auto stateIter = mStates.rbegin(); stateIter != mStates.rend(); ++stateIter) {
    streamPos = flushState(*stateIter, streamPos);
  }
  return streamPos;
};

template <SIMDWidth simdWidth_V>
template <typename Stream_IT>
Stream_IT Encoder<simdWidth_V>::putSymbols(Stream_IT outputIter, ArrayView<const Symbol*, nHardwareStreams> encodeSymbols)
{

  // can't encode symbol with freq=0
  for (const auto& symbol : encodeSymbols) {
    assert(symbol.getFrequency() != 0);
  }

  auto [frequencies, cumulativeFrequencies] = aosToSoa(encodeSymbols);

  const auto [streamPosition, newStates] = ransRenorm<Stream_IT, LOWER_BOUND, STREAM_BITS>(makeCSView(mStates), makeCSView(frequencies), static_cast<uint8_t>(mSymbolTablePrecision), outputIter);
  mStates = ransEncode(makeCSView(newStates), int32ToDouble<simdWidth_V>(makeCSView(frequencies)), int32ToDouble<simdWidth_V>(makeCSView(cumulativeFrequencies)), mNSamples);
  return streamPosition;
}

template <SIMDWidth simdWidth_V>
template <typename Stream_IT>
Stream_IT Encoder<simdWidth_V>::putSymbols(Stream_IT outputIter, ArrayView<const Symbol*, nHardwareStreams> encodeSymbols, size_t mask)
{
  Stream_IT streamPos = outputIter;

  for (size_t i = mask; i-- > 0;) {
    streamPos = putSymbol(streamPos, *encodeSymbols[i], mStates[i]);
  }

  return streamPos;
};

template <SIMDWidth simdWidth_V>
template <typename Stream_IT>
Stream_IT Encoder<simdWidth_V>::putSymbol(Stream_IT outputIter, const Symbol& symbol, state_t& state)
{
  assert(symbol.getFrequency() != 0); // can't encode symbol with freq=0
  // renormalize
  const auto [x, streamPos] = renorm(state, outputIter, symbol.getFrequency());

  // x = C(s,x)
  state = ((x / symbol.getFrequency()) << mSymbolTablePrecision) + (x % symbol.getFrequency()) + symbol.getCumulative();
  return streamPos;
}

template <SIMDWidth simdWidth_V>
template <typename Stream_IT>
Stream_IT Encoder<simdWidth_V>::flushState(state_t& state, Stream_IT iter)
{
  Stream_IT streamPosition = iter;

  ++streamPosition;
  *streamPosition = static_cast<stream_t>(state >> 32);
  ++streamPosition;
  *streamPosition = static_cast<stream_t>(state >> 0);

  state = 0;
  return streamPosition;
}

template <SIMDWidth simdWidth_V>
template <typename Stream_IT>
inline auto Encoder<simdWidth_V>::renorm(state_t state, Stream_IT outputIter, uint32_t frequency) -> std::tuple<state_t, Stream_IT>
{
  state_t maxState = ((LOWER_BOUND >> mSymbolTablePrecision) << STREAM_BITS) * frequency; // this turns into a shift.
  if (state >= maxState) {
    ++outputIter;
    *outputIter = static_cast<stream_t>(state);
    state >>= STREAM_BITS;
    assert(state < maxState);
  }
  return std::make_tuple(state, outputIter);
};

} // namespace simd
} // namespace internal
} // namespace rans
} // namespace o2

#endif /* RANS_INTERNAL_SIMD_ENCODER_H */
