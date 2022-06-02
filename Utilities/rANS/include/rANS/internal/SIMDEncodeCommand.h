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

#ifndef RANS_INTERNAL_SIMDENCODERS_H
#define RANS_INTERNAL_SIMDENCODERS_H

#include <cassert>
#include <cstdint>
#include <tuple>

#include "rANS/internal/backend/cpp/EncoderSymbol.h"

#include "rANS/internal/EncodeCommandInterface.h"
#include "rANS/internal/helper.h"
#include "rANS/internal/backend/simd/Symbol.h"

#include "rANS/internal/backend/simd/Symbol.h"
#include "rANS/internal/backend/simd/types.h"
#include "rANS/internal/backend/simd/kernel.h"
#include "rANS/internal/backend/simd/utils.h"
#include "rANS/internal/helper.h"

namespace o2
{
namespace rans
{
namespace internal
{
template <size_t streamingLowerBound_V, internal::simd::SIMDWidth simdWidth_V>
class SIMDEncoderCommand : public EncodeCommandInterface<internal::simd::UnrolledSymbols,
                                                         SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>>
{
  using base_type = EncodeCommandInterface<internal::simd::UnrolledSymbols, SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>>;

 public:
  using stream_type = typename base_type::stream_type;
  using state_type = typename base_type::state_type;
  using symbol_type = typename base_type::symbol_type;
  using size_type = typename base_type::size_type;
  using difference_type = typename base_type::difference_type;

  [[nodiscard]] inline static constexpr size_type getNstreams() noexcept { return 2 * internal::simd::getElementCount<state_type>(simdWidth_V); };

  SIMDEncoderCommand(size_t symbolTablePrecision);
  SIMDEncoderCommand() : SIMDEncoderCommand{0} {};

  // Flushes the rANS encoder.
  template <typename Stream_IT>
  Stream_IT flush(Stream_IT outputIter);

  template <typename Stream_IT>
  Stream_IT putSymbols(Stream_IT outputIter, const symbol_type& encodeSymbols);

  template <typename Stream_IT>
  Stream_IT putSymbols(Stream_IT outputIter, const symbol_type& encodeSymbols, size_t nActiveStreams);

 private:
  simd::epi64_t<simdWidth_V, 2> mStates;
  size_t mSymbolTablePrecision{};
  simd::pd_t<simdWidth_V> mNSamples{};

  template <typename Stream_IT>
  Stream_IT putSymbol(Stream_IT outputIter, const simd::Symbol& symbol, state_type& state);

  template <typename Stream_IT>
  Stream_IT flushState(state_type& state, Stream_IT outputIter);

  // Renormalize the encoder.
  template <typename Stream_IT>
  std::tuple<state_type, Stream_IT> renorm(state_type state, Stream_IT outputIter, uint32_t frequency);

  // L ('l' in the paper) is the lower bound of our normalization interval.
  // Between this and our byte-aligned emission, we use 31 (not 32!) bits.
  // This is done intentionally because exact reciprocals for 31-bit uints
  // fit in 32-bit uints: this permits some optimizations during encoding.
  inline static constexpr state_type LowerBound = pow2(20); // lower bound of our normalization interval

  inline static constexpr state_type StreamBits = toBits(sizeof(stream_type)); // lower bound of our normalization interval
};

template <size_t streamingLowerBound_V, internal::simd::SIMDWidth simdWidth_V>
SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>::SIMDEncoderCommand(size_t symbolTablePrecision) : mStates{LowerBound}, mSymbolTablePrecision{symbolTablePrecision}, mNSamples{static_cast<double>(pow2(mSymbolTablePrecision))}
{
  if (mSymbolTablePrecision > LowerBound) {
    throw std::runtime_error(fmt::format("[{}]: SymbolTable Precision of {} Bits is larger than allowed by the rANS Encoder (max {} Bits)", __PRETTY_FUNCTION__, mSymbolTablePrecision, LowerBound));
  }
};

template <size_t streamingLowerBound_V, internal::simd::SIMDWidth simdWidth_V>
template <typename Stream_IT>
Stream_IT SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>::flush(Stream_IT iter)
{
  Stream_IT streamPos = iter;
  for (auto stateIter = mStates.rbegin(); stateIter != mStates.rend(); ++stateIter) {
    streamPos = flushState(*stateIter, streamPos);
  }
  return streamPos;
};

template <size_t streamingLowerBound_V, internal::simd::SIMDWidth simdWidth_V>
template <typename Stream_IT>
Stream_IT SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>::putSymbols(Stream_IT outputIter, const symbol_type& symbols)
{

  // can't encode symbol with freq=0
#if !defined(NDEBUG)
  // for (const auto& symbol : symbols) {
  //   assert(symbol->getFrequency() != 0);
  // }
#endif

  auto [streamPosition, renormedStates] = simd::ransRenorm<Stream_IT, LowerBound, StreamBits>(toConstSIMDView(mStates),
                                                                                              toConstSIMDView(symbols.frequencies),
                                                                                              static_cast<uint8_t>(mSymbolTablePrecision),
                                                                                              outputIter);
  ransEncode(toConstSIMDView(renormedStates).template subView<0, 1>(),
             simd::int32ToDouble<simdWidth_V>(toConstSIMDView(symbols.frequencies).template subView<0, 1>()),
             simd::int32ToDouble<simdWidth_V>(toConstSIMDView(symbols.cumulativeFrequencies).template subView<0, 1>()),
             toConstSIMDView(mNSamples),
             toSIMDView(mStates).template subView<0, 1>());
  ransEncode(toConstSIMDView(renormedStates).template subView<1, 1>(),
             simd::int32ToDouble<simdWidth_V>(toConstSIMDView(symbols.frequencies).template subView<1, 1>()),
             simd::int32ToDouble<simdWidth_V>(toConstSIMDView(symbols.cumulativeFrequencies).template subView<1, 1>()),
             toConstSIMDView(mNSamples),
             toSIMDView(mStates).template subView<1, 1>());

  return streamPosition;
}

template <size_t streamingLowerBound_V, internal::simd::SIMDWidth simdWidth_V>
template <typename Stream_IT>
Stream_IT SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>::putSymbols(Stream_IT outputIter, const symbol_type& symbols, size_t nActiveStreams)
{
  Stream_IT streamPos = outputIter;

  for (size_t i = nActiveStreams; i-- > 0;) {
    simd::Symbol encodeSymbol{symbols.frequencies[i], symbols.cumulativeFrequencies[i]};
    streamPos = putSymbol(streamPos, encodeSymbol, mStates[i]);
  }

  return streamPos;
};

template <size_t streamingLowerBound_V, internal::simd::SIMDWidth simdWidth_V>
template <typename Stream_IT>
Stream_IT SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>::putSymbol(Stream_IT outputIter, const simd::Symbol& symbol, state_type& state)
{
  assert(symbol.getFrequency() != 0); // can't encode symbol with freq=0
  // renormalize
  const auto [x, streamPos] = renorm(state, outputIter, symbol.getFrequency());

  // x = C(s,x)
  state = ((x / symbol.getFrequency()) << mSymbolTablePrecision) + (x % symbol.getFrequency()) + symbol.getCumulative();
  return streamPos;
}

template <size_t streamingLowerBound_V, internal::simd::SIMDWidth simdWidth_V>
template <typename Stream_IT>
Stream_IT SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>::flushState(state_type& state, Stream_IT iter)
{
  Stream_IT streamPosition = iter;

  ++streamPosition;
  *streamPosition = static_cast<stream_type>(state >> 32);
  ++streamPosition;
  *streamPosition = static_cast<stream_type>(state >> 0);

  state = 0;
  return streamPosition;
}

template <size_t streamingLowerBound_V, internal::simd::SIMDWidth simdWidth_V>
template <typename Stream_IT>
inline auto SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>::renorm(state_type state, Stream_IT outputIter, uint32_t frequency) -> std::tuple<state_type, Stream_IT>
{
  state_type maxState = ((LowerBound >> mSymbolTablePrecision) << StreamBits) * frequency; // this turns into a shift.
  if (state >= maxState) {
    ++outputIter;
    *outputIter = static_cast<stream_type>(state);
    state >>= StreamBits;
    assert(state < maxState);
  }
  return std::make_tuple(state, outputIter);
};

template <size_t streamingLowerBound_V>
using SSEEncoderCommand = SIMDEncoderCommand<streamingLowerBound_V, simd::SIMDWidth::SSE>;
template <size_t streamingLowerBound_V>
using AVXEncoderCommand = SIMDEncoderCommand<streamingLowerBound_V, simd::SIMDWidth::AVX>;

} // namespace internal
} // namespace rans
} // namespace o2

#endif /* RANS_INTERNAL_SIMDENCODERS_H */