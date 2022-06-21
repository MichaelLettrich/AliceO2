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

#include "rANS/internal/EncodeCommandInterface.h"
#include "rANS/internal/helper.h"
#include "rANS/internal/Symbol.h"

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
template <size_t streamingLowerBound_V, simd::SIMDWidth simdWidth_V>
class SIMDEncoderCommand : public EncodeCommandInterface<simd::UnrolledSymbols,
                                                         SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>>
{
  using base_type = EncodeCommandInterface<simd::UnrolledSymbols, SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>>;

 public:
  using stream_type = typename base_type::stream_type;
  using state_type = typename base_type::state_type;
  using symbol_type = typename base_type::symbol_type;
  using size_type = typename base_type::size_type;
  using difference_type = typename base_type::difference_type;

  static_assert(streamingLowerBound_V <= 20, "SIMD coders are limited to 20 BIT precision because of their used of FP arithmeric");

  [[nodiscard]] inline static constexpr size_type getNstreams() noexcept { return 2 * simd::getElementCount<state_type>(simdWidth_V); };

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
  size_t mSymbolTablePrecision{};
  simd::simdI_t<simdWidth_V> mStates[2]{};
  simd::simdD_t<simdWidth_V> mNSamples{};

  template <typename Stream_IT>
  Stream_IT putSymbol(Stream_IT outputIter, const Symbol& symbol, state_type& state);

  template <typename Stream_IT>
  Stream_IT flushState(state_type& state, Stream_IT outputIter);

  // Renormalize the encoder.
  template <typename Stream_IT>
  std::tuple<state_type, Stream_IT> renorm(state_type state, Stream_IT outputIter, uint32_t frequency);

  inline static constexpr state_type LowerBound = pow2(streamingLowerBound_V); // lower bound of our normalization interval

  inline static constexpr state_type StreamBits = toBits(sizeof(stream_type)); // lower bound of our normalization interval
};

template <size_t streamingLowerBound_V, simd::SIMDWidth simdWidth_V>
SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>::SIMDEncoderCommand(size_t symbolTablePrecision) : mSymbolTablePrecision{symbolTablePrecision}, mStates{}, mNSamples{}
{
  if (mSymbolTablePrecision > LowerBound) {
    throw std::runtime_error(fmt::format("[{}]: SymbolTable Precision of {} Bits is larger than allowed by the rANS Encoder (max {} Bits)", __PRETTY_FUNCTION__, mSymbolTablePrecision, LowerBound));
  }

  mStates[0] = simd::setAll<simdWidth_V>(LowerBound);
  mStates[1] = simd::setAll<simdWidth_V>(LowerBound);

  mNSamples = simd::setAll<simdWidth_V>(static_cast<double>(pow2(mSymbolTablePrecision)));
};

template <size_t streamingLowerBound_V, simd::SIMDWidth simdWidth_V>
template <typename Stream_IT>
Stream_IT SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>::flush(Stream_IT iter)
{
  using namespace simd;
  epi64_t<simdWidth_V, 2> states;
  store(mStates[0], toSIMDView(states).template subView<0, 1>());
  store(mStates[1], toSIMDView(states).template subView<1, 1>());

  Stream_IT streamPos = iter;
  for (auto stateIter = states.rbegin(); stateIter != states.rend(); ++stateIter) {
    streamPos = flushState(*stateIter, streamPos);
  }

  mStates[0] = load(toConstSIMDView(states).template subView<0, 1>());
  mStates[1] = load(toConstSIMDView(states).template subView<1, 1>());

  return streamPos;
};

template <size_t streamingLowerBound_V, simd::SIMDWidth simdWidth_V>
template <typename Stream_IT>
inline Stream_IT SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>::putSymbols(Stream_IT outputIter, const symbol_type& symbols)
{
  using namespace simd;

#if !defined(NDEBUG)
  // for (const auto& symbol : symbols) {
//   //   assert(symbol->getFrequency() != 0);
// }
#endif
  simd::simdI_t<simdWidth_V> renormedStates[2];
  auto streamPosition = ransRenorm<Stream_IT, LowerBound, StreamBits>(mStates,
                                                                      symbols.frequencies,
                                                                      static_cast<uint8_t>(mSymbolTablePrecision),
                                                                      outputIter,
                                                                      renormedStates);
  mStates[0] = ransEncode(renormedStates[0], int32ToDouble<simdWidth_V>(symbols.frequencies[0]), int32ToDouble<simdWidth_V>(symbols.cumulativeFrequencies[0]), mNSamples);
  mStates[1] = ransEncode(renormedStates[1], int32ToDouble<simdWidth_V>(symbols.frequencies[1]), int32ToDouble<simdWidth_V>(symbols.cumulativeFrequencies[1]), mNSamples);

  return streamPosition;
}

template <size_t streamingLowerBound_V, simd::SIMDWidth simdWidth_V>
template <typename Stream_IT>
Stream_IT SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>::putSymbols(Stream_IT outputIter, const symbol_type& symbols, size_t nActiveStreams)
{
  using namespace simd;

  Stream_IT streamPos = outputIter;

  epi64_t<simdWidth_V, 2> states;
  store(mStates[0], toSIMDView(states).template subView<0, 1>());
  store(mStates[1], toSIMDView(states).template subView<1, 1>());

  epi32_t<SIMDWidth::SSE, 2> frequencies;
  epi32_t<SIMDWidth::SSE, 2> cumulativeFrequencies;

  store<uint32_t>(symbols.frequencies[0], toSIMDView(frequencies).template subView<0, 1>());
  store<uint32_t>(symbols.frequencies[1], toSIMDView(frequencies).template subView<1, 1>());
  store<uint32_t>(symbols.cumulativeFrequencies[0], toSIMDView(cumulativeFrequencies).template subView<0, 1>());
  store<uint32_t>(symbols.cumulativeFrequencies[1], toSIMDView(cumulativeFrequencies).template subView<1, 1>());

  for (size_t i = nActiveStreams; i-- > 0;) {
    Symbol encodeSymbol{frequencies[i], cumulativeFrequencies[i]};
    streamPos = putSymbol(streamPos, encodeSymbol, states[i]);
  }

  mStates[0] = load(toConstSIMDView(states).template subView<0, 1>());
  mStates[1] = load(toConstSIMDView(states).template subView<1, 1>());

  return streamPos;
};

template <size_t streamingLowerBound_V, simd::SIMDWidth simdWidth_V>
template <typename Stream_IT>
Stream_IT SIMDEncoderCommand<streamingLowerBound_V, simdWidth_V>::putSymbol(Stream_IT outputIter, const Symbol& symbol, state_type& state)
{
  assert(symbol.getFrequency() != 0); // can't encode symbol with freq=0
  // renormalize
  const auto [x, streamPos] = renorm(state, outputIter, symbol.getFrequency());

  // x = C(s,x)
  state = ((x / symbol.getFrequency()) << mSymbolTablePrecision) + (x % symbol.getFrequency()) + symbol.getCumulative();
  return streamPos;
}

template <size_t streamingLowerBound_V, simd::SIMDWidth simdWidth_V>
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

template <size_t streamingLowerBound_V, simd::SIMDWidth simdWidth_V>
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