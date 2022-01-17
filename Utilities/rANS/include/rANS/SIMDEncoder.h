// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   SIMDEncoder.h
/// @author Michael Lettrich
/// @since  2021-03-18
/// @brief  Encoder - code symbol into a rANS encoded state

#ifndef RANS_SIMDENCODER_H
#define RANS_SIMDENCODER_H

#include <memory>
#include <algorithm>
#include <iomanip>

#include <fairlogger/Logger.h>
#include <stdexcept>

#include "rANS/internal/backend/simd/Encoder.h"
#include "rANS/internal/backend/simd/Symbol.h"
#include "rANS/internal/helper.h"
#include "rANS/RenormedFrequencyTable.h"
#include "rANS/internal/backend/simd/SymbolTable.h"

namespace o2
{
namespace rans
{

template <typename coder_T,
          typename stream_T,
          typename source_T,
          uint8_t nStreams_V = 8>
class SIMDEncoder
{
 protected:
  using encoderSymbolTable_t = internal::simd::SymbolTable;

 public:
  //TODO(milettri): fix once ROOT cling respects the standard http://wg21.link/p1286r2
  SIMDEncoder() noexcept {}; //NOLINT
  SIMDEncoder(const RenormedFrequencyTable& frequencyTable) : mSymbolTable{frequencyTable} {};

  inline size_t getSymbolTablePrecision() const noexcept { return mSymbolTable.getPrecision(); };
  inline size_t getAlphabetRangeBits() const noexcept { return mSymbolTable.getAlphabetRangeBits(); };
  inline symbol_t getMinSymbol() const noexcept { return mSymbolTable.getMinSymbol(); };
  inline symbol_t getMaxSymbol() const noexcept { return mSymbolTable.getMaxSymbol(); };

  template <typename stream_IT, typename source_IT, std::enable_if_t<internal::isCompatibleIter_v<source_T, source_IT>, bool> = true>
  stream_IT process(source_IT inputBegin, source_IT inputEnd, stream_IT outputBegin) const;

 protected:
  encoderSymbolTable_t mSymbolTable{};
  size_t mSymbolTablePrecision{};

  //TODO(milettri): make this depend on hardware
  static constexpr size_t nHardwareStreams_V = 4;
  static constexpr size_t nInterleavedStreams_V = nStreams_V / nHardwareStreams_V;
  using ransCoder_t = typename internal::simd::Encoder<nHardwareStreams_V>;
};

template <typename coder_T, typename stream_T, typename source_T, uint8_t nStreams_V>
template <typename stream_IT,
          typename source_IT,
          std::enable_if_t<internal::isCompatibleIter_v<source_T, source_IT>, bool>>
stream_IT SIMDEncoder<coder_T, stream_T, source_T, nStreams_V>::process(source_IT inputBegin, source_IT inputEnd, stream_IT outputBegin) const
{
  using namespace internal;
  LOG(trace) << "start encoding";
  RANSTimer t;
  t.start();

  if (inputBegin == inputEnd) {
    LOG(warning) << "passed empty message to encoder, skip encoding";
    return outputBegin;
  }

  stream_IT outputIter = outputBegin;
  source_IT inputIT = inputEnd;

  auto maskedEncode = [this](source_IT symbolIter, stream_IT outputIter, ransCoder_t& coder, size_t nActiveStreams = nHardwareStreams_V) {
    std::array<internal::simd::Symbol, nHardwareStreams_V> encoderSymbols{};
    for (auto encSymbolIter = encoderSymbols.rend() - nActiveStreams; encSymbolIter != encoderSymbols.rend(); ++encSymbolIter) {
      const source_T symbol = *(--symbolIter);
      *encSymbolIter = (this->mSymbolTable)[symbol];
    }
    return std::tuple(symbolIter, coder.putSymbols(outputIter, encoderSymbols, nActiveStreams));
  };

  auto encode = [this](source_IT symbolIter, stream_IT outputIter, ransCoder_t& coder) {
    std::array<internal::simd::Symbol, nHardwareStreams_V> encoderSymbols{};
    for (auto encSymbolIter = encoderSymbols.rbegin(); encSymbolIter != encoderSymbols.rend(); ++encSymbolIter) {
      const source_T symbol = *(--symbolIter);
      *encSymbolIter = (this->mSymbolTable)[symbol];
    }
    return std::tuple(symbolIter, coder.putSymbols(outputIter, encoderSymbols));
  };

  // create coders
  std::array<ransCoder_t, nInterleavedStreams_V> interleavedCoders;
  for (auto& coder : interleavedCoders) {
    coder = ransCoder_t(this->getSymbolTablePrecision());
  }

  // calculate sizes and numbers of iterations:
  const auto inputBufferSize = std::distance(inputBegin, inputEnd);
  const size_t nMainLoopIterations = inputBufferSize / nStreams_V;
  const size_t nMainLoopRemainderSymbols = inputBufferSize % nStreams_V;
  const size_t nRemainderLoopIterations = nMainLoopRemainderSymbols / nHardwareStreams_V;
  const size_t nMaskedEncodes = nMainLoopRemainderSymbols % nHardwareStreams_V;

  LOG(trace) << "InputBufferSize: " << inputBufferSize;
  LOG(trace) << "Loops Main: " << nMainLoopIterations << ", Loops Remainder: " << nMainLoopRemainderSymbols << ", Masked Encodes :" << nMaskedEncodes;

  // iterator pointing to the active coder
  auto coderIter = std::rend(interleavedCoders) - nRemainderLoopIterations;

  if (nMaskedEncodes) {
    LOG(trace) << "masked encodes";
    // one more encoding step than nRemainderLoopIterations for masked encoding
    // will not cause out of range
    --coderIter;
    std::tie(inputIT, outputIter) = maskedEncode(inputIT, outputIter, *(coderIter), nMaskedEncodes);
    coderIter++;
  }
  // now encode the rest of the remainder symbols
  LOG(trace) << "remainder";
  for (coderIter; coderIter != std::rend(interleavedCoders); ++coderIter) {
    std::tie(inputIT, outputIter) = encode(inputIT, outputIter, *coderIter);
  }

  // main encoder loop
  LOG(trace) << "main loop";
  for (size_t i = 0; i < nMainLoopIterations; ++i) {
    for (coderIter = std::rbegin(interleavedCoders); coderIter != std::rend(interleavedCoders); ++coderIter) {
      std::tie(inputIT, outputIter) = encode(inputIT, outputIter, *coderIter);
    }
  }

  LOG(trace) << "flushing";
  for (coderIter = std::rbegin(interleavedCoders); coderIter != std::rend(interleavedCoders); ++coderIter) {
    outputIter = coderIter->flush(outputIter);
  }

  // first iterator past the range so that sizes, distances and iterators work correctly.
  ++outputIter;

  t.stop();
  LOG(debug1) << "Encoder::" << __func__ << " {ProcessedBytes: " << inputBufferSize * sizeof(source_T) << ","
              << " inclusiveTimeMS: " << t.getDurationMS() << ","
              << " BandwidthMiBPS: " << std::fixed << std::setprecision(2) << (inputBufferSize * sizeof(source_T) * 1.0) / (t.getDurationS() * 1.0 * (1 << 20)) << "}";

// advanced diagnostics for debug builds
#if !defined(NDEBUG)

  LOG(debug2) << "EncoderProperties: {"
              << "sourceTypeB: " << sizeof(source_T) << ", "
              << "streamTypeB: " << sizeof(stream_T) << ", "
              << "coderTypeB: " << sizeof(coder_T) << ", "
              << "symbolTablePrecision: " << this->mSymbolTablePrecision << ", "
              << "inputBufferSizeB: " << (inputBufferSize * sizeof(source_T)) << "}";
#endif

  LOG(trace) << "done encoding";

  return outputIter;
};

} // namespace rans
} // namespace o2

#endif /* RANS_SIMDENCODER_H */
