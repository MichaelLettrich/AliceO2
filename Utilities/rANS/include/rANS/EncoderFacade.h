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

/// @file   EncoderBase.h
/// @author michael.lettrich@cern.ch
/// @since  Feb 8, 2021
/// @brief

#ifndef INCLUDE_RANS_ENCODERFACADE_H_
#define INCLUDE_RANS_ENCODERFACADE_H_

#include <algorithm>
#include <iomanip>
#include <memory>

#include <fairlogger/Logger.h>
#include <stdexcept>

#include <gsl/span>

#include "rANS/definitions.h"
#include "rANS/internal/SymbolMapper.h"
#include "rANS/internal/helper.h"

namespace o2
{
namespace rans
{

template <class encoder_T, class symbolTable_T, std::size_t nStreams_V>
class EncoderFacade
{
 public:
  using symbolTable_type = symbolTable_T;
  using symbol_type = typename symbolTable_T::value_type;
  using coder_type = encoder_T;
  using source_type = typename symbolTable_type::source_type;
  using stream_type = typename coder_type::stream_type;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  static constexpr size_type NStreams = nStreams_V;

  EncoderFacade() = default;
  template <typename renormedSymbolTable_T>
  EncoderFacade(const renormedSymbolTable_T& renormedFrequencyTable) : mSymbolTable{renormedFrequencyTable} {};

  [[nodiscard]] inline const symbolTable_type& getSymbolTable() const noexcept { return mSymbolTable; };

  [[nodiscard]] inline static constexpr size_type getNStreams() noexcept { return NStreams; };

  template <typename stream_IT, typename source_IT, typename literals_IT = std::nullptr_t, std::enable_if_t<internal::isCompatibleIter_v<source_type, source_IT>, bool> = true>
  stream_IT process(source_IT inputBegin, source_IT inputEnd, stream_IT outputBegin, literals_IT literalsBegin = nullptr) const;

  template <typename stream_IT, typename source_IT, typename literals_IT = std::nullptr_t, std::enable_if_t<internal::isCompatibleIter_v<source_type, source_IT>, bool> = true>
  inline stream_IT process(gsl::span<const source_type> inputStream, gsl::span<stream_type> outputStream, literals_IT literalsBegin = nullptr) const
  {
    return process(inputStream.begin(), inputStream.begin() + inputStream.size(), outputStream.begin(), literalsBegin);
  };

 protected:
  symbolTable_type mSymbolTable{};

  static constexpr size_type NCoderStreams = coder_type::getNstreams();
  static constexpr size_type NCoders = NStreams / NCoderStreams;

  //compile time preconditions:
  static_assert(internal::isPow2(nStreams_V));
  static_assert(coder_type::getNstreams() <= EncoderFacade::getNStreams());
  static_assert(NCoders * NCoderStreams == NStreams);
};

template <class encoder_T, class symbolTable_T, std::size_t nStreams_V>
template <typename stream_IT, typename source_IT, typename literals_IT, std::enable_if_t<internal::isCompatibleIter_v<typename symbolTable_T::source_type, source_IT>, bool>>
stream_IT EncoderFacade<encoder_T, symbolTable_T, nStreams_V>::process(source_IT inputBegin, source_IT inputEnd, stream_IT outputBegin, literals_IT literalsBegin) const
{

  using namespace internal;

  if (inputBegin == inputEnd) {
    LOG(warning) << "passed empty message to encoder, skip encoding";
    return outputBegin;
  }

  std::array<coder_type, NCoders> coders;

  for (auto& coder : coders) {
    coder = coder_type{mSymbolTable.getPrecision()};
  }

  // calculate sizes and numbers of iterations:
  const auto inputBufferSize = std::distance(inputBegin, inputEnd); // size of the input buffer
  const size_t nAllCoderIterations = inputBufferSize / NStreams;    // number
  const size_t nRemainderSymbols = inputBufferSize % NStreams;
  const size_t nPartialCoderIterations = nRemainderSymbols / NCoderStreams;
  const size_t nFractionalEncodes = nRemainderSymbols % NCoderStreams;

  LOGP(info, "inputBufferSize: {}", inputBufferSize);
  LOGP(info, "nAllCoderIterations: {}", nAllCoderIterations);
  LOGP(info, "nRemainderSymbols: {}", nRemainderSymbols);
  LOGP(info, "nPartialCoderIterations: {}", nPartialCoderIterations);
  LOGP(info, "nFractionalEncodes: {}", nFractionalEncodes);

  // from here on, everything runs backwards!
  // We are encoding symbols from the end of the message to the beginning of the message.
  // For consistency, also coders have to run backwards, i.e. coder n+1 runs before coder n.
  // This allows decoding to happen the "natural way"
  // To keep track of iterators we use reverse iterators which makes them appear to run in a forward order again.

  stream_IT outputIter = outputBegin;
  source_IT inputIter = inputEnd;
  --inputIter;
  source_IT inputREnd = inputBegin;
  --inputREnd;
  literals_IT literalsIter = literalsBegin;

  SymbolMapper<symbolTable_type, coder_type, literals_IT> symbolMapper{this->mSymbolTable, literalsIter};

  auto activeCoder = coders.rend() - nPartialCoderIterations;

  uint32_t counter = 0;

  if (nFractionalEncodes) {
    // LOG(trace) << "masked encodes";
    // one more encoding step than nRemainderLoopIterations for masked encoding
    // will not cause out of range
    --activeCoder;
    typename coder_type::symbol_type encoderSymbol;
    inputIter = symbolMapper.unpackSymbols(inputIter, encoderSymbol, nFractionalEncodes);
    outputIter = activeCoder->putSymbols(outputIter, encoderSymbol, nFractionalEncodes);
    ++counter;

    ++activeCoder;
  }

  // we are encoding backwards!
  while (inputIter != inputREnd) {
    // iterate over coders with wrap around
    for (; activeCoder != coders.rend(); ++activeCoder) {
      typename coder_type::symbol_type encoderSymbol;
      inputIter = symbolMapper.unpackSymbols(inputIter, encoderSymbol);
      outputIter = activeCoder->putSymbols(outputIter, encoderSymbol);
      ++counter;
    }
    activeCoder = coders.rbegin();
  }

  // LOG(trace) << "flushing";
  for (activeCoder = std::rbegin(coders); activeCoder != std::rend(coders); ++activeCoder) {
    outputIter = activeCoder->flush(outputIter);
  }

  LOGP(info, "nEncodes: {}", counter);

  // first iterator past the range so that sizes, distances and iterators work correctly.
  ++outputIter;

  return outputIter;
}

}; // namespace rans
}; // namespace o2

#endif /* INCLUDE_RANS_ENCODERFACADE_H_ */
