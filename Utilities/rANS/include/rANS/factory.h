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

/// @file   factories.h
/// @author michael.lettrich@cern.ch
/// @brief  static factory classes for building histograms, encoders and decoders.

#ifndef RANS_FACTORY_H_
#define RANS_FACTORY_H_

#include "rANS/internal/common/typetraits.h"

#include "rANS/internal/containers/Histogram.h"
#include "rANS/internal/containers/RenormedHistogram.h"
#include "rANS/internal/containers/SymbolTable.h"
#include "rANS/internal/containers/Symbol.h"

#include "rANS/internal/encode/Encoder.h"
#include "rANS/internal/encode/SingleStreamEncoderImpl.h"
#include "rANS/internal/encode/SIMDEncoderImpl.h"

#include "rANS/internal/decode/Decoder.h"
#include "rANS/internal/decode/DecoderImpl.h"

namespace o2::rans
{

struct makeHistogram {

  template <typename source_IT>
  [[nodiscard]] inline static decltype(auto) fromSamples(source_IT begin, source_IT end)
  {
    using source_type = typename std::iterator_traits<source_IT>::value_type;
    using histogram_type = Histogram<source_type>;

    histogram_type f{};
    f.addSamples(begin, end);
    return f;
  };

  template <typename source_T>
  [[nodiscard]] inline static decltype(auto) fromSamples(gsl::span<const source_T> range)
  {
    using source_type = typename std::remove_cv_t<source_T>;
    using histogram_type = Histogram<source_type>;

    histogram_type f;
    f.addSamples(range);
    return f;
  };
};

template <CoderTag coderTag_V, size_t nStreams_V = internal::NStreams, size_t renormingLowerBound_V = internal::RenormingLowerBound>
class makeEncoder
{

  using this_type = makeEncoder<coderTag_V, nStreams_V, renormingLowerBound_V>;

 public:
  template <typename source_T>
  [[nodiscard]] inline static constexpr decltype(auto) fromRenormed(const RenormedHistogram<source_T>& renormed)
  {
    using namespace internal;
    constexpr CoderTag coderTag = coderTag_V;
    using source_type = source_T;
    using symbol_type = typename SymbolTraits<coderTag>::type;
    using coder_command = typename CoderTraits<coderTag>::template type<this_type::RenormingLowerBound>;
    using symbolTable_type = SymbolTable<source_type, symbol_type>;
    using encoderType = Encoder<coder_command, symbolTable_type, this_type::NStreams>;

    return encoderType{renormed};
  };

  template <typename source_T>
  [[nodiscard]] inline static decltype(auto) fromHistogram(Histogram<source_T>&& histogram, size_t renormingPrecision = 0)
  {
    const auto renormedHistogram = renormCutoffIncompressible(std::forward<Histogram<source_T>>(histogram), renormingPrecision);
    return this_type::fromRenormed(renormedHistogram);
  };

  template <typename source_IT>
  [[nodiscard]] inline static decltype(auto) fromSamples(source_IT begin, source_IT end, size_t renormingPrecision = 0)
  {
    auto histogram = makeHistogram::fromSamples(begin, end);

    return this_type::fromHistogram(std::move(histogram), renormingPrecision);
  };

  template <typename source_T>
  [[nodiscard]] inline static decltype(auto) fromSamples(gsl::span<const source_T> range, size_t renormingPrecision = 0)
  {
    auto histogram = makeHistogram::template fromSamples(range);
    return this_type::fromHistogram(std::move(histogram), renormingPrecision);
  };

 private:
  static constexpr size_t NStreams = nStreams_V;
  static constexpr size_t RenormingLowerBound = renormingLowerBound_V;
};

template <size_t nStreams_V = internal::NStreams, size_t renormingLowerBound_V = internal::RenormingLowerBound>
class makeDecoder
{

  using this_type = makeDecoder;

 public:
  template <typename source_T>
  [[nodiscard]] inline static constexpr decltype(auto) fromRenormed(const RenormedHistogram<source_T>& renormed)
  {
    using namespace internal;

    using source_type = source_T;
    using coder_type = DecoderImpl<RenormingLowerBound>;
    using symbol_type = typename coder_type::symbol_type;
    using symbolTable_type = SymbolTable<source_type, symbol_type>;
    using decoder_type = Decoder<coder_type, symbolTable_type>;

    return decoder_type{renormed};
  };

  template <typename source_T>
  [[nodiscard]] inline static decltype(auto) fromHistogram(Histogram<source_T>&& histogram, size_t renormingPrecision = 0)
  {
    const auto renormedHistogram = renormCutoffIncompressible(std::forward<Histogram<source_T>>(histogram), renormingPrecision);
    return this_type::fromRenormed(renormedHistogram);
  };

  template <typename source_IT>
  [[nodiscard]] inline static decltype(auto) fromSamples(source_IT begin, source_IT end, size_t renormingPrecision = 0)
  {
    auto histogram = makeHistogram::fromSamples(begin, end);
    return this_type::fromHistogram(std::move(histogram), renormingPrecision);
  };

  template <typename source_T>
  [[nodiscard]] inline static decltype(auto) fromSamples(gsl::span<const source_T> range, size_t renormingPrecision = 0)
  {
    auto histogram = makeHistogram::fromSamples(range);
    return this_type::fromHistogram(std::move(histogram), renormingPrecision);
  };
};

} // namespace o2::rans

#endif /* RANS_FACTORY_H_ */