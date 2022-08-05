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

#ifndef INCLUDE_RANS_TYPETRAITS_H_
#define INCLUDE_RANS_TYPETRAITS_H_

#include "rANS/FrequencyTable.h"
#include "rANS/RenormedFrequencies.h"
#include "rANS/SymbolTable.h"

#include "rANS/internal/SingleStreamEncodeCommand.h"
#include "rANS/internal/SIMDEncodeCommand.h"

#include "rANS/EncoderFacade.h"
#include "rANS/internal/Symbol.h"
#include "rANS/DecoderFacade.h"
#include "rANS/internal/Decoder.h"

namespace o2
{
namespace rans
{

enum class CoderTag : uint8_t { Compat,
                                SingleStream,
                                SSE,
                                AVX2 };

inline constexpr size_t NStreams = 16;
inline constexpr size_t RenormingLowerBound = 20;

inline constexpr size_t LegacyNStreams = 2;
inline constexpr size_t LegacyRenormingLowerBound = 31;

template <typename T>
struct getCoderTag;

template <size_t V>
struct getCoderTag<internal::CompatEncoderCommand<V>> : public std::integral_constant<CoderTag, CoderTag::Compat> {
};

template <size_t V>
struct getCoderTag<internal::SingleStreamEncoderCommand<V>> : public std::integral_constant<CoderTag, CoderTag::SingleStream> {
};

template <size_t V>
struct getCoderTag<internal::SSEEncoderCommand<V>> : public std::integral_constant<CoderTag, CoderTag::SSE> {
};

template <size_t V>
struct getCoderTag<internal::AVXEncoderCommand<V>> : public std::integral_constant<CoderTag, CoderTag::AVX2> {
};

template <class encoderCommand_T, class symbolTable_T, size_t nStreams_V>
struct getCoderTag<EncoderFacade<encoderCommand_T, symbolTable_T, nStreams_V>> : public getCoderTag<encoderCommand_T> {
};

template <typename T>
inline constexpr CoderTag getCoderTag_v = getCoderTag<T>::value;

template <typename T>
struct isSymbolTable : std::false_type {
};

template <typename source_T, typename value_T>
struct isSymbolTable<SymbolTable<source_T, value_T>> : std::true_type {
};

template <typename T>
inline constexpr bool isSymbolTable_v = isSymbolTable<T>::value;

template <typename T>
struct isFrequencyTable : std::false_type {
};

template <typename source_T>
struct isFrequencyTable<FrequencyTable<source_T>> : std::true_type {
};

template <typename T>
inline constexpr bool isFrequencyTable_v = isFrequencyTable<T>::value;

template <typename T>
struct isFrequencyContainer : public std::is_base_of<internal::FrequencyContainer<typename T::source_type, T>, T> {
};

template <typename T>
inline constexpr bool isFrequencyContainer_v = isFrequencyContainer<T>::value;

template <typename T>
struct isRenormedFrequencyTable : std::false_type {
};

template <typename source_T>
struct isRenormedFrequencyTable<RenormedFrequencyTable<source_T>> : std::true_type {
};

template <typename T>
inline constexpr bool isRenormedFrequencyTable_v = isRenormedFrequencyTable<T>::value;

template <CoderTag tag_V>
struct SymbolTraits {
  using type = internal::Symbol;
};

template <>
struct SymbolTraits<CoderTag::SingleStream> {
  using type = internal::PrecomputedSymbol;
};

template <CoderTag tag_V>
struct CoderTraits {
};

template <>
struct CoderTraits<CoderTag::Compat> {

  template <size_t lowerBound_V>
  using type = internal::CompatEncoderCommand<lowerBound_V>;
};

template <>
struct CoderTraits<CoderTag::SingleStream> {

  template <size_t lowerBound_V>
  using type = internal::SingleStreamEncoderCommand<lowerBound_V>;
};

template <>
struct CoderTraits<CoderTag::SSE> {

  template <size_t lowerBound_V>
  using type = internal::SSEEncoderCommand<lowerBound_V>;
};

template <>
struct CoderTraits<CoderTag::AVX2> {

  template <size_t lowerBound_V>
  using type = internal::AVXEncoderCommand<lowerBound_V>;
};

struct makeFrequencyTable {

  template <typename source_IT>
  [[nodiscard]] inline static decltype(auto) fromSamples(source_IT begin, source_IT end)
  {
    using source_type = typename std::iterator_traits<source_IT>::value_type;
    using frequencyTable_type = FrequencyTable<source_type>;

    frequencyTable_type f{};
    f.addSamples(begin, end);
    return f;
  };

  template <typename source_T>
  [[nodiscard]] inline static decltype(auto) fromSamples(gsl::span<const source_T> range)
  {
    using source_type = typename std::remove_cv_t<source_T>;
    using frequencyTable_type = FrequencyTable<source_type>;

    frequencyTable_type f;
    f.addSamples(range);
    return f;
  };
};

template <CoderTag coderTag_V, size_t nStreams_V = NStreams, size_t renormingLowerBound_V = RenormingLowerBound>
class makeEncoder
{

  using this_type = makeEncoder<coderTag_V, nStreams_V, renormingLowerBound_V>;

 public:
  template <typename source_T>
  [[nodiscard]] inline static constexpr decltype(auto) fromRenormed(const RenormedFrequencyTable<source_T>& renormed)
  {
    constexpr CoderTag coderTag = coderTag_V;
    using source_type = source_T;
    using symbol_type = typename SymbolTraits<coderTag>::type;
    using coder_command = typename CoderTraits<coderTag>::template type<this_type::RenormingLowerBound>;
    using symbolTable_type = SymbolTable<source_type, symbol_type>;
    using encoderType = EncoderFacade<coder_command, symbolTable_type, this_type::NStreams>;

    return encoderType{renormed};
  };

  template <typename source_T>
  [[nodiscard]] inline static decltype(auto) fromFrequencyTable(FrequencyTable<source_T>&& frequencyTable, size_t renormingPrecision = 0)
  {
    const auto renormedFrequencies = renormCutoffIncompressible(std::forward<FrequencyTable<source_T>>(frequencyTable), renormingPrecision);
    return this_type::fromRenormed(renormedFrequencies);
  };

  template <typename source_IT>
  [[nodiscard]] inline static decltype(auto) fromSamples(source_IT begin, source_IT end, size_t renormingPrecision = 0)
  {
    auto frequencyTable = makeFrequencyTable::fromSamples(begin, end);

    return this_type::fromFrequencyTable(std::move(frequencyTable), renormingPrecision);
  };

  template <typename source_T>
  [[nodiscard]] inline static decltype(auto) fromSamples(gsl::span<const source_T> range, size_t renormingPrecision = 0)
  {
    auto frequencyTable = makeFrequencyTable::template fromSamples(range);
    return this_type::fromFrequencyTable(std::move(frequencyTable), renormingPrecision);
  };

 private:
  static constexpr size_t NStreams = nStreams_V;
  static constexpr size_t RenormingLowerBound = renormingLowerBound_V;
};

template <size_t nStreams_V = NStreams, size_t renormingLowerBound_V = RenormingLowerBound>
class makeDecoder
{

  using this_type = makeDecoder;

 public:
  template <typename source_T>
  [[nodiscard]] inline static constexpr decltype(auto) fromRenormed(const RenormedFrequencyTable<source_T>& renormed)
  {
    using source_type = source_T;
    using coder_type = internal::Decoder<RenormingLowerBound>;
    using symbol_type = typename coder_type::symbol_type;
    using symbolTable_type = SymbolTable<source_type, symbol_type>;
    using decoder_type = DecoderFacade<coder_type, symbolTable_type>;

    return decoder_type{renormed};
  };

  template <typename source_T>
  [[nodiscard]] inline static decltype(auto) fromFrequencyTable(FrequencyTable<source_T>&& frequencyTable, size_t renormingPrecision = 0)
  {
    const auto renormedFrequencies = renormCutoffIncompressible(std::forward<FrequencyTable<source_T>>(frequencyTable), renormingPrecision);
    return this_type::fromRenormed(renormedFrequencies);
  };

  template <typename source_IT>
  [[nodiscard]] inline static decltype(auto) fromSamples(source_IT begin, source_IT end, size_t renormingPrecision = 0)
  {
    auto frequencyTable = makeFrequencyTable::fromSamples(begin, end);
    return this_type::fromFrequencyTable(std::move(frequencyTable), renormingPrecision);
  };

  template <typename source_T>
  [[nodiscard]] inline static decltype(auto) fromSamples(gsl::span<const source_T> range, size_t renormingPrecision = 0)
  {
    auto frequencyTable = makeFrequencyTable::fromSamples(range);
    return this_type::fromFrequencyTable(std::move(frequencyTable), renormingPrecision);
  };
};

} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_TYPETRAITS_H_ */