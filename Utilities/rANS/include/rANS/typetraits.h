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

#include "rANS/StaticFrequencyTable.h"
#include "rANS/DynamicFrequencyTable.h"
#include "rANS/HashFrequencyTable.h"

#include "rANS/RenormedFrequencies.h"

#include "rANS/StaticSymbolTable.h"
#include "rANS/DynamicSymbolTable.h"
#include "rANS/HashSymbolTable.h"

#include "rANS/internal/SingleStreamEncodeCommand.h"
#include "rANS/internal/SIMDEncodeCommand.h"

#include "rANS/EncoderFacade.h"
#include "rANS/internal/Symbol.h"

namespace o2
{
namespace rans
{

enum class ContainerTag : uint8_t { Dynamic,
                                    Static,
                                    Hash };

enum class CoderTag : uint8_t { Compat,
                                SingleStream,
                                SSE,
                                AVX2 };

inline constexpr size_t NStreams = 16;
inline constexpr size_t RenormingLowerBound = 20;

inline constexpr size_t LegacyNStreams = 2;
inline constexpr size_t LegacyRenormingLowerBound = 31;

template <typename T>
struct getContainerTag;

template <typename source_T>
struct getContainerTag<StaticFrequencyTable<source_T>> : public std::integral_constant<ContainerTag, ContainerTag::Dynamic> {
};
template <typename source_T>
struct getContainerTag<DynamicFrequencyTable<source_T>> : public std::integral_constant<ContainerTag, ContainerTag::Dynamic> {
};
template <typename source_T>
struct getContainerTag<HashFrequencyTable<source_T>> : public std::integral_constant<ContainerTag, ContainerTag::Hash> {
};
template <typename source_T>
struct getContainerTag<RenormedStaticFrequencyTable<source_T>> : public std::integral_constant<ContainerTag, ContainerTag::Static> {
};
template <typename source_T>
struct getContainerTag<RenormedDynamicFrequencyTable<source_T>> : public std::integral_constant<ContainerTag, ContainerTag::Dynamic> {
};
template <typename source_T>
struct getContainerTag<RenormedHashFrequencyTable<source_T>> : public std::integral_constant<ContainerTag, ContainerTag::Hash> {
};
template <typename source_T, typename symbol_T>
struct getContainerTag<StaticSymbolTable<source_T, symbol_T>> : public std::integral_constant<ContainerTag, ContainerTag::Static> {
};
template <typename source_T, typename symbol_T>
struct getContainerTag<DynamicSymbolTable<source_T, symbol_T>> : public std::integral_constant<ContainerTag, ContainerTag::Dynamic> {
};
template <typename source_T, typename symbol_T>
struct getContainerTag<HashSymbolTable<source_T, symbol_T>> : public std::integral_constant<ContainerTag, ContainerTag::Hash> {
};

template <class encoderCommand_T, class symbolTable_T, size_t nStreams_V>
struct getContainerTag<EncoderFacade<encoderCommand_T, symbolTable_T, nStreams_V>> : public getContainerTag<symbolTable_T> {
};

template <typename T>
inline constexpr ContainerTag getContainerTag_v = getContainerTag<T>::value;

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

template <ContainerTag tag_V>
struct ContainerTraits {
};

template <>
struct ContainerTraits<ContainerTag::Dynamic> {

  template <typename source_T>
  using frequencyTable_type = DynamicFrequencyTable<source_T>;
  template <typename source_T>
  using renormedFrequencyTable_type = RenormedDynamicFrequencyTable<source_T>;
  template <typename source_T, class symbol_T>
  using symbolTable_type = DynamicSymbolTable<source_T, symbol_T>;
};

template <>
struct ContainerTraits<ContainerTag::Static> {

  template <typename source_T>
  using frequencyTable_type = StaticFrequencyTable<source_T>;
  template <typename source_T>
  using renormedFrequencyTable_type = RenormedStaticFrequencyTable<source_T>;
  template <typename source_T, class symbol_T>
  using symbolTable_type = StaticSymbolTable<source_T, symbol_T>;
};

template <>
struct ContainerTraits<ContainerTag::Hash> {

  template <typename source_T>
  using frequencyTable_type = HashFrequencyTable<source_T>;
  template <typename source_T>
  using renormedFrequencyTable_type = RenormedHashFrequencyTable<source_T>;
  template <typename source_T, class symbol_T>
  using symbolTable_type = HashSymbolTable<source_T, symbol_T>;
};

template <typename T>
struct isSymbolTableContainer : public std::is_base_of<internal::SymbolTableContainer<typename T::source_type,
                                                                                      typename T::index_type,
                                                                                      typename T::value_type,
                                                                                      typename T::container_type, T>,
                                                       T> {
};

template <typename T>
inline constexpr bool isSymbolTableContainer_v = isSymbolTableContainer<T>::value;

template <typename T>
struct isFrequencyTable : std::false_type {
};

template <typename source_T>
struct isFrequencyTable<DynamicFrequencyTable<source_T>> : std::true_type {
};

template <typename source_T>
struct isFrequencyTable<StaticFrequencyTable<source_T>> : std::true_type {
};

template <typename source_T>
struct isFrequencyTable<HashFrequencyTable<source_T>> : std::true_type {
};

template <typename T>
inline constexpr bool isFrequencyTable_v = isFrequencyTable<T>::value;

template <typename T>
struct isFrequencyContainer : public std::is_base_of<internal::FrequencyContainer<typename T::source_type,
                                                                                  typename T::index_type,
                                                                                  typename T::value_type,
                                                                                  typename T::container_type, T>,
                                                     T> {
};

template <typename T>
inline constexpr bool isFrequencyContainer_v = isFrequencyContainer<T>::value;

template <typename T>
struct isRenormedFrequencyTable : std::false_type {
};

template <typename source_T>
struct isRenormedFrequencyTable<RenormedStaticFrequencyTable<source_T>> : std::true_type {
};

template <typename source_T>
struct isRenormedFrequencyTable<RenormedDynamicFrequencyTable<source_T>> : std::true_type {
};

template <typename source_T>
struct isRenormedFrequencyTable<RenormedHashFrequencyTable<source_T>> : std::true_type {
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

template <typename source_T>
inline constexpr ContainerTag getPreferedContainerTag()
{
  if constexpr (sizeof(source_T) < 4) {
    return ContainerTag::Static;
  } else {
    return ContainerTag::Dynamic;
  }
};

template <ContainerTag tag_V>
struct makeFrequencyTable {

  template <typename source_IT>
  [[nodiscard]] inline static decltype(auto) fromSamples(source_IT begin, source_IT end)
  {
    using source_type = typename std::iterator_traits<source_IT>::value_type;
    using frequencyTable_type = typename ContainerTraits<tag_V>::template frequencyTable_type<source_type>;

    frequencyTable_type f{};
    f.addSamples(begin, end);
    return f;
  };

  template <typename source_T>
  [[nodiscard]] inline static decltype(auto) fromSamples(gsl::span<const source_T> range)
  {
    using source_type = typename std::remove_cv_t<source_T>;
    using frequencyTable_type = typename ContainerTraits<tag_V>::template frequencyTable_type<source_type>;

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
  template <typename renormedFrequencyTable_T>
  [[nodiscard]] inline static constexpr decltype(auto) fromRenormed(const renormedFrequencyTable_T& renormed)
  {
    constexpr ContainerTag containerTag = getContainerTag_v<renormedFrequencyTable_T>;
    constexpr CoderTag coderTag = coderTag_V;
    using source_type = typename renormedFrequencyTable_T::source_type;
    using symbol_type = typename SymbolTraits<coderTag>::type;
    using coder_command = typename CoderTraits<coderTag>::template type<this_type::RenormingLowerBound>;
    using symbolTable_type = typename ContainerTraits<containerTag>::template symbolTable_type<source_type, symbol_type>;
    using encoderType = EncoderFacade<coder_command, symbolTable_type, this_type::NStreams>;

    return encoderType{renormed};
  };

  template <typename frequencyTable_T>
  [[nodiscard]] inline static decltype(auto) fromFrequencyTable(frequencyTable_T&& frequencyTable, size_t renormingPrecision = 0)
  {
    const auto renormedFrequencies = renormCutoffIncompressible(std::forward<frequencyTable_T>(frequencyTable), renormingPrecision);
    return this_type::fromRenormed(renormedFrequencies);
  };

  template <typename source_IT, ContainerTag containerTag_V = getPreferedContainerTag<typename std::iterator_traits<source_IT>::value_type>()>
  [[nodiscard]] inline static decltype(auto) fromSamples(source_IT begin, source_IT end, size_t renormingPrecision = 0)
  {
    auto frequencyTable = makeFrequencyTable<containerTag_V>::template fromSamples(begin, end);

    return this_type::fromFrequencyTable(std::move(frequencyTable), renormingPrecision);
  };

  template <typename source_T, ContainerTag containerTag_V = getPreferedContainerTag<source_T>()>
  [[nodiscard]] inline static decltype(auto) fromSamples(gsl::span<const source_T> range, size_t renormingPrecision = 0)
  {
    auto frequencyTable = makeFrequencyTable<containerTag_V>::template fromSamples(range);
    return this_type::fromFrequencyTable(std::move(frequencyTable), renormingPrecision);
  };

 private:
  static constexpr size_t NStreams = nStreams_V;
  static constexpr size_t RenormingLowerBound = renormingLowerBound_V;
};

} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_TYPETRAITS_H_ */