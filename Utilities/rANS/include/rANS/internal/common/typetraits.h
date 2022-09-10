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

/// @file   typetraits.h
/// @author michael.lettrich@cern.ch
/// @brief  manipulation of types at compile time

#ifndef RANS_INTERNAL_COMMON_TYPETRAITS_H_
#define RANS_INTERNAL_COMMON_TYPETRAITS_H_

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

enum class CoderTag : uint8_t { Compat,
                                SingleStream,
                                SSE,
                                AVX2 };

namespace internal
{

inline constexpr size_t NStreams = 16;
inline constexpr size_t RenormingLowerBound = 20;

inline constexpr size_t LegacyNStreams = 2;
inline constexpr size_t LegacyRenormingLowerBound = 31;

template <typename T>
struct getCoderTag;

template <size_t V>
struct getCoderTag<CompatEncoderImpl<V>> : public std::integral_constant<CoderTag, CoderTag::Compat> {
};

template <size_t V>
struct getCoderTag<SingleStreamEncoderImpl<V>> : public std::integral_constant<CoderTag, CoderTag::SingleStream> {
};

template <size_t V>
struct getCoderTag<SSEEncoderImpl<V>> : public std::integral_constant<CoderTag, CoderTag::SSE> {
};

template <size_t V>
struct getCoderTag<AVXEncoderImpl<V>> : public std::integral_constant<CoderTag, CoderTag::AVX2> {
};

template <class encoderImpl_T, class symbolTable_T, size_t nStreams_V>
struct getCoderTag<Encoder<encoderImpl_T, symbolTable_T, nStreams_V>> : public getCoderTag<encoderImpl_T> {
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
struct isHistogram : std::false_type {
};

template <typename source_T>
struct isHistogram<Histogram<source_T>> : std::true_type {
};

template <typename T>
inline constexpr bool isHistogram_v = isHistogram<T>::value;

template <typename T>
struct isCountingContainer : public std::is_base_of<CountingContainer<typename T::source_type, T>, T> {
};

template <typename T>
inline constexpr bool isCountingContainer_v = isCountingContainer<T>::value;

template <typename T>
struct isRenormedHistogram : std::false_type {
};

template <typename source_T>
struct isRenormedHistogram<RenormedHistogram<source_T>> : std::true_type {
};

template <typename T>
inline constexpr bool isRenormedHistogram_v = isRenormedHistogram<T>::value;

template <CoderTag tag_V>
struct SymbolTraits {
  using type = Symbol;
};

template <>
struct SymbolTraits<CoderTag::SingleStream> {
  using type = PrecomputedSymbol;
};

template <CoderTag tag_V>
struct CoderTraits {
};

template <>
struct CoderTraits<CoderTag::Compat> {

  template <size_t lowerBound_V>
  using type = CompatEncoderImpl<lowerBound_V>;
};

template <>
struct CoderTraits<CoderTag::SingleStream> {

  template <size_t lowerBound_V>
  using type = SingleStreamEncoderImpl<lowerBound_V>;
};

template <>
struct CoderTraits<CoderTag::SSE> {

  template <size_t lowerBound_V>
  using type = SSEEncoderImpl<lowerBound_V>;
};

template <>
struct CoderTraits<CoderTag::AVX2> {

  template <size_t lowerBound_V>
  using type = AVXEncoderImpl<lowerBound_V>;
};

} // namespace internal
} // namespace o2::rans

#endif /* RANS_INTERNAL_COMMON_TYPETRAITS_H_ */