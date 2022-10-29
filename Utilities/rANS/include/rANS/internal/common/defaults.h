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
/// @brief  sane compile time defaults for encoders/decoders

#ifndef RANS_INTERNAL_COMMON_DEFAULTS_H_
#define RANS_INTERNAL_COMMON_DEFAULTS_H_

#include <cstdint>
#include <cstring>

#include "rANS/internal/common/defines.h"
#include "rANS/internal/common/typetraits.h"

namespace o2::rans
{

enum class CoderTag : uint8_t { Compat,
                                SingleStream,
                                SSE,
                                AVX2 };

namespace defaults
{

#if defined(RANS_AVX2)
inline constexpr o2::rans::CoderTag DefaultTag = o2::rans::CoderTag::AVX2;
#elif defined(RANS_SSE_ONLY)
inline constexpr o2::rans::CoderTag DefaultTag = o2::rans::CoderTag::SSE;
#elif (defined(RANS_SINGLE_STREAM) && !defined(RANS_SIMD))
inline constexpr o2::rans::CoderTag DefaultTag = o2::rans::CoderTag::SingleStream;
#elif (defined(RANS_COMPAT) && !defined(RANS_SINGLE_STREAM) && !defined(RANS_SIMD))
inline constexpr o2::rans::CoderTag DefaultTag = o2::rans::CoderTag::Compat;
#else
#error your hardware or compiler settings do not support librans
#endif

namespace internal
{
inline constexpr size_t RenormingLowerBound = 20;
}

template <CoderTag tag_V = DefaultTag>
struct EncoderImpl;

template <>
struct EncoderImpl<CoderTag::Compat> {
  inline static constexpr size_t nStreams = 2;
  inline static constexpr size_t renormingLowerBound = internal::RenormingLowerBound;
};

#ifdef RANS_SINGLE_STREAM
template <>
struct EncoderImpl<CoderTag::SingleStream> {
  inline static constexpr size_t nStreams = 2;
  inline static constexpr size_t renormingLowerBound = internal::RenormingLowerBound;
};
#endif

#ifdef RANS_SSE
template <>
struct EncoderImpl<CoderTag::SSE> {
  inline static constexpr size_t nStreams = 16;
  inline static constexpr size_t renormingLowerBound = internal::RenormingLowerBound;
};
#endif

#ifdef RANS_AVX2
template <>
struct EncoderImpl<CoderTag::AVX2> {
  inline static constexpr size_t nStreams = 16;
  inline static constexpr size_t renormingLowerBound = internal::RenormingLowerBound;
};
#endif

} // namespace defaults
} // namespace o2::rans

#endif /* RANS_INTERNAL_COMMON_DEFAULTS_H_ */