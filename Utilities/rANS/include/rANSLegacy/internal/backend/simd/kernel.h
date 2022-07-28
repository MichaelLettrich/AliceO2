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
/// @brief

#ifndef RANSLEGACY_INTERNAL_SIMD_KERNEL_H
#define RANSLEGACY_INTERNAL_SIMD_KERNEL_H

#include <immintrin.h>
#include <cfenv>

#include <array>
#include <tuple>

#include "rANSLegacy/internal/backend/simd/SymbolTable.h"
#include "rANSLegacy/internal/backend/simd/types.h"
#include "rANSLegacy/internal/backend/simd/utils.h"
#include "rANSLegacy/internal/helper.h"
#include "rANSLegacy/internal/Symbol.h"

namespace o2
{
namespace ranslegacy
{
namespace internal
{
namespace simd
{

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline __m128i load(SIMDView<const T, SIMDWidth::SSE, 1> v) noexcept
{
  return _mm_load_si128(reinterpret_cast<const __m128i*>(v.data()));
};

inline __m128d load(pdcV_t<SIMDWidth::SSE> v) noexcept
{
  return _mm_load_pd(v.data());
};

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline __m128i load(const AlignedArray<T, SIMDWidth::SSE>& v) noexcept
{
  return _mm_load_si128(reinterpret_cast<const __m128i*>(v.data()));
};

inline __m128d load(const pd_t<SIMDWidth::SSE>& v) noexcept
{
  return _mm_load_pd(v.data());
};

#ifdef __AVX2__

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline __m256i load(const AlignedArray<T, SIMDWidth::AVX>& v) noexcept
{
  return _mm256_load_si256(reinterpret_cast<const __m256i*>(v.data()));
};

inline __m256d load(const pd_t<SIMDWidth::AVX> v) noexcept
{
  return _mm256_load_pd(v.data());
};

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline __m256i load(SIMDView<const T, SIMDWidth::AVX, 1> v) noexcept
{
  return _mm256_load_si256(reinterpret_cast<const __m256i*>(v.data()));
};

inline __m256d load(pdcV_t<SIMDWidth::AVX> v) noexcept
{
  return _mm256_load_pd(v.data());
};

#endif /* __AVX2__ */

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline AlignedArray<T, SIMDWidth::SSE> store(__m128i inVec) noexcept
{
  AlignedArray<T, SIMDWidth::SSE> out;
  _mm_store_si128(reinterpret_cast<__m128i*>(out.data()), inVec);
  return out;
};

inline AlignedArray<double_t, SIMDWidth::SSE> store(__m128d inVec) noexcept
{
  AlignedArray<double_t, SIMDWidth::SSE> out;
  _mm_store_pd(out.data(), inVec);
  return out;
};

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline void store(__m128i inVec, SIMDView<T, SIMDWidth::SSE, 1, true> v) noexcept
{
  _mm_store_si128(reinterpret_cast<__m128i*>(v.data()), inVec);
};

inline void store(__m128d inVec, SIMDView<double_t, SIMDWidth::SSE, 1, true> v) noexcept
{
  _mm_store_pd(v.data(), inVec);
};

#ifdef __AVX2__

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline AlignedArray<T, SIMDWidth::AVX> store(__m256i inVec) noexcept
{
  AlignedArray<T, SIMDWidth::AVX> out;
  _mm256_store_si256(reinterpret_cast<__m256i*>(out.data()), inVec);
  return out;
};

inline AlignedArray<double_t, SIMDWidth::AVX> store(__m256d inVec) noexcept
{
  AlignedArray<double_t, SIMDWidth::AVX> out;
  _mm256_store_pd(out.data(), inVec);
  return out;
};

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline void store(__m256i inVec, SIMDView<T, SIMDWidth::AVX, 1, true> v) noexcept
{
  _mm256_store_si256(reinterpret_cast<__m256i*>(v.data()), inVec);
};

inline void store(__m256d inVec, SIMDView<double_t, SIMDWidth::AVX, 1, true> v) noexcept
{
  _mm256_store_pd(v.data(), inVec);
};

#endif /* __AVX2__ */

template <SIMDWidth width_V>
inline decltype(auto) setAll(uint64_t value) noexcept
{
  if constexpr (width_V == SIMDWidth::SSE) {
    return _mm_set1_epi64x(value);
  } else {
    return _mm256_set1_epi64x(value);
  }
};

template <SIMDWidth width_V>
inline decltype(auto) setAll(uint32_t value) noexcept
{
  if constexpr (width_V == SIMDWidth::SSE) {
    return _mm_set1_epi32(value);
  } else {
    return _mm256_set1_epi32(value);
  }
};

template <SIMDWidth width_V>
inline decltype(auto) setAll(uint16_t value) noexcept
{
  if constexpr (width_V == SIMDWidth::SSE) {
    return _mm_set1_epi16(value);
  } else {
    return _mm256_set1_epi16(value);
  }
};

template <SIMDWidth width_V>
inline decltype(auto) setAll(uint8_t value) noexcept
{
  if constexpr (width_V == SIMDWidth::SSE) {
    return _mm_set1_epi8(value);
  } else {
    return _mm256_set1_epi8(value);
  }
};

template <SIMDWidth width_V>
inline decltype(auto) setAll(double_t value) noexcept
{
  if constexpr (width_V == SIMDWidth::SSE) {
    return _mm_set1_pd(value);
  } else {
    return _mm256_set1_pd(value);
  }
};

//
// uint32 -> double
//
template <SIMDWidth width_V>
inline auto int32ToDouble(__m128i in) noexcept
{
  if constexpr (width_V == SIMDWidth::SSE) {
    return _mm_cvtepi32_pd(in);
  } else if constexpr (width_V == SIMDWidth::AVX) {
#ifdef __AVX2__
    return _mm256_cvtepi32_pd(in);
#endif /* __AVX2__ */
  }
};

//
// uint64 -> double
// Only works for inputs in the range: [0, 2^52)
//

inline __m128d uint64ToDouble(__m128i in) noexcept
{
#if !defined(NDEBUG)
  auto vec = store<uint64_t>(in);
  for (auto i : vec) {
    assert(i < pow2(52));
  }
#endif
  in = _mm_or_si128(in, _mm_castpd_si128(_mm_set1_pd(detail::AlignMantissaMagic)));
  __m128d out = _mm_sub_pd(_mm_castsi128_pd(in), _mm_set1_pd(detail::AlignMantissaMagic));
  return out;
};

#ifdef __AVX2__
//
// uint64 -> double
// Only works for inputs in the range: [0, 2^52)
//
inline __m256d uint64ToDouble(__m256i in) noexcept
{
#if !defined(NDEBUG)
  auto vec = store<uint64_t>(in);
  for (auto i : vec) {
    assert(i < pow2(52));
  }
#endif
  in = _mm256_or_si256(in, _mm256_castpd_si256(_mm256_set1_pd(detail::AlignMantissaMagic)));
  __m256d out = _mm256_sub_pd(_mm256_castsi256_pd(in), _mm256_set1_pd(detail::AlignMantissaMagic));
  return out;
};
#endif /* __AVX2__ */

inline __m128i doubleToUint64(__m128d in) noexcept
{
#if !defined(NDEBUG)
  auto vec = store(in);
  for (auto i : vec) {
    assert(i < pow2(52));
  }
#endif
  in = _mm_add_pd(in, _mm_set1_pd(detail::AlignMantissaMagic));
  __m128i out = _mm_xor_si128(_mm_castpd_si128(in),
                              _mm_castpd_si128(_mm_set1_pd(detail::AlignMantissaMagic)));
  return out;
}

#ifdef __AVX2__

inline __m256i doubleToUint64(__m256d in) noexcept
{
#if !defined(NDEBUG)
  auto vec = store(in);
  for (auto i : vec) {
    assert(i < pow2(52));
  }
#endif

  in = _mm256_add_pd(in, _mm256_set1_pd(detail::AlignMantissaMagic));
  __m256i out = _mm256_xor_si256(_mm256_castpd_si256(in),
                                 _mm256_castpd_si256(_mm256_set1_pd(detail::AlignMantissaMagic)));
  return out;
}

#endif /* __AVX2__ */

template <SIMDWidth>
struct DivMod;

template <>
struct DivMod<SIMDWidth::SSE> {
  __m128d div;
  __m128d mod;
};

// calculate both floor(a/b) and a%b
inline DivMod<SIMDWidth::SSE>
  divMod(__m128d numerator, __m128d denominator) noexcept
{
  __m128d div = _mm_floor_pd(_mm_div_pd(numerator, denominator));
  __m128d mod = _mm_fnmadd_pd(div, denominator, numerator);
  return {div, mod};
}

#ifdef __AVX2__

template <>
struct DivMod<SIMDWidth::AVX> {
  __m256d div;
  __m256d mod;
};

// calculate both floor(a/b) and a%b
inline DivMod<SIMDWidth::AVX> divMod(__m256d numerator, __m256d denominator) noexcept
{
  __m256d div = _mm256_floor_pd(_mm256_div_pd(numerator, denominator));
  __m256d mod = _mm256_fnmadd_pd(div, denominator, numerator);
  return {div, mod};
}
#endif /* __AVX2__ */

//
// rans Encode
//
inline __m128i ransEncode(__m128i state, __m128d frequency, __m128d cumulative, __m128d normalization) noexcept
{
#if !defined(NDEBUG)
  auto vec = store<uint64_t>(state);
  for (auto i : vec) {
    assert(i < pow2(52));
  }
#endif

  auto [div, mod] = divMod(uint64ToDouble(state), frequency);
  auto newState = _mm_fmadd_pd(normalization, div, cumulative);
  newState = _mm_add_pd(newState, mod);

  return doubleToUint64(newState);
};
#ifdef __AVX2__

//
// rans Encode
//
inline __m256i ransEncode(__m256i state, __m256d frequency, __m256d cumulative, __m256d normalization) noexcept
{
#if !defined(NDEBUG)
  auto vec = store<uint64_t>(state);
  for (auto i : vec) {
    assert(i < pow2(52));
  }
#endif

  auto [div, mod] = divMod(uint64ToDouble(state), frequency);
  auto newState = _mm256_fmadd_pd(normalization, div, cumulative);
  newState = _mm256_add_pd(newState, mod);

  return doubleToUint64(newState);
};

#endif /* __AVX2__ */

inline void aosToSoa(ArrayView<const Symbol*, 2> in, __m128i* __restrict__ frequency, __m128i* __restrict__ cumulatedFrequency) noexcept
{
  __m128i in0Reg = _mm_loadu_si128(reinterpret_cast<__m128i const*>(in[0]->data()));
  __m128i in1Reg = _mm_loadu_si128(reinterpret_cast<__m128i const*>(in[1]->data()));

  *frequency = _mm_unpacklo_epi32(in0Reg, in1Reg);
  *cumulatedFrequency = _mm_shuffle_epi32(*frequency, _MM_SHUFFLE(0, 0, 3, 2));
}

inline void aosToSoa(ArrayView<const Symbol*, 4> in, __m128i* __restrict__ frequency, __m128i* __restrict__ cumulatedFrequency) noexcept
{
  __m128i in0Reg = _mm_loadu_si128(reinterpret_cast<__m128i const*>(in[0]->data()));
  __m128i in1Reg = _mm_loadu_si128(reinterpret_cast<__m128i const*>(in[1]->data()));
  __m128i in2Reg = _mm_loadu_si128(reinterpret_cast<__m128i const*>(in[2]->data()));
  __m128i in3Reg = _mm_loadu_si128(reinterpret_cast<__m128i const*>(in[3]->data()));

  __m128i merged0Reg = _mm_unpacklo_epi32(in0Reg, in1Reg);
  __m128i merged1Reg = _mm_unpacklo_epi32(in2Reg, in3Reg);
  *frequency = _mm_unpacklo_epi64(merged0Reg, merged1Reg);
  *cumulatedFrequency = _mm_unpackhi_epi64(merged0Reg, merged1Reg);
};

inline __m128i cmpgeq_epi64(__m128i a, __m128i b) noexcept
{
  __m128i cmpGreater = _mm_cmpgt_epi64(a, b);
  __m128i cmpEqual = _mm_cmpeq_epi64(a, b);
  return _mm_or_si128(cmpGreater, cmpEqual);
};

#ifdef __AVX2__
inline __m256i cmpgeq_epi64(__m256i a, __m256i b) noexcept
{
  __m256i cmpGreater = _mm256_cmpgt_epi64(a, b);
  __m256i cmpEqual = _mm256_cmpeq_epi64(a, b);
  return _mm256_or_si256(cmpGreater, cmpEqual);
};
#endif /* __AVX2__ */

template <SIMDWidth width_V, uint64_t lowerBound_V, uint8_t streamBits_V>
inline auto computeMaxState(__m128i frequencyVec, uint8_t symbolTablePrecisionBits) noexcept
{
  const uint64_t xmax = (lowerBound_V >> symbolTablePrecisionBits) << streamBits_V;
  const uint8_t shift = log2UIntNZ(xmax);
  if constexpr (width_V == SIMDWidth::SSE) {
    __m128i frequencyVecEpi64 = _mm_cvtepi32_epi64(frequencyVec);
    return _mm_slli_epi64(frequencyVecEpi64, shift);
  }
  if constexpr (width_V == SIMDWidth::AVX) {
#ifdef __AVX2__
    __m256i frequencyVecEpi64 = _mm256_cvtepi32_epi64(frequencyVec);
    return _mm256_slli_epi64(frequencyVecEpi64, shift);
#endif /* __AVX2__ */
  }
};

template <uint8_t streamBits_V>
inline __m128i computeNewState(__m128i stateVec, __m128i cmpVec) noexcept
{
  //newState = (state >= maxState) ? state >> streamBits_V : state
  __m128i newStateVec = _mm_srli_epi64(stateVec, streamBits_V);
  newStateVec = _mm_blendv_epi8(stateVec, newStateVec, cmpVec);
  return newStateVec;
};

#ifdef __AVX2__
template <uint8_t streamBits_V>
inline __m256i computeNewState(__m256i stateVec, __m256i cmpVec) noexcept
{
  //newState = (state >= maxState) ? state >> streamBits_V : state
  __m256i newStateVec = _mm256_srli_epi64(stateVec, streamBits_V);
  newStateVec = _mm256_blendv_epi8(stateVec, newStateVec, cmpVec);
  return newStateVec;
};

#endif /* __AVX2__ */

inline constexpr std::array<epi8_t<SIMDWidth::SSE>, 16>
  SSEStreamOutLUT{{
    {0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0000   0xFFFFu, //0b0000
    {0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0001   0x0FFFu, //0b0001
    {0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0010   0x1FFFu, //0b0010
    {0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0011   0x10FFu, //0b0011
    {0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0100   0x2FFFu, //0b0100
    {0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0101   0x20FFu, //0b0101
    {0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0110   0x12FFu, //0b0110
    {0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0111   0x120Fu, //0b0111
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b1000   0x3FFFu, //0b1000
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b1001   0x30FFu, //0b1001
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b1010   0x31FFu, //0b1010
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b1011   0x310Fu, //0b1011
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b1100   0x32FFu, //0b1100
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b1101   0x320Fu, //0b1101
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b1110   0x312Fu, //0b1110
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8}  //0b1111   0x3120u, //0b1111
  }};

inline constexpr std::array<uint32_t, 256> AVXStreamOutLUT{
  0xFFFFFFFFu, //0b00000000
  0x0FFFFFFFu, //0b00000001
  0x1FFFFFFFu, //0b00000010
  0x10FFFFFFu, //0b00000011
  0x2FFFFFFFu, //0b00000100
  0x20FFFFFFu, //0b00000101
  0x12FFFFFFu, //0b00000110
  0x120FFFFFu, //0b00000111
  0x3FFFFFFFu, //0b00001000
  0x30FFFFFFu, //0b00001001
  0x31FFFFFFu, //0b00001010
  0x310FFFFFu, //0b00001011
  0x32FFFFFFu, //0b00001100
  0x320FFFFFu, //0b00001101
  0x312FFFFFu, //0b00001110
  0x3120FFFFu, //0b00001111
  0x4FFFFFFFu, //0b00010000
  0x40FFFFFFu, //0b00010001
  0x14FFFFFFu, //0b00010010
  0x140FFFFFu, //0b00010011
  0x42FFFFFFu, //0b00010100
  0x420FFFFFu, //0b00010101
  0x142FFFFFu, //0b00010110
  0x1420FFFFu, //0b00010111
  0x34FFFFFFu, //0b00011000
  0x340FFFFFu, //0b00011001
  0x314FFFFFu, //0b00011010
  0x3140FFFFu, //0b00011011
  0x342FFFFFu, //0b00011100
  0x3420FFFFu, //0b00011101
  0x3142FFFFu, //0b00011110
  0x31420FFFu, //0b00011111
  0x5FFFFFFFu, //0b00100000
  0x50FFFFFFu, //0b00100001
  0x51FFFFFFu, //0b00100010
  0x510FFFFFu, //0b00100011
  0x52FFFFFFu, //0b00100100
  0x520FFFFFu, //0b00100101
  0x512FFFFFu, //0b00100110
  0x5120FFFFu, //0b00100111
  0x53FFFFFFu, //0b00101000
  0x530FFFFFu, //0b00101001
  0x531FFFFFu, //0b00101010
  0x5310FFFFu, //0b00101011
  0x532FFFFFu, //0b00101100
  0x5320FFFFu, //0b00101101
  0x5312FFFFu, //0b00101110
  0x53120FFFu, //0b00101111
  0x54FFFFFFu, //0b00110000
  0x540FFFFFu, //0b00110001
  0x514FFFFFu, //0b00110010
  0x5140FFFFu, //0b00110011
  0x542FFFFFu, //0b00110100
  0x5420FFFFu, //0b00110101
  0x5142FFFFu, //0b00110110
  0x51420FFFu, //0b00110111
  0x534FFFFFu, //0b00111000
  0x5340FFFFu, //0b00111001
  0x5314FFFFu, //0b00111010
  0x53140FFFu, //0b00111011
  0x5342FFFFu, //0b00111100
  0x53420FFFu, //0b00111101
  0x53142FFFu, //0b00111110
  0x531420FFu, //0b00111111
  0x6FFFFFFFu, //0b01000000
  0x60FFFFFFu, //0b01000001
  0x16FFFFFFu, //0b01000010
  0x160FFFFFu, //0b01000011
  0x62FFFFFFu, //0b01000100
  0x620FFFFFu, //0b01000101
  0x162FFFFFu, //0b01000110
  0x1620FFFFu, //0b01000111
  0x36FFFFFFu, //0b01001000
  0x360FFFFFu, //0b01001001
  0x316FFFFFu, //0b01001010
  0x3160FFFFu, //0b01001011
  0x362FFFFFu, //0b01001100
  0x3620FFFFu, //0b01001101
  0x3162FFFFu, //0b01001110
  0x31620FFFu, //0b01001111
  0x64FFFFFFu, //0b01010000
  0x640FFFFFu, //0b01010001
  0x164FFFFFu, //0b01010010
  0x1640FFFFu, //0b01010011
  0x642FFFFFu, //0b01010100
  0x6420FFFFu, //0b01010101
  0x1642FFFFu, //0b01010110
  0x16420FFFu, //0b01010111
  0x364FFFFFu, //0b01011000
  0x3640FFFFu, //0b01011001
  0x3164FFFFu, //0b01011010
  0x31640FFFu, //0b01011011
  0x3642FFFFu, //0b01011100
  0x36420FFFu, //0b01011101
  0x31642FFFu, //0b01011110
  0x316420FFu, //0b01011111
  0x56FFFFFFu, //0b01100000
  0x560FFFFFu, //0b01100001
  0x516FFFFFu, //0b01100010
  0x5160FFFFu, //0b01100011
  0x562FFFFFu, //0b01100100
  0x5620FFFFu, //0b01100101
  0x5162FFFFu, //0b01100110
  0x51620FFFu, //0b01100111
  0x536FFFFFu, //0b01101000
  0x5360FFFFu, //0b01101001
  0x5316FFFFu, //0b01101010
  0x53160FFFu, //0b01101011
  0x5362FFFFu, //0b01101100
  0x53620FFFu, //0b01101101
  0x53162FFFu, //0b01101110
  0x531620FFu, //0b01101111
  0x564FFFFFu, //0b01110000
  0x5640FFFFu, //0b01110001
  0x5164FFFFu, //0b01110010
  0x51640FFFu, //0b01110011
  0x5642FFFFu, //0b01110100
  0x56420FFFu, //0b01110101
  0x51642FFFu, //0b01110110
  0x516420FFu, //0b01110111
  0x5364FFFFu, //0b01111000
  0x53640FFFu, //0b01111001
  0x53164FFFu, //0b01111010
  0x531640FFu, //0b01111011
  0x53642FFFu, //0b01111100
  0x536420FFu, //0b01111101
  0x531642FFu, //0b01111110
  0x5316420Fu, //0b01111111
  0x7FFFFFFFu, //0b10000000
  0x70FFFFFFu, //0b10000001
  0x71FFFFFFu, //0b10000010
  0x710FFFFFu, //0b10000011
  0x72FFFFFFu, //0b10000100
  0x720FFFFFu, //0b10000101
  0x712FFFFFu, //0b10000110
  0x7120FFFFu, //0b10000111
  0x73FFFFFFu, //0b10001000
  0x730FFFFFu, //0b10001001
  0x731FFFFFu, //0b10001010
  0x7310FFFFu, //0b10001011
  0x732FFFFFu, //0b10001100
  0x7320FFFFu, //0b10001101
  0x7312FFFFu, //0b10001110
  0x73120FFFu, //0b10001111
  0x74FFFFFFu, //0b10010000
  0x740FFFFFu, //0b10010001
  0x714FFFFFu, //0b10010010
  0x7140FFFFu, //0b10010011
  0x742FFFFFu, //0b10010100
  0x7420FFFFu, //0b10010101
  0x7142FFFFu, //0b10010110
  0x71420FFFu, //0b10010111
  0x734FFFFFu, //0b10011000
  0x7340FFFFu, //0b10011001
  0x7314FFFFu, //0b10011010
  0x73140FFFu, //0b10011011
  0x7342FFFFu, //0b10011100
  0x73420FFFu, //0b10011101
  0x73142FFFu, //0b10011110
  0x731420FFu, //0b10011111
  0x75FFFFFFu, //0b10100000
  0x750FFFFFu, //0b10100001
  0x751FFFFFu, //0b10100010
  0x7510FFFFu, //0b10100011
  0x752FFFFFu, //0b10100100
  0x7520FFFFu, //0b10100101
  0x7512FFFFu, //0b10100110
  0x75120FFFu, //0b10100111
  0x753FFFFFu, //0b10101000
  0x7530FFFFu, //0b10101001
  0x7531FFFFu, //0b10101010
  0x75310FFFu, //0b10101011
  0x7532FFFFu, //0b10101100
  0x75320FFFu, //0b10101101
  0x75312FFFu, //0b10101110
  0x753120FFu, //0b10101111
  0x754FFFFFu, //0b10110000
  0x7540FFFFu, //0b10110001
  0x7514FFFFu, //0b10110010
  0x75140FFFu, //0b10110011
  0x7542FFFFu, //0b10110100
  0x75420FFFu, //0b10110101
  0x75142FFFu, //0b10110110
  0x751420FFu, //0b10110111
  0x7534FFFFu, //0b10111000
  0x75340FFFu, //0b10111001
  0x75314FFFu, //0b10111010
  0x753140FFu, //0b10111011
  0x75342FFFu, //0b10111100
  0x753420FFu, //0b10111101
  0x753142FFu, //0b10111110
  0x7531420Fu, //0b10111111
  0x76FFFFFFu, //0b11000000
  0x760FFFFFu, //0b11000001
  0x716FFFFFu, //0b11000010
  0x7160FFFFu, //0b11000011
  0x762FFFFFu, //0b11000100
  0x7620FFFFu, //0b11000101
  0x7162FFFFu, //0b11000110
  0x71620FFFu, //0b11000111
  0x736FFFFFu, //0b11001000
  0x7360FFFFu, //0b11001001
  0x7316FFFFu, //0b11001010
  0x73160FFFu, //0b11001011
  0x7362FFFFu, //0b11001100
  0x73620FFFu, //0b11001101
  0x73162FFFu, //0b11001110
  0x731620FFu, //0b11001111
  0x764FFFFFu, //0b11010000
  0x7640FFFFu, //0b11010001
  0x7164FFFFu, //0b11010010
  0x71640FFFu, //0b11010011
  0x7642FFFFu, //0b11010100
  0x76420FFFu, //0b11010101
  0x71642FFFu, //0b11010110
  0x716420FFu, //0b11010111
  0x7364FFFFu, //0b11011000
  0x73640FFFu, //0b11011001
  0x73164FFFu, //0b11011010
  0x731640FFu, //0b11011011
  0x73642FFFu, //0b11011100
  0x736420FFu, //0b11011101
  0x731642FFu, //0b11011110
  0x7316420Fu, //0b11011111
  0x756FFFFFu, //0b11100000
  0x7560FFFFu, //0b11100001
  0x7516FFFFu, //0b11100010
  0x75160FFFu, //0b11100011
  0x7562FFFFu, //0b11100100
  0x75620FFFu, //0b11100101
  0x75162FFFu, //0b11100110
  0x751620FFu, //0b11100111
  0x7536FFFFu, //0b11101000
  0x75360FFFu, //0b11101001
  0x75316FFFu, //0b11101010
  0x753160FFu, //0b11101011
  0x75362FFFu, //0b11101100
  0x753620FFu, //0b11101101
  0x753162FFu, //0b11101110
  0x7531620Fu, //0b11101111
  0x7564FFFFu, //0b11110000
  0x75640FFFu, //0b11110001
  0x75164FFFu, //0b11110010
  0x751640FFu, //0b11110011
  0x75642FFFu, //0b11110100
  0x756420FFu, //0b11110101
  0x751642FFu, //0b11110110
  0x7516420Fu, //0b11110111
  0x75364FFFu, //0b11111000
  0x753640FFu, //0b11111001
  0x753164FFu, //0b11111010
  0x7531640Fu, //0b11111011
  0x753642FFu, //0b11111100
  0x7536420Fu, //0b11111101
  0x7531642Fu, //0b11111110
  0x75316420u  //0b11111111
};

template <SIMDWidth>
struct StreamOutResult;

template <>
struct StreamOutResult<SIMDWidth::SSE> {
  uint32_t nElemens;
  __m128i streamOutVec;
};

// inline StreamOutResult<SIMDWidth::SSE> streamOut(const __m128i* __restrict__ stateVec, const __m128i* __restrict__ cmpVec) noexcept
// {
//   //  std::cout << "streamOut\n";
//   //  std::cout << asHex(store<uint32_t>(stateVec[0])) << asHex(store<uint32_t>(stateVec[1])) << "\n";
//   //  std::cout << asHex(store<uint32_t>(cmpVec[0])) << asHex(store<uint32_t>(cmpVec[1])) << "\n";
//   //shift by one
//   auto shifted1 = _mm_slli_epi64(stateVec[1], 32);
//   //  std::cout << "stateVec[1] "  << asHex(store<uint32_t>(shifted1)) << "\n";

//   __m128i statesFused = _mm_blend_epi32(stateVec[0], shifted1, 0b1010);
//   __m128i cmpFused = _mm_blend_epi32(cmpVec[0], cmpVec[1], 0b1010);
//   //  std::cout <<"cmpFused "<< asHex(store<uint32_t>(cmpFused)) << "\n";
//   //  std::cout <<"statesFused "<< asHex(store<uint32_t>(statesFused)) << "\n";
//   statesFused = _mm_and_si128(statesFused, cmpFused);
//   //  std::cout <<"statesFused "<< asHex(store<uint32_t>(statesFused)) << "\n";
//   const uint32_t id = _mm_movemask_ps(_mm_castsi128_ps(cmpFused));

//   __m128i permutationMask = _mm_set1_epi32(AVXInterleavedPermutationLUT[id]);
//   //  std::cout <<"permutationMask "<< asHex(store<uint32_t>(permutationMask)) << "\n";
//   constexpr epi32_t<SIMDWidth::SSE> mask{0xF0000000u, 0x0F000000u, 0x00F00000u, 0x000F0000u};
//   permutationMask = _mm_and_si128(permutationMask, load(toConstSIMDView(mask)));
//   //  std::cout <<"permutationMask "<< asHex(store<uint32_t>(permutationMask)) << "\n";

//   constexpr epi32_t<SIMDWidth::SSE> shift{28u, 24u, 20u, 16u};
//   permutationMask = _mm_srlv_epi32(permutationMask, load(toConstSIMDView(shift)));
//   //  std::cout <<"permutationMask "<< asHex(store<uint32_t>(permutationMask)) << "\n";
//   auto streamOutVec = _mm_shuffle_epi8(statesFused, permutationMask);

//   //  std::cout << "streamOut end\n";
//   return {static_cast<uint32_t>(_mm_popcnt_u32(id)), streamOutVec};
// };

inline StreamOutResult<SIMDWidth::SSE> streamOut(const __m128i* __restrict__ stateVec, const __m128i* __restrict__ cmpVec) noexcept
{
  // std::cout << "streamOut\n";
  // std::cout << asHex(store<uint32_t>(stateVec[0])) << asHex(store<uint32_t>(stateVec[1])) << "\n";
  // std::cout << asHex(store<uint32_t>(cmpVec[0])) << asHex(store<uint32_t>(cmpVec[1])) << "\n";
  //shift by one
  auto shifted1 = _mm_slli_epi64(stateVec[1], 32);
  // std::cout << "stateVec[1] " << asHex(store<uint32_t>(shifted1)) << "\n";

  __m128i statesFused = _mm_blend_epi32(stateVec[0], shifted1, 0b1010);
  __m128i cmpFused = _mm_blend_epi32(cmpVec[0], cmpVec[1], 0b1010);
  // std::cout << "cmpFused " << asHex(store<uint32_t>(cmpFused)) << "\n";
  // std::cout << "statesFused " << asHex(store<uint32_t>(statesFused)) << "\n";
  const uint32_t id = _mm_movemask_ps(_mm_castsi128_ps(cmpFused));
  // std::cout << fmt::format("id: {:#0b}\n", id);

  __m128i permutationMask = load(toConstSIMDView(SSEStreamOutLUT[id]));
  // std::cout << "permutationMask " << asHex(store<uint32_t>(permutationMask)) << "\n";
  auto streamOutVec = _mm_shuffle_epi8(statesFused, permutationMask);

  // std::cout << "streamOut end\n";
  return {static_cast<uint32_t>(_mm_popcnt_u32(id)), streamOutVec};
};

#ifdef __AVX2__
template <>
struct StreamOutResult<SIMDWidth::AVX> {
  uint32_t nElemens;
  __m256i streamOutVec;
};

// inline constexpr epi32_t<SIMDWidth::AVX> PermuteAVX{0x7u, 0x5u, 0x3u, 0x1u, 0x6u, 0x4u, 0x2u, 0x0u};

// inline StreamOutResult<SIMDWidth::AVX> streamOut(const __m256i* __restrict__ stateVec, const __m256i* __restrict__ cmpVec) noexcept
// {
//   // std::cout << "streamOut\n";
//   // std::cout << asHex(store<uint32_t>(stateVec[0])) << asHex(store<uint32_t>(stateVec[1])) << "\n";
//   // std::cout << asHex(store<uint32_t>(cmpVec[0])) << asHex(store<uint32_t>(cmpVec[1])) << "\n";
//   //shift by one
//   auto shifted1 = _mm256_slli_epi64(stateVec[1], 32);
//   // std::cout << "stateVec[1] " << asHex(store<uint32_t>(shifted1)) << "\n";

//   __m256i statesFused = _mm256_blend_epi32(stateVec[0], shifted1, 0b10101010);
//   __m256i cmpFused = _mm256_blend_epi32(cmpVec[0], cmpVec[1], 0b10101010);
//   // std::cout << "cmpFused " << asHex(store<uint32_t>(cmpFused)) << "\n";
//   // std::cout << "statesFused " << asHex(store<uint32_t>(statesFused)) << "\n";

//   statesFused = _mm256_permutevar8x32_epi32(statesFused, load(toConstSIMDView(PermuteAVX)));
//   // std::cout << "statesFused2 " << asHex(store<uint32_t>(statesFused2)) << "\n";
//   cmpFused = _mm256_permutevar8x32_epi32(cmpFused, load(toConstSIMDView(PermuteAVX)));
//   // std::cout << "cmpFused2 " << asHex(store<uint32_t>(cmpFused2)) << "\n";
//   statesFused = _mm256_and_si256(statesFused, cmpFused);
//   const uint32_t mask = _mm256_movemask_ps(_mm256_castsi256_ps(cmpFused));
//   // std::cout << fmt::format("mask: {:#0b}\n", mask);
//   uint64_t expanded_mask = _pdep_u64(mask, 0x0101010101010101); // unpack each bit to a byte
//   expanded_mask *= 0xFFU;                                       // mask |= mask<<1 | mask<<2 | ... | mask<<7;
//   // ABC... -> AAAAAAAABBBBBBBBCCCCCCCC...: replicate each bit to fill its byte

//   const uint64_t identity_indices = 0x0706050403020100; // the identity shuffle for vpermps, packed to one index per byte
//   uint64_t wanted_indices = _pext_u64(identity_indices, expanded_mask);

//   __m128i bytevec = _mm_cvtsi64_si128(wanted_indices);
//   __m256i shufmask = _mm256_cvtepu8_epi32(bytevec);

//   auto streamOutVec = _mm256_permutevar8x32_epi32(statesFused, shufmask);
//   // std::cout << "permuted " << asHex(store<uint32_t>(streamOutVec2)) << "\n";

//   // statesFused = _mm256_and_si256(statesFused, cmpFused);
//   // std::cout << "statesFused " << asHex(store<uint32_t>(statesFused)) << "\n";
//   // const uint32_t id = _mm256_movemask_ps(_mm256_castsi256_ps(cmpFused));

//   // __m256i permutationMask = _mm256_set1_epi32(AVXInterleavedPermutationLUT[id]);
//   // std::cout << "permutationMask " << asHex(store<uint32_t>(permutationMask)) << "\n";
//   // constexpr epi32_t<SIMDWidth::AVX> mask{0xF0000000u, 0x0F000000u, 0x00F00000u, 0x000F0000u, 0x0000F000u, 0x00000F00u, 0x000000F0u, 0x0000000Fu};
//   // permutationMask = _mm256_and_si256(permutationMask, load(mask));
//   // std::cout << "permutationMask " << asHex(store<uint32_t>(permutationMask)) << "\n";
//   // constexpr epi32_t<SIMDWidth::AVX> shift{28u, 24u, 20u, 16u, 12u, 8u, 4u, 0u};
//   // permutationMask = _mm256_srlv_epi32(permutationMask, load(shift));
//   // std::cout << "permutationMask " << asHex(store<uint32_t>(permutationMask)) << "\n";
//   // auto streamOutVec = _mm256_permutevar8x32_epi32(statesFused, permutationMask);

//   // std::cout << "streamOut end\n";
//   return {static_cast<uint32_t>(_mm_popcnt_u32(mask)), streamOutVec};
// };

inline StreamOutResult<SIMDWidth::AVX> streamOut(const __m256i* __restrict__ stateVec, const __m256i* __restrict__ cmpVec) noexcept
{
  //  std::cout << "streamOut\n";
  //  std::cout << asHex(store<uint32_t>(stateVec[0])) << asHex(store<uint32_t>(stateVec[1])) << "\n";
  //  std::cout << asHex(store<uint32_t>(cmpVec[0])) << asHex(store<uint32_t>(cmpVec[1])) << "\n";
  //shift by one
  auto shifted1 = _mm256_slli_epi64(stateVec[1], 32);
  //  std::cout << "stateVec[1] "  << asHex(store<uint32_t>(shifted1)) << "\n";

  __m256i statesFused = _mm256_blend_epi32(stateVec[0], shifted1, 0b10101010);
  __m256i cmpFused = _mm256_blend_epi32(cmpVec[0], cmpVec[1], 0b10101010);
  //  std::cout <<"cmpFused "<< asHex(store<uint32_t>(cmpFused)) << "\n";
  //  std::cout <<"statesFused "<< asHex(store<uint32_t>(statesFused)) << "\n";
  statesFused = _mm256_and_si256(statesFused, cmpFused);
  //  std::cout <<"statesFused "<< asHex(store<uint32_t>(statesFused)) << "\n";
  const uint32_t id = _mm256_movemask_ps(_mm256_castsi256_ps(cmpFused));

  __m256i permutationMask = _mm256_set1_epi32(AVXStreamOutLUT[id]);
  //  std::cout <<"permutationMask "<< asHex(store<uint32_t>(permutationMask)) << "\n";
  constexpr epi32_t<SIMDWidth::AVX> mask{0xF0000000u, 0x0F000000u, 0x00F00000u, 0x000F0000u, 0x0000F000u, 0x00000F00u, 0x000000F0u, 0x0000000Fu};
  permutationMask = _mm256_and_si256(permutationMask, load(toConstSIMDView(mask)));
  //  std::cout <<"permutationMask "<< asHex(store<uint32_t>(permutationMask)) << "\n";
  constexpr epi32_t<SIMDWidth::AVX> shift{28u, 24u, 20u, 16u, 12u, 8u, 4u, 0u};
  permutationMask = _mm256_srlv_epi32(permutationMask, load(toConstSIMDView(shift)));
  //  std::cout <<"permutationMask "<< asHex(store<uint32_t>(permutationMask)) << "\n";
  auto streamOutVec = _mm256_permutevar8x32_epi32(statesFused, permutationMask);

  //  std::cout << "streamOut end\n";
  return {static_cast<uint32_t>(_mm_popcnt_u32(id)), streamOutVec};
};

#endif /* __AVX2__ */

template <SIMDWidth, typename output_IT>
struct RenormResult;

template <typename output_IT>
struct RenormResult<SIMDWidth::SSE, output_IT> {
  output_IT outputIter;
  __m128i newState;
};

#ifdef __AVX2__
template <typename output_IT>
struct RenormResult<SIMDWidth::AVX, output_IT> {
  output_IT outputIter;
  __m256i newState;
};
#endif /* __AVX2__ */

template <typename output_IT, uint64_t lowerBound_V, uint8_t streamBits_V>
inline output_IT ransRenorm(const __m128i* __restrict__ state, const __m128i* __restrict__ frequency, uint8_t symbolTablePrecisionBits, output_IT outputIter, __m128i* __restrict__ newState) noexcept
{
  __m128i maxState[2];
  __m128i cmp[2];

  // calculate maximum state
  maxState[0] = computeMaxState<SIMDWidth::SSE, lowerBound_V, streamBits_V>(frequency[0], symbolTablePrecisionBits);
  maxState[1] = computeMaxState<SIMDWidth::SSE, lowerBound_V, streamBits_V>(frequency[1], symbolTablePrecisionBits);
  //cmp = (state >= maxState)
  cmp[0] = cmpgeq_epi64(state[0], maxState[0]);
  cmp[1] = cmpgeq_epi64(state[1], maxState[1]);
  //newState = (state >= maxState) ? state >> streamBits_V : state
  newState[0] = computeNewState<streamBits_V>(state[0], cmp[0]);
  newState[1] = computeNewState<streamBits_V>(state[1], cmp[1]);

  auto [nStreamOutWords, streamOutResult] = streamOut(state, cmp);
  if constexpr (std::is_pointer_v<output_IT>) {
    _mm_storeu_si128(reinterpret_cast<__m128i*>(outputIter + 1), streamOutResult);
    outputIter += nStreamOutWords;
  } else {
    auto result = store<uint32_t>(streamOutResult);
    for (size_t i = 0; i < nStreamOutWords; ++i) {
      *(++outputIter) = result[i];
    }
  }

  return outputIter;
};

#ifdef __AVX2__
template <typename output_IT, uint64_t lowerBound_V, uint8_t streamBits_V>
inline output_IT ransRenorm(const __m256i* __restrict__ state, const __m128i* __restrict__ frequency, uint8_t symbolTablePrecisionBits, output_IT outputIter, __m256i* __restrict__ newState) noexcept
{
  __m256i maxState[2];
  __m256i cmp[2];

  // calculate maximum state
  maxState[0] = computeMaxState<SIMDWidth::AVX, lowerBound_V, streamBits_V>(frequency[0], symbolTablePrecisionBits);
  maxState[1] = computeMaxState<SIMDWidth::AVX, lowerBound_V, streamBits_V>(frequency[1], symbolTablePrecisionBits);
  //cmp = (state >= maxState)
  cmp[0] = cmpgeq_epi64(state[0], maxState[0]);
  cmp[1] = cmpgeq_epi64(state[1], maxState[1]);
  //newState = (state >= maxState) ? state >> streamBits_V : state
  newState[0] = computeNewState<streamBits_V>(state[0], cmp[0]);
  newState[1] = computeNewState<streamBits_V>(state[1], cmp[1]);

  auto [nStreamOutWords, streamOutResult] = streamOut(state, cmp);
  if constexpr (std::is_pointer_v<output_IT>) {
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(outputIter + 1), streamOutResult);
    outputIter += nStreamOutWords;
  } else {
    auto result = store<uint32_t>(streamOutResult);
    for (size_t i = 0; i < nStreamOutWords; ++i) {
      *(++outputIter) = result[i];
    }
  }

  return outputIter;
};
#endif /* __AVX2__ */

template <typename source_IT>
inline const Symbol* lookupSymbol(source_IT iter, const simd::SymbolTable& symbolTable, std::vector<typename std::iterator_traits<source_IT>::value_type>& literals) noexcept
{
  const auto symbol = *iter;
  const auto* encoderSymbol = &(symbolTable[symbol]);
  if (symbolTable.isEscapeSymbol(*encoderSymbol)) {
    literals.push_back(symbol);
  }
  return encoderSymbol;
};

struct UnrolledSymbols {
  __m128i frequencies[2];
  __m128i cumulativeFrequencies[2];
};

template <typename source_IT, SIMDWidth width_V>
std::tuple<source_IT, UnrolledSymbols> getSymbols(source_IT symbolIter, const SymbolTable& symbolTable, std::vector<typename std::iterator_traits<source_IT>::value_type>& literals)
{
  UnrolledSymbols unrolledSymbols;

  if constexpr (width_V == SIMDWidth::SSE) {
    AlignedArray<const Symbol*, simd::SIMDWidth::SSE, 4> ret;
    ret[3] = lookupSymbol(symbolIter - 1, symbolTable, literals);
    ret[2] = lookupSymbol(symbolIter - 2, symbolTable, literals);
    ret[1] = lookupSymbol(symbolIter - 3, symbolTable, literals);
    ret[0] = lookupSymbol(symbolIter - 4, symbolTable, literals);

    aosToSoa(ArrayView{ret}.template subView<0, 2>(), &unrolledSymbols.frequencies[0], &unrolledSymbols.cumulativeFrequencies[0]);
    aosToSoa(ArrayView{ret}.template subView<2, 2>(), &unrolledSymbols.frequencies[1], &unrolledSymbols.cumulativeFrequencies[1]);
    return {symbolIter - 4, unrolledSymbols};
  } else {
    AlignedArray<const Symbol*, simd::SIMDWidth::SSE, 8> ret;
    ret[7] = lookupSymbol(symbolIter - 1, symbolTable, literals);
    ret[6] = lookupSymbol(symbolIter - 2, symbolTable, literals);
    ret[5] = lookupSymbol(symbolIter - 3, symbolTable, literals);
    ret[4] = lookupSymbol(symbolIter - 4, symbolTable, literals);
    ret[3] = lookupSymbol(symbolIter - 5, symbolTable, literals);
    ret[2] = lookupSymbol(symbolIter - 6, symbolTable, literals);
    ret[1] = lookupSymbol(symbolIter - 7, symbolTable, literals);
    ret[0] = lookupSymbol(symbolIter - 8, symbolTable, literals);

    aosToSoa(ArrayView{ret}.template subView<0, 4>(), &unrolledSymbols.frequencies[0], &unrolledSymbols.cumulativeFrequencies[0]);
    aosToSoa(ArrayView{ret}.template subView<4, 4>(), &unrolledSymbols.frequencies[1], &unrolledSymbols.cumulativeFrequencies[1]);
    return {symbolIter - 8, unrolledSymbols};
  }
};

} // namespace simd
} // namespace internal
} // namespace ranslegacy
} // namespace o2

#endif /* RANSLEGACY_INTERNAL_SIMD_KERNEL_H */