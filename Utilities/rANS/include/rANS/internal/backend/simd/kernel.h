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

#ifndef RANS_INTERNAL_SIMD_KERNEL_H
#define RANS_INTERNAL_SIMD_KERNEL_H

#include <immintrin.h>
#include <cfenv>

#include <array>
#include <tuple>

#include "rANS/internal/backend/simd/Symbol.h"
#include "rANS/internal/backend/simd/types.h"
#include "rANS/internal/backend/simd/utils.h"
#include "rANS/internal/helper.h"

namespace o2
{
namespace rans
{
namespace internal
{
namespace simd
{
template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline __m128i load(const AlignedArray<T, SIMDWidth::SSE>& v) noexcept
{
  return _mm_load_si128(reinterpret_cast<const __m128i*>(v.data()));
};

inline __m128d load(const AlignedArray<double, SIMDWidth::SSE>& v) noexcept
{
  return _mm_load_pd(v.data());
};

#ifdef __AVX2__

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline __m256i load(const AlignedArray<T, SIMDWidth::AVX>& v) noexcept
{
  return _mm256_load_si256(reinterpret_cast<const __m256i*>(v.data()));
};

inline __m256d load(const AlignedArray<double, SIMDWidth::AVX>& v) noexcept
{
  return _mm256_load_pd(v.data());
};

#endif /* __AVX2__ */
#ifdef __AVX512F__

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline __m512i load(const AlignedArray<T, SIMDWidth::AVX512>& v) noexcept
{
  return _mm512_load_si512(reinterpret_cast<const __m512i*>(v.data()));
};

inline __m512d load(const AlignedArray<double, SIMDWidth::AVX512>& v) noexcept
{
  return _mm512_load_pd(v.data());
};

#endif /* __AVX_512F_ */

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline AlignedArray<T, SIMDWidth::SSE> store(__m128i inVec) noexcept
{
  AlignedArray<T, SIMDWidth::SSE> out;
  _mm_store_si128(reinterpret_cast<__m128i*>(out.data()), inVec);
  return out;
};

inline AlignedArray<double, SIMDWidth::SSE> store(__m128d inVec) noexcept
{
  AlignedArray<double, SIMDWidth::SSE> out;
  _mm_store_pd(out.data(), inVec);
  return out;
};

#ifdef __AVX2__

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline AlignedArray<T, SIMDWidth::AVX> store(__m256i inVec) noexcept
{
  AlignedArray<T, SIMDWidth::AVX> out;
  _mm256_store_si256(reinterpret_cast<__m256i*>(out.data()), inVec);
  return out;
};

inline AlignedArray<double, SIMDWidth::AVX> store(__m256d inVec) noexcept
{
  AlignedArray<double, SIMDWidth::AVX> out;
  _mm256_store_pd(out.data(), inVec);
  return out;
};

#endif /* __AVX2__ */
#ifdef __AVX512F__

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline AlignedArray<T, SIMDWidth::AVX512> store(__m512i inVec) noexcept
{
  AlignedArray<T, SIMDWidth::AVX512> out;
  _mm512_store_si512(reinterpret_cast<__m512i*>(out.data()), inVec);
  return out;
};

inline AlignedArray<double, SIMDWidth::AVX512> store(__m512d inVec) noexcept
{
  AlignedArray<double, SIMDWidth::AVX512> out;
  _mm512_store_pd(out.data(), inVec);
  return out;
};

#endif /* __AVX_512F_ */

//
// uint32 -> double
//
template <SIMDWidth width_V>
inline auto int32ToDouble(const epi32_t<SIMDWidth::SSE>& in) noexcept
{
  auto inReg = load(in);
  if constexpr (width_V == SIMDWidth::SSE) {
    return store(_mm_cvtepi32_pd(inReg));
  } else if constexpr (width_V == SIMDWidth::AVX) {
#ifdef __AVX2__
    return store(_mm256_cvtepi32_pd(inReg));
#endif /* __AVX2__ */
  }
};

template <SIMDWidth width_V>
inline auto int32ToDouble(const epi32_t<SIMDWidth::AVX>& in) noexcept
{
  if constexpr (width_V == SIMDWidth::AVX) {
#ifdef __AVX2__
    __m128i inReg = _mm_load_si128(reinterpret_cast<const __m128i*>(in.data()));
    return store(_mm256_cvtepi32_pd(inReg));
#endif /* __AVX2__ */
  } else if constexpr (width_V == SIMDWidth::AVX512) {
#ifdef __AVX512F__
    auto inReg = load(in);
    return store(_mm512_cvtepi32_pd(inReg));
#endif /* __AVX_512F_ */
  }
};

//
// uint64 -> double
// Only works for inputs in the range: [0, 2^52)
//
inline pd_t<SIMDWidth::SSE> uint64ToDouble(const epi64_t<SIMDWidth::SSE>& in) noexcept
{
#if !defined(NDEBUG)
  for (auto i : in) {
    assert(i < pow2(52));
  }
#endif

  __m128i inReg = load(in);
  inReg = _mm_or_si128(inReg, _mm_castpd_si128(_mm_set1_pd(detail::AlignMantissaMagic)));
  __m128d outReg = _mm_sub_pd(_mm_castsi128_pd(inReg), _mm_set1_pd(detail::AlignMantissaMagic));
  return store(outReg);
};

#ifdef __AVX2__

//
// uint64 -> double
// Only works for inputs in the range: [0, 2^52)
//
inline pd_t<SIMDWidth::AVX> uint64ToDouble(const epi64_t<SIMDWidth::AVX>& in) noexcept
{
#if !defined(NDEBUG)
  for (auto i : in) {
    assert(i < pow2(52));
  }
#endif
  __m256i inReg = load(in);
  inReg = _mm256_or_si256(inReg, _mm256_castpd_si256(_mm256_set1_pd(detail::AlignMantissaMagic)));
  __m256d outReg = _mm256_sub_pd(_mm256_castsi256_pd(inReg), _mm256_set1_pd(detail::AlignMantissaMagic));
  return store(outReg);
};

#endif /* __AVX2__ */
#ifdef __AVX512F__

//
// uint64 -> double
//
inline pd_t<SIMDWidth::AVX512> uint64ToDouble(const epi64_t<SIMDWidth::AVX512>& in)
{
  auto inReg = load(in);
  auto outReg = _mm512_cvtepu64_pd(inReg);
  return store(outReg);
}
#endif /* __AVX_512F_ */

//
// double -> uint64
// Only works for inputs in the range: [0, 2^52)
//
inline epi64_t<SIMDWidth::SSE> doubleToUint64(const pd_t<SIMDWidth::SSE>& in) noexcept
{
#if !defined(NDEBUG)
  for (auto i : in) {
    assert(i < pow2(52));
  }
#endif
  __m128d inReg = load(in);
  inReg = _mm_add_pd(inReg, _mm_set1_pd(detail::AlignMantissaMagic));
  __m128i outReg = _mm_xor_si128(_mm_castpd_si128(inReg),
                                 _mm_castpd_si128(_mm_set1_pd(detail::AlignMantissaMagic)));
  return store<uint64_t>(outReg);
}

#ifdef __AVX2__

//
// double -> uint64
// Only works for inputs in the range: [0, 2^52)
//
inline epi64_t<SIMDWidth::AVX> doubleToUint64(const pd_t<SIMDWidth::AVX>& in) noexcept
{
#if !defined(NDEBUG)
  for (auto i : in) {
    assert(i < pow2(52));
  }
#endif

  __m256d inReg = load(in);
  inReg = _mm256_add_pd(inReg, _mm256_set1_pd(detail::AlignMantissaMagic));
  __m256i outReg = _mm256_xor_si256(_mm256_castpd_si256(inReg),
                                    _mm256_castpd_si256(_mm256_set1_pd(detail::AlignMantissaMagic)));
  return store<uint64_t>(outReg);
}

#endif /* __AVX2__ */
#ifdef __AVX512F__

//
// double -> uint64
//
inline epi64_t<SIMDWidth::AVX512> doubleToUint64(const pd_t<SIMDWidth::AVX512>& in)
{
  auto inReg = load(in);
  auto outReg = _mm512_cvttpd_epu64(inReg);
  return store<uint64_t>(outReg);
}

#endif /* __AVX_512F_ */

// calculate both floor(a/b) and a%b
inline std::tuple<pd_t<SIMDWidth::SSE>, pd_t<SIMDWidth::SSE>> divMod(const pd_t<SIMDWidth::SSE>& numerator, const pd_t<SIMDWidth::SSE>& denominator) noexcept
{
  __m128d numeratorReg = load(numerator);
  __m128d denominatorReg = load(denominator);
  __m128d div = _mm_floor_pd(_mm_div_pd(numeratorReg, denominatorReg));
  __m128d mod = _mm_fnmadd_pd(div, denominatorReg, numeratorReg);
  return std::make_tuple(store(div), store(mod));
}

#ifdef __AVX2__

// calculate both floor(a/b) and a%b
inline std::tuple<pd_t<SIMDWidth::AVX>, pd_t<SIMDWidth::AVX>> divMod(const pd_t<SIMDWidth::AVX>& numerator, const pd_t<SIMDWidth::AVX>& denominator) noexcept
{
  __m256d numeratorReg = load(numerator);
  __m256d denominatorReg = load(denominator);
  __m256d div = _mm256_floor_pd(_mm256_div_pd(numeratorReg, denominatorReg));
  __m256d mod = _mm256_fnmadd_pd(div, denominatorReg, numeratorReg);
  return std::make_tuple(store(div), store(mod));
}

#endif /* __AVX2__ */
#ifdef __AVX512F__

// calculate both floor(a/b) and a%b
inline std::tuple<pd_t<SIMDWidth::AVX512>, pd_t<SIMDWidth::AVX512>> divMod(const pd_t<SIMDWidth::AVX512>& numerator, const pd_t<SIMDWidth::AVX512>& denominator) noexcept
{
  __m512d numeratorReg = load(numerator);
  __m512d denominatorReg = load(denominator);
  __m512d div = _mm512_floor_pd(_mm512_div_pd(numeratorReg, denominatorReg));
  __m512d mod = _mm512_fnmadd_pd(div, denominatorReg, numeratorReg);
  return std::make_tuple(store(div), store(mod));
}

#endif /* __AVX_512F_ */

//
// rans Encode
//
inline epi64_t<SIMDWidth::SSE> ransEncode(const epi64_t<SIMDWidth::SSE>& state, const pd_t<SIMDWidth::SSE>& frequency, const pd_t<SIMDWidth::SSE>& cumulative, const pd_t<SIMDWidth::SSE>& normalization) noexcept
{
#if !defined(NDEBUG)
  for (auto i : state) {
    assert(i < pow2(52));
  }
#endif

  auto [div, mod] = divMod(uint64ToDouble(state), frequency);
  auto divReg = load(div);
  auto modReg = load(mod);
  auto cumulativeReg = load(cumulative);
  auto normalizationReg = load(normalization);
  auto newState = _mm_fmadd_pd(normalizationReg, divReg, cumulativeReg);
  newState = _mm_add_pd(newState, modReg);

  return doubleToUint64(store(newState));
};

#ifdef __AVX2__

//
// rans Encode
//
inline epi64_t<SIMDWidth::AVX> ransEncode(const epi64_t<SIMDWidth::AVX>& state, const pd_t<SIMDWidth::AVX>& frequency, const pd_t<SIMDWidth::AVX>& cumulative, const pd_t<SIMDWidth::AVX>& normalization) noexcept
{
#if !defined(NDEBUG)
  for (auto i : state) {
    assert(i < pow2(52));
  }
#endif

  auto [div, mod] = divMod(uint64ToDouble(state), frequency);
  auto divReg = load(div);
  auto modReg = load(mod);
  auto cumulativeReg = load(cumulative);
  auto normalizationReg = load(normalization);
  auto newState = _mm256_fmadd_pd(normalizationReg, divReg, cumulativeReg);
  newState = _mm256_add_pd(newState, modReg);

  return doubleToUint64(store(newState));
};

#endif /* __AVX2__ */
#ifdef __AVX512F__

//
// rans Encode
//
inline epi64_t<SIMDWidth::AVX512> ransEncode(const epi64_t<SIMDWidth::AVX512>& state, const pd_t<SIMDWidth::AVX512>& frequency, const pd_t<SIMDWidth::AVX512>& cumulative, const pd_t<SIMDWidth::AVX512>& normalization) noexcept
{
  auto [div, mod] = divMod(uint64ToDouble(state), frequency);
  auto divReg = load(div);
  auto modReg = load(mod);
  auto cumulativeReg = load(cumulative);
  auto normalizationReg = load(normalization);
  auto newState = _mm512_fmadd_pd(normalizationReg, divReg, cumulativeReg);
  newState = _mm512_add_pd(newState, modReg);

  return doubleToUint64(store(newState));
};

#endif /* __AVX_512F_ */

inline std::tuple<simd::epi32_t<simd::SIMDWidth::SSE>, simd::epi32_t<simd::SIMDWidth::SSE>> aosToSoa(const std::array<const Symbol*, 2>& in) noexcept
{
  __m128i in0Reg = _mm_load_si128(reinterpret_cast<__m128i const*>(in[0]->data()));
  __m128i in1Reg = _mm_load_si128(reinterpret_cast<__m128i const*>(in[1]->data()));

  __m128i res0Reg = _mm_unpacklo_epi32(in0Reg, in1Reg);
  __m128i res1Reg = _mm_shuffle_epi32(res0Reg, _MM_SHUFFLE(0, 0, 3, 2));

  return {store<uint32_t>(res0Reg), store<uint32_t>(res1Reg)};
};

#ifdef __AVX2__

inline std::tuple<epi32_t<SIMDWidth::SSE>, epi32_t<SIMDWidth::SSE>> aosToSoa(const std::array<const Symbol*, 4>& in) noexcept
{
  __m128i in0Reg = _mm_load_si128(reinterpret_cast<__m128i const*>(in[0]->data()));
  __m128i in1Reg = _mm_load_si128(reinterpret_cast<__m128i const*>(in[1]->data()));
  __m128i in2Reg = _mm_load_si128(reinterpret_cast<__m128i const*>(in[2]->data()));
  __m128i in3Reg = _mm_load_si128(reinterpret_cast<__m128i const*>(in[3]->data()));

  __m128i merged0Reg = _mm_unpacklo_epi32(in0Reg, in1Reg);
  __m128i merged1Reg = _mm_unpacklo_epi32(in2Reg, in3Reg);
  __m128i res0Reg = _mm_unpacklo_epi64(merged0Reg, merged1Reg);
  __m128i res1Reg = _mm_unpackhi_epi64(merged0Reg, merged1Reg);

  return {store<uint32_t>(res0Reg), store<uint32_t>(res1Reg)};
};

inline std::tuple<epi32_t<SIMDWidth::AVX>, epi32_t<SIMDWidth::AVX>> aosToSoa(const std::array<const Symbol*, 8>& in) noexcept
{
  std::array<const Symbol*, 4> first{in[0], in[1], in[2], in[3]};
  std::array<const Symbol*, 4> second{in[4], in[5], in[6], in[7]};

  auto [arr00, arr01] = aosToSoa(first);
  auto [arr10, arr11] = aosToSoa(second);

  auto arr00Reg = load(arr00);
  auto arr01Reg = load(arr01);
  auto arr10Reg = load(arr10);
  auto arr11Reg = load(arr11);

  auto res0Reg = _mm256_set_m128i(arr10Reg, arr00Reg);
  auto res1Reg = _mm256_set_m128i(arr11Reg, arr01Reg);

  return {store<uint32_t>(res0Reg), store<uint32_t>(res1Reg)};
};
#endif /* __AVX2__ */

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

template <SIMDWidth width_V>
inline epi64_t<width_V> cmpgeq_epi64(epi64_t<width_V> a, epi64_t<width_V> b) noexcept
{
  auto aVec = load(a);
  auto bVec = load(b);
  return store<uint64_t>(cmpgeq_epi64(aVec, bVec));
};

template <SIMDWidth width_V, uint64_t lowerBound_V, uint8_t streamBits_V>
inline auto computeMaxState(__m128i frequencyVec, uint8_t symbolTablePrecisionBits) noexcept
{
  const uint64_t xmax = (lowerBound_V >> symbolTablePrecisionBits) << streamBits_V;
  const uint8_t shift = log2UInt(xmax);
  if constexpr (width_V == SIMDWidth::SSE) {
    __m128i frequencyVecEpi64 = _mm_cvtepi32_epi64(frequencyVec);
    return _mm_slli_epi64(frequencyVecEpi64, shift);
  }
  if constexpr (width_V == SIMDWidth::AVX) {
    __m256i frequencyVecEpi64 = _mm256_cvtepi32_epi64(frequencyVec);
    return _mm256_slli_epi64(frequencyVecEpi64, shift);
  }
};

template <uint8_t streamBits_V>
inline epi64_t<SIMDWidth::SSE> computeNewState(__m128i stateVec, __m128i cmpVec) noexcept
{
  //newState = (state >= maxState) ? state >> streamBits_V : state
  __m128i newStateVec = _mm_srli_epi64(stateVec, streamBits_V);
  newStateVec = _mm_blendv_epi8(stateVec, newStateVec, cmpVec);
  return store<uint64_t>(newStateVec);
};

#ifdef __AVX2__
template <uint8_t streamBits_V>
inline epi64_t<SIMDWidth::AVX> computeNewState(__m256i stateVec, __m256i cmpVec) noexcept
{
  //newState = (state >= maxState) ? state >> streamBits_V : state
  __m256i newStateVec = _mm256_srli_epi64(stateVec, streamBits_V);
  newStateVec = _mm256_blendv_epi8(stateVec, newStateVec, cmpVec);
  return store<uint64_t>(newStateVec);
};
#endif /* __AVX2__ */

inline constexpr std::array<epi8_t<SIMDWidth::SSE>, 6>
  SSEPermutationLUT{{
    {0xFF_u8},                                                                                                                                        //0b0000 valid
    {0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0001 valid
    {0xFF_u8},                                                                                                                                        //0b0010 invalid
    {0xFF_u8},                                                                                                                                        //0b0011 invalid
    {0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0100 valid
    {0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}  //0b0101  valid
  }};

inline constexpr std::array<epi32_t<SIMDWidth::AVX>, 16>
  AVXPermutationLUT{{
    {0x0u, 0x1u, 0x2u, 0x3u, 0x0u, 0x1u, 0x2u, 0x3u}, //0b0000
    {0x0u, 0x1u, 0x2u, 0x3u, 0x0u, 0x1u, 0x2u, 0x3u}, //0b0001
    {0x1u, 0x0u, 0x2u, 0x3u, 0x0u, 0x1u, 0x2u, 0x3u}, //0b0010
    {0x0u, 0x1u, 0x2u, 0x3u, 0x0u, 0x1u, 0x2u, 0x3u}, //0b0011
    {0x2u, 0x1u, 0x0u, 0x3u, 0x0u, 0x1u, 0x2u, 0x3u}, //0b0100
    {0x0u, 0x2u, 0x1u, 0x3u, 0x0u, 0x1u, 0x2u, 0x3u}, //0b0101
    {0x1u, 0x2u, 0x0u, 0x3u, 0x0u, 0x1u, 0x2u, 0x3u}, //0b0110
    {0x0u, 0x1u, 0x2u, 0x3u, 0x0u, 0x1u, 0x2u, 0x3u}, //0b0111
    {0x3u, 0x0u, 0x0u, 0x0u, 0x0u, 0x1u, 0x2u, 0x3u}, //0b1000
    {0x0u, 0x3u, 0x1u, 0x2u, 0x0u, 0x1u, 0x2u, 0x3u}, //0b1001
    {0x1u, 0x3u, 0x0u, 0x2u, 0x0u, 0x1u, 0x2u, 0x3u}, //0b1010
    {0x0u, 0x1u, 0x3u, 0x2u, 0x0u, 0x1u, 0x2u, 0x3u}, //0b1011
    {0x2u, 0x3u, 0x0u, 0x1u, 0x0u, 0x1u, 0x2u, 0x3u}, //0b1100
    {0x0u, 0x2u, 0x3u, 0x1u, 0x0u, 0x1u, 0x2u, 0x3u}, //0b1101
    {0x1u, 0x2u, 0x3u, 0x0u, 0x0u, 0x1u, 0x2u, 0x3u}, //0b1110
    {0x0u, 0x1u, 0x2u, 0x3u, 0x0u, 0x1u, 0x2u, 0x3u}  //0b1111
  }};

inline std::tuple<uint32_t, epi32_t<SIMDWidth::SSE>> streamOut(__m128i stateVec, __m128i cmpVec) noexcept
{
  constexpr epi32_t<SIMDWidth::SSE> extractionMask{0xFFFFFFFFu, 0x0u, 0xFFFFFFFFu, 0x0u};
  __m128i extractionMaskVec = load(extractionMask);
  __m128i streamOutMaskVec = _mm_and_si128(cmpVec, extractionMaskVec);
  const uint32_t id = _mm_movemask_ps(_mm_castsi128_ps(streamOutMaskVec));

  __m128i shuffleMaskVec = load(SSEPermutationLUT[id]);
  __m128i streamOutVec = _mm_shuffle_epi8(stateVec, shuffleMaskVec);

  return {_mm_popcnt_u32(id), store<uint32_t>(streamOutVec)};
};

#ifdef __AVX2__
inline std::tuple<uint32_t, epi32_t<SIMDWidth::AVX>> streamOut(__m256i stateVec, __m256i cmpVec) noexcept
{
  constexpr epi32_t<SIMDWidth::AVX> extractionMask{0xFFFFFFFFu, 0x0u, 0xFFFFFFFFu, 0x0u, 0xFFFFFFFFu, 0x0u, 0xFFFFFFFFu, 0x0u};
  __m256i extractionMaskVec = load(extractionMask);
  __m256i streamOutMaskVec = _mm256_and_si256(cmpVec, extractionMaskVec);

  // take [0,2,4,6]th element from uint32_t streamOutVector [0,1,2,3,4,5,6,7]
  constexpr epi32_t<SIMDWidth::AVX> permutationMask{0x6u, 0x4u, 0x2u, 0x0u, 0x7u, 0x7u, 0x7u, 0x7u};
  __m256i streamOutVec = _mm256_permutevar8x32_epi32(stateVec, load(permutationMask));
  streamOutMaskVec = _mm256_permutevar8x32_epi32(streamOutMaskVec, load(permutationMask));
  // condense streamOut Vector mask from SIMD vector into an integer ID
  const uint32_t id = _mm256_movemask_ps(_mm256_castsi256_ps(streamOutMaskVec));
  // use ID to fetch right reordering from permutation LUT and move data into the right order
  streamOutVec = _mm256_castps_si256(_mm256_permutevar_ps(_mm256_castsi256_ps(streamOutVec), load(AVXPermutationLUT[id])));

  return {_mm_popcnt_u32(id), store<uint32_t>(streamOutVec)};
};
#endif /* __AVX2__ */

namespace impl
{
template <SIMDWidth width_V, typename output_IT, uint64_t lowerBound_V, uint8_t streamBits_V>
inline auto ransRenorm(const epi64_t<width_V>& states, const epi32_t<SIMDWidth::SSE>& frequencies, uint8_t symbolTablePrecisionBits, output_IT outputIter) noexcept -> std::tuple<output_IT, epi64_t<width_V>>
{
  auto stateVec = load(states);
  __m128i frequencyVec = load(frequencies);

  // calculate maximum state
  auto maxStateVec = computeMaxState<width_V, lowerBound_V, streamBits_V>(frequencyVec, symbolTablePrecisionBits);
  //cmpVec = (state >= maxState)
  auto cmpVec = cmpgeq_epi64(stateVec, maxStateVec);
  //newState = (state >= maxState) ? state >> streamBits_V : state
  epi64_t<width_V> newState = computeNewState<streamBits_V>(stateVec, cmpVec);

  const auto [nStreamOutWords, streamOutResult] = streamOut(stateVec, cmpVec);

  for (size_t i = 0; i < nStreamOutWords; ++i) {
    *(++outputIter) = streamOutResult[i];
  }

  return std::make_tuple(outputIter, newState);
};

} // namespace impl

template <typename output_IT, uint64_t lowerBound_V, uint8_t streamBits_V>
inline auto ransRenorm(const epi64_t<SIMDWidth::SSE>& states, const epi32_t<SIMDWidth::SSE>& frequencies, uint8_t symbolTablePrecisionBits, output_IT outputIter) noexcept -> std::tuple<output_IT, epi64_t<SIMDWidth::SSE>>
{
  return impl::ransRenorm<SIMDWidth::SSE, output_IT, lowerBound_V, streamBits_V>(states, frequencies, symbolTablePrecisionBits, outputIter);
};

#ifdef __AVX2__
template <typename output_IT, uint64_t lowerBound_V, uint8_t streamBits_V>
inline auto ransRenorm(const epi64_t<SIMDWidth::AVX>& states, const epi32_t<SIMDWidth::SSE>& frequencies, uint8_t symbolTablePrecisionBits, output_IT outputIter) noexcept -> std::tuple<output_IT, epi64_t<SIMDWidth::AVX>>
{
  return impl::ransRenorm<SIMDWidth::AVX, output_IT, lowerBound_V, streamBits_V>(states, frequencies, symbolTablePrecisionBits, outputIter);
};
#endif /* __AVX2__ */

} // namespace simd
} // namespace internal
} // namespace rans
} // namespace o2

#endif /* RANS_INTERNAL_SIMD_KERNEL_H */