// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   simdops.h
/// @author Michael Lettrich
/// @brief wrapper around basic SIMD operations

#ifndef RANS_INTERNAL_COMMON_SIMD_H_
#define RANS_INTERNAL_COMMON_SIMD_H_

#include <immintrin.h>

#include <gsl/span>

#include "rANS/internal/common/utils.h"
#include "rANS/internal/common/simdtypes.h"
#include "rANS/internal/containers/AlignedArray.h"

namespace o2::rans::internal::simd
{

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline __m128i load(gsl::span<const T, getElementCount<T>(SIMDWidth::SSE)> v) noexcept
{
  return _mm_load_si128(reinterpret_cast<const __m128i*>(v.data()));
};

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline __m128i load(gsl::span<T, getElementCount<T>(SIMDWidth::SSE)> v) noexcept
{
  return _mm_load_si128(reinterpret_cast<const __m128i*>(v.data()));
};

inline __m128d load(gsl::span<const double_t, getElementCount<double_t>(SIMDWidth::SSE)> v) noexcept
{
  return _mm_load_pd(v.data());
};

inline __m128d load(gsl::span<double_t, getElementCount<double_t>(SIMDWidth::SSE)> v) noexcept
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
inline __m256i load(gsl::span<const T, getElementCount<T>(SIMDWidth::AVX)> v) noexcept
{
  return _mm256_load_si256(reinterpret_cast<const __m256i*>(v.data()));
};

template <typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline __m256i load(gsl::span<T, getElementCount<T>(SIMDWidth::AVX)> v) noexcept
{
  return _mm256_load_si256(reinterpret_cast<const __m256i*>(v.data()));
};

inline __m256d load(gsl::span<double_t, getElementCount<double_t>(SIMDWidth::AVX)> v) noexcept
{
  return _mm256_load_pd(v.data());
};

inline __m256d load(gsl::span<const double_t, getElementCount<double_t>(SIMDWidth::AVX)> v) noexcept
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
inline void store(__m128i inVec, gsl::span<T, getElementCount<T>(SIMDWidth::SSE)> v) noexcept
{
  _mm_store_si128(reinterpret_cast<__m128i*>(v.data()), inVec);
};

inline void store(__m128d inVec, gsl::span<double_t, getElementCount<double>(SIMDWidth::SSE)> v) noexcept
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
inline void store(__m256i inVec, gsl::span<T, getElementCount<T>(SIMDWidth::AVX)> v) noexcept
{
  _mm256_store_si256(reinterpret_cast<__m256i*>(v.data()), inVec);
};

inline void store(__m256d inVec, gsl::span<double_t, getElementCount<double>(SIMDWidth::AVX)> v) noexcept
{
  _mm256_store_pd(v.data(), inVec);
};

#endif /* __AVX2__ */

template <SIMDWidth width_V>
inline auto setAll(uint64_t value) noexcept
{
  if constexpr (width_V == SIMDWidth::SSE) {
    return _mm_set1_epi64x(value);
  } else {
    return _mm256_set1_epi64x(value);
  }
};

template <SIMDWidth width_V>
inline auto setAll(uint32_t value) noexcept
{
  if constexpr (width_V == SIMDWidth::SSE) {
    return _mm_set1_epi32(value);
  } else {
    return _mm256_set1_epi32(value);
  }
};

template <SIMDWidth width_V>
inline auto setAll(uint16_t value) noexcept
{
  if constexpr (width_V == SIMDWidth::SSE) {
    return _mm_set1_epi16(value);
  } else {
    return _mm256_set1_epi16(value);
  }
};

template <SIMDWidth width_V>
inline auto setAll(uint8_t value) noexcept
{
  if constexpr (width_V == SIMDWidth::SSE) {
    return _mm_set1_epi8(value);
  } else {
    return _mm256_set1_epi8(value);
  }
};

template <SIMDWidth width_V>
inline auto setAll(double_t value) noexcept
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

// adding 2^53 to any IEEE754 double precision floating point number in the range of [0 - 2^52]
// zeros out the exponent and sign bits and the mantissa becomes precisely the integer representation.
inline constexpr double AlignMantissaMagic = 0x0010000000000000; // 2^53

inline __m128d uint64ToDouble(__m128i in) noexcept
{
#if !defined(NDEBUG)
  auto vec = store<uint64_t>(in);
  for (auto i : gsl::make_span(vec)) {
    assert(i < pow2(52));
  }
#endif
  in = _mm_or_si128(in, _mm_castpd_si128(_mm_set1_pd(AlignMantissaMagic)));
  __m128d out = _mm_sub_pd(_mm_castsi128_pd(in), _mm_set1_pd(AlignMantissaMagic));
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
  for (auto i : gsl::make_span(vec)) {
    assert(i < pow2(52));
  }
#endif
  in = _mm256_or_si256(in, _mm256_castpd_si256(_mm256_set1_pd(AlignMantissaMagic)));
  __m256d out = _mm256_sub_pd(_mm256_castsi256_pd(in), _mm256_set1_pd(AlignMantissaMagic));
  return out;
};
#endif /* __AVX2__ */

inline __m128i doubleToUint64(__m128d in) noexcept
{
#if !defined(NDEBUG)
  auto vec = store(in);
  for (auto i : gsl::make_span(vec)) {
    assert(i < pow2(52));
  }
#endif
  in = _mm_add_pd(in, _mm_set1_pd(AlignMantissaMagic));
  __m128i out = _mm_xor_si128(_mm_castpd_si128(in),
                              _mm_castpd_si128(_mm_set1_pd(AlignMantissaMagic)));
  return out;
}

#ifdef __AVX2__

inline __m256i doubleToUint64(__m256d in) noexcept
{
#if !defined(NDEBUG)
  auto vec = store(in);
  for (auto i : gsl::make_span(vec)) {
    assert(i < pow2(52));
  }
#endif

  in = _mm256_add_pd(in, _mm256_set1_pd(AlignMantissaMagic));
  __m256i out = _mm256_xor_si256(_mm256_castpd_si256(in),
                                 _mm256_castpd_si256(_mm256_set1_pd(AlignMantissaMagic)));
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

} // namespace o2::rans::internal::simd

#endif /* RANS_INTERNAL_COMMON_SIMD_H_ */