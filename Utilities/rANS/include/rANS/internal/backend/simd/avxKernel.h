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

#ifndef RANS_INTERNAL_SIMD_AVXKERNEL_H
#define RANS_INTERNAL_SIMD_AVXKERNEL_H

//#ifdef __AVX2__

#include <immintrin.h>
#include <cfenv>

#include <tuple>
#include <array>

#include "rANS/internal/backend/simd/utils.h"
#include "rANS/internal/backend/simd/types.h"
#include "rANS/internal/backend/simd/EncoderSymbol.h"

namespace o2
{
namespace rans
{
namespace internal
{
namespace simd
{

namespace detail
{
//  Only works for inputs in the range: [0, 2^52)
inline __m256i _mm256_cvttpd_epu64(__m256d src)
{
  src = _mm256_add_pd(src, _mm256_set1_pd(detail::AlignMantissaMagic));
  return _mm256_xor_si256(
    _mm256_castpd_si256(src),
    _mm256_castpd_si256(_mm256_set1_pd(detail::AlignMantissaMagic)));
}

//  Only works for inputs in the range: [0, 2^52)
inline __m256d _mm256_cvtepu64_pd(__m256i x)
{
  x = _mm256_or_si256(x, _mm256_castpd_si256(_mm256_set1_pd(detail::AlignMantissaMagic)));
  return _mm256_sub_pd(_mm256_castsi256_pd(x), _mm256_set1_pd(detail::AlignMantissaMagic));
}

inline void _mm256_moddiv_pd(__m256d numerator, __m256d denominator, double* div, double* mod)
{
  __m256d divResult = _mm256_floor_pd(_mm256_div_pd(numerator, denominator));
  _mm256_store_pd(mod, _mm256_fnmadd_pd(divResult, denominator, numerator));
  _mm256_store_pd(div, divResult);
}

inline __m256i _mm256_ransencode_pd(__m256i state, __m256d frequency, __m256d cumulative, double normalization)
{
  __m256d div;
  __m256d mod;
  _mm256_moddiv_pd(detail::_mm256_cvtepu64_pd(state), frequency, reinterpret_cast<double*>(&div), reinterpret_cast<double*>(&mod));

  __m256d newState = _mm256_fmadd_pd(_mm256_set1_pd(normalization), div, cumulative);
  newState = _mm256_add_pd(newState, mod);

  return detail::_mm256_cvttpd_epu64(newState);
};

} // namespace detail

//
// uint64 -> double
//
inline pd_t<SIMDWidth::AVX> uint64ToDouble(const epi64_t<SIMDWidth::AVX>& in)
{
  pd_t<SIMDWidth::AVX> out;
  __m256i inReg = _mm256_load_si256(reinterpret_cast<const __m256i*>(in.data()));
  _mm256_store_pd(out.data(), detail::_mm256_cvtepu64_pd(inReg));
  return out;
}

//
// uint32 -> double
//
inline pd_t<SIMDWidth::AVX> int32ToDouble(const epi32_t<SIMDWidth::AVX>& in)
{
  pd_t<SIMDWidth::AVX> out;
  __m128i inReg = _mm_load_si128(reinterpret_cast<const __m128i*>(in.data()));
  _mm256_store_pd(out.data(), _mm256_cvtepi32_pd(inReg));
  return out;
};

//
// double -> uint64
//
inline epi64_t<SIMDWidth::AVX> doubleToUint64(const pd_t<SIMDWidth::AVX>& in)
{
  epi64_t<SIMDWidth::AVX> out;
  __m256d inReg = _mm256_load_pd(in.data());
  _mm256_store_si256(reinterpret_cast<__m256i*>(out.data()), detail::_mm256_cvttpd_epu64(inReg));
  return out;
}

//
// rans Encode
//
inline epi64_t<SIMDWidth::AVX> ransEncode(const epi64_t<SIMDWidth::AVX>& state, const pd_t<SIMDWidth::AVX>& frequency, const pd_t<SIMDWidth::AVX>& cumulative, double normalization)
{
  epi64_t<SIMDWidth::AVX> newState;
  __m256i stateReg = _mm256_load_si256(reinterpret_cast<const __m256i*>(state.data()));
  __m256d frequencyReg = _mm256_load_pd(frequency.data());
  __m256d cumulativeReg = _mm256_load_pd(cumulative.data());
  _mm256_store_si256(reinterpret_cast<__m256i*>(newState.data()), detail::_mm256_ransencode_pd(stateReg, frequencyReg, cumulativeReg, normalization));
  return newState;
};

} // namespace simd
} // namespace internal
} // namespace rans
} // namespace o2

//#endif /* __AVX2__ */
#endif /* RANS_INTERNAL_SIMD_AVXKERNEL_H */