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

#ifndef RANS_INTERNAL_SIMD_SSEKERNEL_H
#define RANS_INTERNAL_SIMD_SSEKERNEL_H

#ifdef __SSE__

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
inline __m128i _mm_cvttpd_epu64(__m128d src) noexcept
{
  src = _mm_add_pd(src, _mm_set1_pd(AlignMantissaMagic));
  return _mm_xor_si128(
    _mm_castpd_si128(src),
    _mm_castpd_si128(_mm_set1_pd(AlignMantissaMagic)));
}

//  Only works for inputs in the range: [0, 2^52)
inline __m128d _mm_cvtepu64_pd(__m128i x) noexcept
{
  x = _mm_or_si128(x, _mm_castpd_si128(_mm_set1_pd(AlignMantissaMagic)));
  return _mm_sub_pd(_mm_castsi128_pd(x), _mm_set1_pd(AlignMantissaMagic));
}

inline void _mm_moddiv_pd(__m128d numerator, __m128d denominator, double* div, double* mod) noexcept
{
  __m128d divResult = _mm_floor_pd(_mm_div_pd(numerator, denominator));
  _mm_store_pd(mod, _mm_fnmadd_pd(divResult, denominator, numerator));
  _mm_store_pd(div, divResult);
}

inline __m128i _mm_ransencode_pd(__m128i state, __m128d frequency, __m128d cumulative, double normalization) noexcept
{
  __m128d div;
  __m128d mod;
  _mm_moddiv_pd(detail::_mm_cvtepu64_pd(state), frequency, reinterpret_cast<double*>(&div), reinterpret_cast<double*>(&mod));

  __m128d newState = _mm_fmadd_pd(_mm_set1_pd(normalization), div, cumulative);
  newState = _mm_add_pd(newState, mod);

  return detail::_mm_cvttpd_epu64(newState);
};
} // namespace detail

//
// uint32 -> double
//
inline pd_t<SIMDWidth::SSE> int32ToDouble(const epi32_t<SIMDWidth::SSE>& in) noexcept
{
  pd_t<SIMDWidth::SSE> out;
  __m128i inReg = _mm_load_si128(reinterpret_cast<const __m128i*>(in.data()));
  _mm_store_pd(out.data(), _mm_cvtepi32_pd(inReg));
  return out;
};

//
// uint64 -> double
// Only works for inputs in the range: [0, 2^52)
//
inline pd_t<SIMDWidth::SSE> uint64ToDouble(const epi64_t<SIMDWidth::SSE>& in) noexcept
{
  pd_t<SIMDWidth::SSE> out;
  __m128i inReg = _mm_load_si128(reinterpret_cast<const __m128i*>(in.data()));
  _mm_store_pd(out.data(), detail::_mm_cvtepu64_pd(inReg));
  return out;
}

//
// double -> uint64
//
inline epi64_t<SIMDWidth::SSE> doubleToUint64(const pd_t<SIMDWidth::SSE>& in) noexcept
{
  epi64_t<SIMDWidth::SSE> out;
  __m128d inReg = _mm_load_pd(in.data());
  _mm_store_si128(reinterpret_cast<__m128i*>(out.data()), detail::_mm_cvttpd_epu64(inReg));
  return out;
}

//
// rans Encode
//
inline epi64_t<SIMDWidth::SSE> ransEncode(const epi64_t<SIMDWidth::SSE>& state, const pd_t<SIMDWidth::SSE>& frequency, const pd_t<SIMDWidth::SSE>& cumulative, uint32_t normalization) noexcept
{
  epi64_t<SIMDWidth::SSE> newState;
  __m128i stateReg = _mm_load_si128(reinterpret_cast<const __m128i*>(state.data()));
  __m128d frequencyReg = _mm_load_pd(frequency.data());
  __m128d cumulativeReg = _mm_load_pd(cumulative.data());
  _mm_store_si128(reinterpret_cast<__m128i*>(newState.data()), detail::_mm_ransencode_pd(stateReg, frequencyReg, cumulativeReg, normalization));
  return newState;
};

} // namespace simd
} // namespace internal
} // namespace rans
} // namespace o2

#endif /* __SSE__ */
#endif /* RANS_INTERNAL_SIMD_SSEKERNEL_H */