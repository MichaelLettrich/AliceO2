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

#ifndef RANS_INTERNAL_SIMD_AVX512KERNEL_H
#define RANS_INTERNAL_SIMD_AVX512KERNEL_H

#ifdef __AVX512F__

#include <immintrin.h>
#include <cfenv>
#include <cmath>
#include <cassert>

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
inline void _mm512_moddiv_pd(__m512d numerator, __m512d denominator, double* div, double* mod)
{
#pragma STDC FENV_ACCESS ON
  RoundingGuard rounding{RoundingMode::TowardsZero};

  __m512d divResult = _mm512_floor_pd(_mm512_div_pd(numerator, denominator));
  _mm512_store_pd(mod, _mm512_fnmadd_pd(divResult, denominator, numerator));
  _mm512_store_pd(div, divResult);
}

inline __m512i _mm512_ransencode_pd(__m512i state, __m512d frequency, __m512d cumulative, double normalization)
{
  __m512d div;
  __m512d mod;
  _mm512_moddiv_pd(_mm512_cvtepu64_pd(state), frequency, reinterpret_cast<double*>(&div), reinterpret_cast<double*>(&mod));

  __m512d newState = _mm512_fmadd_pd(_mm512_set1_pd(normalization), div, cumulative);
  newState = _mm512_add_pd(newState, mod);

  return _mm512_cvttpd_epu64(newState);
};

} // namespace detail

//
// uint32 -> double
//
inline pd_t<SIMDWidth::AVX512> int32ToDouble(const epi32_t<SIMDWidth::AVX512>& in)
{
  pd_t<SIMDWidth::AVX512> out;
  __m256i inReg = _mm256_load_si256(reinterpret_cast<const __m256i*>(in.data()));
  _mm512_store_pd(out.data(), _mm512_cvtepi32_pd(inReg));
  return out;
};

//
// uint64 -> double
//
inline pd_t<SIMDWidth::AVX512> uint64ToDouble(const epi64_t<SIMDWidth::AVX512>& in)
{
  pd_t<SIMDWidth::AVX512> out;
  __m512i inReg = _mm512_load_si512(reinterpret_cast<const __m512i*>(in.data()));
  _mm512_store_pd(out.data(), _mm512_cvtepu64_pd(inReg));
  return out;
}

//
// double -> uint64
//
inline epi64_t<SIMDWidth::AVX512> doubleToUint64(const pd_t<SIMDWidth::AVX512>& in)
{
  epi64_t<SIMDWidth::AVX512> out;
  __m512d inReg = _mm512_load_pd(in.data());
  _mm512_store_si512(reinterpret_cast<__m512i*>(out.data()), _mm512_cvttpd_epu64(inReg));
  return out;
}

//
// rans Encode
//
inline epi64_t<SIMDWidth::AVX512> ransEncode(const epi64_t<SIMDWidth::AVX512>& state, const pd_t<SIMDWidth::AVX512>& frequency, const pd_t<SIMDWidth::AVX512>& cumulative, double normalization)
{
  epi64_t<SIMDWidth::AVX512> newState;
  __m512d stateReg = _mm512_load_si512(reinterpret_cast<const __m512i*>(state.data()));
  __m512d frequencyReg = _mm512_load_pd(frequency.data());
  __m512d cumulativeReg = _mm512_load_pd(cumulative.data());
  _mm512_store_si512(reinterpret_cast<const __m512i*>(newState.data()), detail::_mm512_ransencode_pd(stateReg, frequencyReg, cumulativeReg, normalization));
  return newState;
};

} // namespace simd
} // namespace internal
} // namespace rans
} // namespace o2

#endif /* __AVX_512F_ */
#endif /* RANS_INTERNAL_SIMD_AVX512KERNEL_H */