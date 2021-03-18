// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   tyspes.h
/// @author Michael Lettrich
/// @since  2021-06-12
/// @brief

#ifndef RANS_INTERNAL_SIMD_TYPES_H
#define RANS_INTERNAL_SIMD_TYPES_H

#include <cstring>
#include <cassert>

#include "rANS/internal/backend/simd/utils.h"

namespace o2
{
namespace rans
{
namespace internal
{
namespace simd
{

inline constexpr size_t simdAlign(size_t SIMDwidth) noexcept
{
  return SIMDwidth * 8;
}

enum SIMDWidth : unsigned int { SSE = 2,
                                AVX = 4,
                                AVX512 = 8 };

enum Alignment : unsigned int { SSEaligned = 16,
                                AVXaligned = 32,
                                AVX512aligned = 64 };

template <typename T, size_t SIMDWidth_V>
class alignas(simdAlign(SIMDWidth_V)) AlignedArray
{
 public:
  using value_type = T;
  using pointer = T*;
  using reference = T&;
  using iterator = T*;
  using const_iterator = const T*;

  inline constexpr AlignedArray() noexcept {}; //NOLINT

  inline constexpr AlignedArray(T elem) noexcept
  {
#pragma omp simd
    for (size_t i = 0; i < SIMDWidth_V; ++i) {
      mData[i] = elem;
    }
  };

  template <typename... Args, std::enable_if_t<(sizeof...(Args) == SIMDWidth_V), bool> = true>
  inline constexpr AlignedArray(Args... args) noexcept : mData{args...} {};

  inline constexpr size_t size() const noexcept { return SIMDWidth_V; };
  inline constexpr const T* data() const noexcept { return mData; };
  inline constexpr T* data() noexcept { return const_cast<T*>(static_cast<const AlignedArray&>(*this).data()); };
  inline constexpr const_iterator begin() const noexcept { return data(); };
  inline constexpr const_iterator end() const noexcept { return data() + size(); };
  inline constexpr iterator begin() noexcept { return const_cast<iterator>(static_cast<const AlignedArray&>(*this).begin()); };
  inline constexpr iterator end() noexcept { return const_cast<iterator>(static_cast<const AlignedArray&>(*this).end()); };
  inline constexpr const T& operator[](size_t i) const
  {
    assert(i < SIMDWidth_V);
    return mData[i];
  };
  inline constexpr T& operator[](size_t i) { return const_cast<T&>(static_cast<const AlignedArray&>(*this)[i]); };

 private:
  T mData[SIMDWidth_V]{};
};

template <typename T>
struct array_size : std::integral_constant<size_t, 0> {
};
template <typename T, size_t N>
struct array_size<std::array<T, N>> : std::integral_constant<size_t, N> {
};
template <typename T, size_t N>
struct array_size<AlignedArray<T, N>> : std::integral_constant<size_t, N> {
};

template <typename T>
using getSIMDWidth = array_size<T>;

template <size_t simdWidth_V>
using simdpd_t = AlignedArray<double, simdWidth_V>;
template <size_t simdWidth_V>
using simdepi64_t = AlignedArray<uint64_t, simdWidth_V>;
template <size_t simdWidth_V>
using simdepi32_t = AlignedArray<uint32_t, simdWidth_V>;

} // namespace simd
} // namespace internal
} // namespace rans
} // namespace o2

#endif /* RANS_INTERNAL_SIMD_TYPES_H */