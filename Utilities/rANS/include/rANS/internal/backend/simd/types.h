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
#include <cstdint>
#include <type_traits>

#include <sstream>
#include <fmt/format.h>

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

enum class SIMDWidth : uint32_t { SSE = 128,
                                  AVX = 256,
                                  AVX512 = 512 };

inline constexpr size_t getLaneWidthBits(SIMDWidth width) noexcept { return static_cast<size_t>(width); };

inline constexpr size_t getLaneWidthBytes(SIMDWidth width) noexcept { return toBytes(static_cast<size_t>(width)); };

inline constexpr size_t getAlignment(SIMDWidth width) noexcept { return getLaneWidthBytes(width); };

template <typename T>
inline constexpr size_t getElementCount(SIMDWidth width) noexcept
{
  return getLaneWidthBytes(width) / sizeof(T);
};

template <typename T, SIMDWidth simdWidth_V>
class alignas(getAlignment(simdWidth_V)) AlignedArray
{
 public:
  using value_type = T;
  using pointer = T*;
  using reference = T&;
  using iterator = T*;
  using const_iterator = const T*;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  inline constexpr AlignedArray() noexcept {}; //NOLINT

  inline constexpr AlignedArray(T elem) noexcept
  {
    //#pragma omp simd
    for (size_t i = 0; i < size(); ++i) {
      mData[i] = elem;
    }
  };

  template <typename... Args, std::enable_if_t<(sizeof...(Args) == getElementCount<T>(simdWidth_V)), bool> = true>
  inline constexpr AlignedArray(Args... args) noexcept : mData{args...} {};

  inline constexpr size_t size() const noexcept { return getElementCount<T>(simdWidth_V); };
  inline constexpr const T* data() const noexcept { return mData; };
  inline constexpr T* data() noexcept { return const_cast<T*>(static_cast<const AlignedArray&>(*this).data()); };
  inline constexpr const_iterator begin() const noexcept { return data(); };
  inline constexpr const_iterator end() const noexcept { return data() + size(); };
  inline constexpr iterator begin() noexcept { return const_cast<iterator>(static_cast<const AlignedArray&>(*this).begin()); };
  inline constexpr iterator end() noexcept { return const_cast<iterator>(static_cast<const AlignedArray&>(*this).end()); };
  inline constexpr const_reverse_iterator rbegin() const noexcept { return std::reverse_iterator(this->end()); };
  inline constexpr const_reverse_iterator rend() const noexcept { return std::reverse_iterator(this->begin()); };
  inline constexpr reverse_iterator rbegin() noexcept { return std::reverse_iterator(this->end()); };
  inline constexpr reverse_iterator rend() noexcept { return std::reverse_iterator(this->begin()); };
  inline constexpr const T& operator[](size_t i) const
  {
    assert(i < size());
    return mData[i];
  };
  inline constexpr T& operator[](size_t i) { return const_cast<T&>(static_cast<const AlignedArray&>(*this)[i]); };

 private:
  T mData[getElementCount<T>(simdWidth_V)]{};
};

template <typename T>
struct simdWidth;

template <typename T, SIMDWidth simd_V>
struct simdWidth<AlignedArray<T, simd_V>> : public std::integral_constant<SIMDWidth, simd_V> {
};

template <typename T>
inline constexpr SIMDWidth getSimdWidth_v = simdWidth<T>::value;

template <typename T>
struct elementCount;

template <typename T, SIMDWidth simd_V>
struct elementCount<AlignedArray<T, simd_V>> : public std::integral_constant<size_t, getElementCount<T>(simd_V)> {
};

template <typename T>
inline constexpr size_t elementCount_v = elementCount<T>::value;

template <SIMDWidth simd_V>
using pd_t = AlignedArray<double, simd_V>;
template <SIMDWidth simd_V>
using epi64_t = AlignedArray<uint64_t, simd_V>;
template <SIMDWidth simd_V>
using epi32_t = AlignedArray<uint32_t, simd_V>;
template <SIMDWidth simd_V>
using epi16_t = AlignedArray<uint16_t, simd_V>;
template <SIMDWidth simd_V>
using epi8_t = AlignedArray<uint8_t, simd_V>;

inline constexpr std::uint8_t operator"" _u8(unsigned long long int value) { return static_cast<uint8_t>(value); };
inline constexpr std::int8_t operator"" _i8(unsigned long long int value) { return static_cast<int8_t>(value); };

inline constexpr std::uint16_t operator"" _u16(unsigned long long int value) { return static_cast<uint16_t>(value); };
inline constexpr std::int16_t operator"" _i16(unsigned long long int value) { return static_cast<int16_t>(value); };

template <typename T, SIMDWidth width_V>
std::ostream& operator<<(std::ostream& stream, const AlignedArray<T, width_V>& vec)
{
  stream << "[";
  for (const auto& elem : vec) {
    stream << elem << ", ";
  }
  stream << "]";
  return stream;
};

template <typename T, SIMDWidth width_V>
std::string asHex(const AlignedArray<T, width_V>& vec)
{
  std::stringstream ss;
  ss << "[";
  for (const auto& elem : vec) {
    ss << fmt::format("{:#0x}, ", elem);
  }
  ss << "]";
  return ss.str();
};

} // namespace simd
} // namespace internal
} // namespace rans
} // namespace o2

#endif /* RANS_INTERNAL_SIMD_TYPES_H */