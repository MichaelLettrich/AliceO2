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

#include <gsl/span>

#include "rANS/internal/backend/simd/utils.h"
#include "rANS/internal/helper.h"

namespace o2::rans::internal::simd::impl
{
class IdentityFormatingFunctor
{
 public:
  template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
  inline std::string operator()(const T& value)
  {
    return fmt::format("{}", value);
  }

  inline std::string operator()(const uint8_t& value)
  {
    return fmt::format("{}", static_cast<uint32_t>(value));
  }
};

class HexFormatingFunctor
{
 public:
  template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
  inline std::string operator()(const T& value)
  {
    return fmt::format("{:#0x}, ", value);
  }
};
} // namespace o2::rans::internal::simd::impl

namespace std
{
template <typename T, size_t extent_V, class formatingFunctor = o2::rans::internal::simd::impl::IdentityFormatingFunctor>
std::ostream& operator<<(std::ostream& stream, const gsl::span<T, extent_V>& span)
{

  if (span.empty()) {
    stream << "[]";
    return stream;
  } else {
    formatingFunctor formater;

    stream << "[";
    for (size_t i = 0; i < span.size() - 1; ++i) {
      stream << formater(span[i]) << ", ";
    }
    stream << formater(*(--span.end())) << "]";
    return stream;
  }
  return stream;
}
} // namespace std

namespace o2
{
namespace rans
{
namespace internal
{
namespace simd
{

enum class SIMDWidth : uint32_t { SSE = 128u,
                                  AVX = 256u };

inline constexpr size_t getLaneWidthBits(SIMDWidth width) noexcept { return static_cast<size_t>(width); };

inline constexpr size_t getLaneWidthBytes(SIMDWidth width) noexcept { return toBytes(static_cast<size_t>(width)); };

inline constexpr size_t getAlignment(SIMDWidth width) noexcept { return getLaneWidthBytes(width); };

template <class T, size_t N>
inline constexpr T* assume_aligned(T* ptr) noexcept
{
  return reinterpret_cast<T*>(__builtin_assume_aligned(ptr, N, 0));
};

template <class T, SIMDWidth width_V>
inline constexpr T* assume_aligned(T* ptr) noexcept
{
  constexpr size_t alignment = getAlignment(width_V);
  return assume_aligned<T, alignment>(ptr);
};

template <typename T, SIMDWidth width_V>
inline constexpr bool isAligned(T* ptr)
{
  // only aligned iff ptr is divisible by alignment
  constexpr size_t alignment = getAlignment(width_V);
  return !(reinterpret_cast<uintptr_t>(ptr) % alignment);
};

template <typename T>
inline constexpr size_t getElementCount(SIMDWidth width) noexcept
{
  return getLaneWidthBytes(width) / sizeof(T);
};

template <typename T>
inline constexpr SIMDWidth getSimdWidth(size_t nHardwareStreams) noexcept
{
  return static_cast<SIMDWidth>(nHardwareStreams * toBits(sizeof(T)));
};

template <typename array_T>
class AlignedArrayIterator
{
  template <typename T>
  struct getValueType {
    using type = std::conditional_t<std::is_const_v<T>, const typename T::value_type, typename T::value_type>;
  };

 public:
  using difference_type = std::ptrdiff_t;
  using value_type = gsl::span<typename getValueType<array_T>::type, array_T::nElementsPerLane()>;
  using pointer = value_type*;
  using reference = value_type&;
  using iterator_category = std::random_access_iterator_tag;

  inline constexpr AlignedArrayIterator() noexcept = default;

  inline constexpr AlignedArrayIterator(array_T* array, size_t index) noexcept : mIndex{index} {};
  inline constexpr AlignedArrayIterator(const AlignedArrayIterator& iter) noexcept = default;
  inline constexpr AlignedArrayIterator(AlignedArrayIterator&& iter) noexcept = default;
  inline constexpr AlignedArrayIterator& operator=(const AlignedArrayIterator& other) noexcept = default;
  inline constexpr AlignedArrayIterator& operator=(AlignedArrayIterator&& other) noexcept = default;
  inline ~AlignedArrayIterator() noexcept = default;

  // pointer arithmetics
  inline constexpr AlignedArrayIterator& operator++() noexcept
  {
    ++mIndex;
    return *this;
  };

  inline constexpr AlignedArrayIterator operator++(int) noexcept
  {
    auto res = *this;
    ++(*this);
    return res;
  };

  inline constexpr AlignedArrayIterator& operator--() noexcept
  {
    --mIndex;
    return *this;
  };

  inline constexpr AlignedArrayIterator operator--(int) noexcept
  {
    auto res = *this;
    --(*this);
    return res;
  };

  inline constexpr AlignedArrayIterator& operator+=(difference_type i) noexcept
  {
    mIndex += i;
    return *this;
  };

  inline constexpr AlignedArrayIterator operator+(difference_type i) const noexcept
  {
    auto tmp = *const_cast<AlignedArrayIterator*>(this);
    return tmp += i;
  }

  inline constexpr AlignedArrayIterator& operator-=(difference_type i) noexcept
  {
    mIndex -= i;
    return *this;
  };

  inline constexpr AlignedArrayIterator operator-(difference_type i) const noexcept
  {
    auto tmp = *const_cast<AlignedArrayIterator*>(this);
    return tmp -= i;
  };

  inline constexpr difference_type operator-(const AlignedArrayIterator& other) const noexcept
  {
    return this->mIter - other.mIter;
  };

  // comparison
  inline constexpr bool operator==(const AlignedArrayIterator& other) const noexcept { return this->mIndex == other.mIndex; };
  inline constexpr bool operator!=(const AlignedArrayIterator& other) const noexcept { return this->mIndex != other.mIndex; };
  inline constexpr bool operator<(const AlignedArrayIterator& other) const noexcept { return this->mIndex < other->mIndex; };
  inline constexpr bool operator>(const AlignedArrayIterator& other) const noexcept { return this->mIndex > other->mIndex; };
  inline constexpr bool operator>=(const AlignedArrayIterator& other) const noexcept { return this->mIndex >= other->mIndex; };
  inline constexpr bool operator<=(const AlignedArrayIterator& other) const noexcept { return this->mIndex <= other->mIndex; };

  // dereference
  inline constexpr value_type operator*() const noexcept { return (*mContainer)[mIndex]; };

  inline constexpr value_type operator[](difference_type i) const noexcept { return (*mContainer)[mIndex + i]; };

 private:
  array_T* mContainer;
  size_t mIndex{};
}; // namespace simd

template <typename T, SIMDWidth width_V, size_t size_V = 1>
class alignas(getAlignment(width_V)) AlignedArray
{
 public:
  using value_type = T;
  using pointer = value_type*;
  using reference = value_type&;
  using iterator = AlignedArrayIterator<AlignedArray>;
  using const_iterator = AlignedArrayIterator<const AlignedArray>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  static inline constexpr size_t size() noexcept { return size_V; };
  static inline constexpr size_t nElementsPerLane() noexcept { return getElementCount<value_type>(width_V); };
  static inline constexpr size_t nElements() noexcept { return size() * nElementsPerLane(); };

  inline constexpr AlignedArray() noexcept {}; //NOLINT

  template <typename elem_T, std::enable_if_t<std::is_convertible_v<elem_T, value_type>, bool> = true>
  inline constexpr AlignedArray(elem_T value) noexcept
  {
    for (auto& elem : mData) {
      elem = static_cast<value_type>(value);
    }
  };

  template <typename... Args, std::enable_if_t<(sizeof...(Args) == AlignedArray<T, width_V, size_V>::nElements()) && std::is_convertible_v<std::common_type_t<Args...>, value_type>, bool> = true>
  inline constexpr AlignedArray(Args... args) noexcept : mData{static_cast<value_type>(args)...} {};

  inline constexpr const T* data() const noexcept { return mData; };
  inline constexpr T* data() noexcept { return const_cast<T*>(static_cast<const AlignedArray&>(*this).data()); };
  inline constexpr const_iterator begin() const noexcept { return {this, 0}; };
  inline constexpr const_iterator end() const noexcept { return {this, size()}; };
  inline constexpr iterator begin() noexcept { return {this, 0}; };
  inline constexpr iterator end() noexcept { return {this, size()}; };
  inline constexpr const_reverse_iterator rbegin() const noexcept { return std::reverse_iterator(this->end()); };
  inline constexpr const_reverse_iterator rend() const noexcept { return std::reverse_iterator(this->begin()); };
  inline constexpr reverse_iterator rbegin() noexcept { return std::reverse_iterator(this->end()); };
  inline constexpr reverse_iterator rend() noexcept { return std::reverse_iterator(this->begin()); };

  inline constexpr gsl::span<T, nElementsPerLane()> operator[](size_t idx) { return gsl::span<T, nElementsPerLane()>{mData + idx * nElementsPerLane(), nElementsPerLane()}; };

  inline constexpr gsl::span<const T, nElementsPerLane()> operator[](size_t idx) const { return gsl::span<const T, nElementsPerLane()>{mData + idx * nElementsPerLane(), nElementsPerLane()}; };

  inline constexpr const T& operator()(size_t idx, size_t elem) const
  {
    return (*this)[idx][elem];
  };

  inline constexpr T& operator()(size_t idx, size_t elem) { return const_cast<T&>(static_cast<const AlignedArray&>(*this)(idx, elem)); };

  inline constexpr const T& operator()(size_t idx) const
  {
    return *(this->data() + idx);
  };

  inline constexpr T& operator()(size_t idx) { return const_cast<T&>(static_cast<const AlignedArray&>(*this)(idx)); };

 private:
  T mData[nElements()]{};
};

template <typename T>
struct simdWidth;

template <typename T, SIMDWidth simd_V>
struct simdWidth<AlignedArray<T, simd_V>> : public std::integral_constant<SIMDWidth, simd_V> {
};

template <typename T>
inline constexpr SIMDWidth simdWidth_v = simdWidth<T>::value;

template <typename T>
struct elementCount;

template <typename T, SIMDWidth simd_V>
struct elementCount<AlignedArray<T, simd_V>> : public std::integral_constant<size_t, getElementCount<T>(simd_V)> {
};

template <typename T>
inline constexpr size_t elementCount_v = elementCount<T>::value;

template <SIMDWidth width_V, size_t size_V = 1>
using pd_t = AlignedArray<double_t, width_V, size_V>;
template <SIMDWidth width_V, size_t size_V = 1>
using epi64_t = AlignedArray<uint64_t, width_V, size_V>;
template <SIMDWidth width_V, size_t size_V = 1>
using epi32_t = AlignedArray<uint32_t, width_V, size_V>;
template <SIMDWidth width_V, size_t size_V = 1>
using epi16_t = AlignedArray<uint16_t, width_V, size_V>;
template <SIMDWidth width_V, size_t size_V = 1>
using epi8_t = AlignedArray<uint8_t, width_V, size_V>;

inline constexpr std::uint8_t operator"" _u8(unsigned long long int value) { return static_cast<uint8_t>(value); };
inline constexpr std::int8_t operator"" _i8(unsigned long long int value) { return static_cast<int8_t>(value); };

inline constexpr std::uint16_t operator"" _u16(unsigned long long int value) { return static_cast<uint16_t>(value); };
inline constexpr std::int16_t operator"" _i16(unsigned long long int value) { return static_cast<int16_t>(value); };

template <SIMDWidth>
struct SimdInt;

template <>
struct SimdInt<SIMDWidth::SSE> {
  using value_type = __m128i;
};

template <>
struct SimdInt<SIMDWidth::AVX> {
  using value_type = __m256i;
};

template <SIMDWidth width_V>
using simdI_t = typename SimdInt<width_V>::value_type;

using simdIsse_t = simdI_t<SIMDWidth::SSE>;
using simdIavx_t = simdI_t<SIMDWidth::AVX>;

template <SIMDWidth>
struct SimdDouble;

template <>
struct SimdDouble<SIMDWidth::SSE> {
  using value_type = __m128d;
};

template <>
struct SimdDouble<SIMDWidth::AVX> {
  using value_type = __m256d;
};

template <SIMDWidth width_V>
using simdD_t = typename SimdDouble<width_V>::value_type;

using simdDsse_t = simdD_t<SIMDWidth::SSE>;
using simdDavx_t = simdD_t<SIMDWidth::AVX>;

// alignment atributes cause gcc warnings, but we don't need them, so disable for this specific case.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
template <typename T>
struct toSIMDWidth;

template <>
struct toSIMDWidth<__m128> : public std::integral_constant<SIMDWidth, SIMDWidth::SSE> {
};
template <>
struct toSIMDWidth<__m128i> : public std::integral_constant<SIMDWidth, SIMDWidth::SSE> {
};
template <>
struct toSIMDWidth<__m128d> : public std::integral_constant<SIMDWidth, SIMDWidth::SSE> {
};

template <>
struct toSIMDWidth<__m256> : public std::integral_constant<SIMDWidth, SIMDWidth::AVX> {
};
template <>
struct toSIMDWidth<__m256i> : public std::integral_constant<SIMDWidth, SIMDWidth::AVX> {
};
template <>
struct toSIMDWidth<__m256d> : public std::integral_constant<SIMDWidth, SIMDWidth::AVX> {
};

template <typename T>
inline constexpr SIMDWidth toSIMDWidth_v = toSIMDWidth<T>::value;

#pragma GCC diagnostic pop

template <typename T, SIMDWidth width_V, size_t size_V, class formater_T = impl::IdentityFormatingFunctor>
std::ostream& operator<<(std::ostream& stream, const AlignedArray<T, width_V, size_V>& array)
{
  stream << "[";
  for (auto subspan : array) {
    operator<<<const T, getElementCount<T>(width_V), formater_T>(stream, subspan);
    stream << ", ";
  }
  stream << "]";
  return stream;
};

template <typename T, SIMDWidth width_V, size_t size_V>
std::string asHex(const AlignedArray<T, width_V, size_V>& array)
{
  std::ostringstream stream;
  operator<<<T, width_V, size_V, impl::HexFormatingFunctor>(stream, array);
  return stream.str();
};

} // namespace simd
} // namespace internal
} // namespace rans
} // namespace o2

namespace gsl
{
template <typename T, o2::rans::internal::simd::SIMDWidth width_V, size_t size_V>
auto make_span(const o2::rans::internal::simd::AlignedArray<T, width_V, size_V>& array)
{
  return gsl::span<const T, o2::rans::internal::simd::AlignedArray<T, width_V, size_V>::nElements()>(array.data(), array.nElements());
}

template <typename T, o2::rans::internal::simd::SIMDWidth width_V, size_t size_V>
auto make_span(o2::rans::internal::simd::AlignedArray<T, width_V, size_V>& array)
{
  return gsl::span<T, o2::rans::internal::simd::AlignedArray<T, width_V, size_V>::nElements()>(array.data(), array.nElements());
}

} // namespace gsl

#endif /* RANS_INTERNAL_SIMD_TYPES_H */