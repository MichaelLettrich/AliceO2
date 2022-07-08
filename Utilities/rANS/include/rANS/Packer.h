// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   Packing.h
/// @author Michael Lettrich
/// @since  2021-03-18
/// @brief  Encoder - code symbol into a rANS encoded state

#ifndef RANS_PACKING_H
#define RANS_PACKING_H

#include <cstdint>
#include <cstring>
#include <array>

#include <immintrin.h>

#include <rANS/internal/helper.h>

namespace o2
{
namespace rans
{

namespace internal
{
inline constexpr std::array<uint64_t, 65> All1BackTill = []() constexpr
{
  std::array<uint64_t, 65> ret{};
  for (size_t i = 0; i < 64; ++i) {
    ret[i] = (1ull << i) - 1;
  }
  ret[64] = ~0;
  return ret;
}
();

inline constexpr std::array<uint64_t, 65> All1FrontTill = []() constexpr
{
  std::array<uint64_t, 65> ret{};
  for (size_t i = 0; i < 65; ++i) {
    ret[i] = ~0ull ^ All1BackTill[i];
  }
  return ret;
}
();

[[nodiscard]] inline uint64_t bitExtract(uint64_t data, uint32_t start, uint32_t length) noexcept
{
#ifdef __BMI__
  return _bextr_u64(data, start, length);
#else
  const uint64_t mask = All1BackTill[start + length] ^ All1BackTill[start];
  return (data & mask) >> start;
#endif
};
} // namespace internal

void pack(uint64_t*& __restrict__ buffer, uint32_t& bitOffset, uint64_t data, uint32_t packingWidth)
{
  constexpr uint32_t I64Bits = internal::toBits(sizeof(uint64_t));

  const uint32_t newBitOffset = bitOffset + packingWidth;
  if (newBitOffset <= I64Bits) {
    *buffer |= data << bitOffset;
    bitOffset = newBitOffset;
  } else {
    const uint32_t leftTail = I64Bits - bitOffset;
    *buffer |= internal::bitExtract(data, 0, leftTail) << bitOffset;
    *(++buffer) |= data >> leftTail;
    bitOffset = newBitOffset - I64Bits;
  }
};

template <typename source_IT>
inline uint64_t* packStream(const source_IT srcBufferBegin, const source_IT srcBufferEnd, uint64_t* __restrict__ dstBufferBegin, uint32_t packingWidth)
{
  auto dstIter = dstBufferBegin;
  uint32_t bitOffset = 0;
  for (auto srcIter = srcBufferBegin; srcIter != srcBufferEnd; ++srcIter) {
    o2::rans::pack(dstIter, bitOffset, static_cast<uint64_t>(*srcIter), packingWidth);
  };
  return ++dstIter;
};

[[nodiscard]] inline uint64_t unpackNext(const uint64_t*& __restrict__ buffer, uint32_t& bitOffset, uint32_t packingWidth)
{
  constexpr uint32_t I64Bits = internal::toBits(sizeof(uint64_t));
  const uint32_t bitEnd = bitOffset + packingWidth;
  if (bitEnd <= I64Bits) {
    uint64_t ret = internal::bitExtract(*buffer, bitOffset, packingWidth);
    bitOffset = bitEnd;
    return ret;
  } else {
    // first part
    uint64_t value = internal::bitExtract(*buffer, bitOffset, I64Bits - bitOffset);
    //second part
    bitOffset = bitEnd - I64Bits;
    value |= internal::bitExtract(*(++buffer), 0, bitOffset) << (packingWidth - bitOffset);
    return value;
  }
};

[[nodiscard]] inline uint64_t rUnpackNext(uint64_t*& __restrict__ buffer, uint32_t& bitOffset, uint32_t packingWidth)
{
  constexpr uint32_t I64Bits = internal::toBits(sizeof(uint64_t));

  if (packingWidth < bitOffset) {
    bitOffset -= packingWidth;
    return internal::bitExtract(*buffer, bitOffset, packingWidth);
  } else {
    // first part
    const uint32_t remainderBits = packingWidth - bitOffset;
    uint64_t value = internal::bitExtract(*buffer, 0, bitOffset) << remainderBits;
    //second part
    bitOffset = 64 - remainderBits;
    value |= internal::bitExtract(*(--buffer), bitOffset, remainderBits);
    return value;
  }
};

[[nodiscard]] inline uint64_t unpackByIndex(const uint64_t* __restrict__ buffer, size_t index, uint32_t packingWidth)
{
  constexpr uint32_t I64Bits = internal::toBits(sizeof(uint64_t));
  constexpr uint32_t I64Log2 = 6;

  const size_t StartBitPos = index * packingWidth;
  const size_t ArrayIdx = StartBitPos >> I64Log2;   // integer div pow2
  uint32_t beginElem = StartBitPos & (I64Bits - 1); // integer mod pow2

  const uint64_t* iter = buffer + ArrayIdx;

  return unpackNext(iter, beginElem, packingWidth);
};

template <typename source_T>
inline void unpackStream(const uint64_t* __restrict__ packingBufferBegin, source_T* __restrict__ destBufferBegin, size_t nElements, uint32_t packingWidth)
{
  auto iter = packingBufferBegin;
  uint32_t bitOffset = 0;

  for (size_t i = 0; i < nElements; ++i) {
    *(destBufferBegin++) = static_cast<source_T>(unpackNext(iter, bitOffset, packingWidth));
  }
};

void eliasDeltaEncode(uint64_t*& __restrict__ buffer, uint32_t& bitOffset, uint32_t data)
{
  const uint32_t highestPow2 = internal::log2UIntNZ(data);
  const uint32_t nLeadingZeros = internal::log2UIntNZ(highestPow2 + 1);
  uint64_t value = highestPow2 + 1;
  value = value << highestPow2 | internal::bitExtract(data, 0, highestPow2);
  uint32_t packingOffset = nLeadingZeros + nLeadingZeros + 1 + highestPow2;

  pack(buffer, bitOffset, value, packingOffset);
};

[[nodiscard]] inline uint32_t eliasDeltaDecode(uint64_t*& __restrict__ buffer, uint32_t& bitOffset)
{
  constexpr uint32_t I64Bits = internal::toBits(sizeof(uint64_t));
  constexpr uint32_t maxSymbolBits = 42; // a 32 bit integer can at max use 42 Bits with delta coding;

  // unpack maxSymbolBits bits of data
  uint64_t* unpackIter = buffer;
  uint32_t unpackPos = bitOffset;
  uint64_t rawData = rUnpackNext(unpackIter, unpackPos, maxSymbolBits);

  // do delta decoding algorithm
  rawData <<= I64Bits - maxSymbolBits;
  const uint32_t nLeadingZeros = __builtin_clzl(rawData);
  uint32_t packingOffset = 2 * nLeadingZeros + 1;
  const uint32_t highestPow2 = internal::bitExtract(rawData, I64Bits - packingOffset, nLeadingZeros + 1) - 1;
  packingOffset += highestPow2;
  uint32_t value = (1 << highestPow2) + internal::bitExtract(rawData, I64Bits - packingOffset, highestPow2);

  // set out parameters
  if (packingOffset > bitOffset) {
    bitOffset = I64Bits + bitOffset - packingOffset;
    --buffer;
  } else {
    bitOffset -= packingOffset;
  }
  return value;
};

} // namespace rans
} // namespace o2

#endif /* RANS_PACKING_H */