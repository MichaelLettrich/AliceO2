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

#include <rANS/internal/helper.h>

namespace o2
{
namespace rans
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

inline void pack(uint64_t*& __restrict__ buffer, size_t& bitOffest, uint64_t data, size_t packingWidth)
{
  constexpr size_t I64Bits = internal::toBits(sizeof(uint64_t));

  const size_t newBitOffset = bitOffest + packingWidth;
  if (newBitOffset <= I64Bits) {
    *buffer |= data << bitOffest;
    bitOffest = newBitOffset;
  } else {
    const size_t leftTail = I64Bits - bitOffest;
    *buffer |= (data & All1BackTill[leftTail]) << bitOffest;
    *(++buffer) |= (data >> leftTail);
    bitOffest = newBitOffset - I64Bits;
  }
};

template <typename source_IT>
inline uint64_t* packStream(const source_IT srcBufferBegin, const source_IT srcBufferEnd, uint64_t* __restrict__ dstBufferBegin, size_t packingWidth)
{
  auto dstIter = dstBufferBegin;
  size_t pos = 0;
  for (auto srcIter = srcBufferBegin; srcIter != srcBufferEnd; ++srcIter) {
    o2::rans::pack(dstIter, pos, static_cast<uint64_t>(*srcIter), packingWidth);
  };
  return ++dstIter;
};

inline uint64_t unpackNext(const uint64_t*& __restrict__ buffer, size_t& bitOffest, size_t packingWidth)
{
  constexpr size_t I64Bits = internal::toBits(sizeof(uint64_t));
  const size_t bitEnd = bitOffest + packingWidth;

  if (bitEnd <= I64Bits) {
    const uint64_t mask = All1FrontTill[bitOffest] & All1BackTill[bitEnd];
    uint64_t ret = (*buffer & mask) >> bitOffest;
    bitOffest = bitEnd;
    return ret;
  } else {
    // first part
    uint64_t mask = All1FrontTill[bitOffest];
    uint64_t value = (*buffer & mask) >> bitOffest;
    //second part
    bitOffest = bitEnd - I64Bits;
    mask = All1BackTill[bitOffest];
    value |= (*(++buffer) & mask) << (packingWidth - bitOffest);
    return value;
  }
};

inline uint64_t unpackByIndex(const uint64_t* __restrict__ buffer, size_t index, size_t packingWidth)
{
  constexpr size_t I64Bits = internal::toBits(sizeof(uint64_t));
  constexpr size_t I64Log2 = 6;

  const size_t StartBitPos = index * packingWidth;
  const size_t ArrayIdx = StartBitPos >> I64Log2; // integer div pow2
  size_t beginElem = StartBitPos & (I64Bits - 1); // integer mod pow2

  const uint64_t* iter = buffer + ArrayIdx;

  return unpackNext(iter, beginElem, packingWidth);
};

template <typename source_T>
inline void unpackStream(const uint64_t* __restrict__ packingBufferBegin, source_T* __restrict__ destBufferBegin, size_t nElements, size_t packingWidth)
{
  auto iter = packingBufferBegin;
  size_t bitOffset = 0;

  for (size_t i = 0; i < nElements; ++i) {
    *(destBufferBegin++) = static_cast<source_T>(unpackNext(iter, bitOffset, packingWidth));
  }
};

} // namespace rans
} // namespace o2

#endif /* RANS_PACKING_H */