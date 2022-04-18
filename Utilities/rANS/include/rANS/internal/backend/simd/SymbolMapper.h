// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   SymbolTable.h
/// @author Michael Lettrich
/// @since  2019-06-21
/// @brief  Container for information needed to encode/decode a symbol of the alphabet

#ifndef RANS_INTERNAL_SIMD_SYMBOLMAPPER_H
#define RANS_INTERNAL_SIMD_SYMBOLMAPPER_H

#include "rANS/internal/backend/simd/Symbol.h"
#include "rANS/internal/backend/simd/SymbolTable.h"
#include "rANS/internal/backend/simd/types.h"
#include "rANS/internal/backend/simd/kernel.h"
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

// struct UnrolledSymbols {
//   epi32_t<SIMDWidth::SSE, 2> frequencies;
//   epi32_t<SIMDWidth::SSE, 2> cumulativeFrequencies;
// };

// template <typename source_IT>
// inline const internal::simd::Symbol* lookupSymbol(source_IT iter, const simd::SymbolTable& symbolTable, std::vector<typename std::iterator_traits<source_IT>::value_type>& literals) noexcept
// {
//   const auto symbol = *iter;
//   const auto* encoderSymbol = &(symbolTable[symbol]);
//   if (symbolTable.isEscapeSymbol(*encoderSymbol)) {
//     literals.push_back(symbol);
//   }
//   return encoderSymbol;
// };

// template <typename source_IT, SIMDWidth width_V>
// std::tuple<source_IT, UnrolledSymbols> getSymbols(source_IT symbolIter, const o2::rans::internal::simd::SymbolTable& symbolTable, std::vector<typename std::iterator_traits<source_IT>::value_type>& literals)
// {
//   UnrolledSymbols unrolledSymbols;

//   if constexpr (width_V == SIMDWidth::SSE) {
//     AlignedArray<const Symbol*, simd::SIMDWidth::SSE, 4> ret;
//     ret[3] = lookupSymbol(symbolIter - 1, symbolTable, literals);
//     ret[2] = lookupSymbol(symbolIter - 2, symbolTable, literals);
//     ret[1] = lookupSymbol(symbolIter - 3, symbolTable, literals);
//     ret[0] = lookupSymbol(symbolIter - 4, symbolTable, literals);

//     AlignedArray<symbol_t, SIMDWidth::SSE, 4> symbols{
//       static_cast<symbol_t>(*(symbolIter - 4)),
//       static_cast<symbol_t>(*(symbolIter - 3)),
//       static_cast<symbol_t>(*(symbolIter - 2)),
//       static_cast<symbol_t>(*(symbolIter - 1)),
//     };

//     const __m128i symbolsVec = load(toConstSIMDView(symbols));
//     const __m128i minVec = _mm_set1_epi32(symbolTable.getMinSymbol());
//     const __m128i maxVec = _mm_set1_epi32(symbolTable.getMaxSymbol());
//     const __m128i escapeIdxVec = _mm_set1_epi32(symbolTable.size() - 1);

//     // mask
//     const __m128i isGT = _mm_cmpgt_epi32(symbolsVec, maxVec);
//     const __m128i isLT = _mm_cmplt_epi32(symbolsVec, minVec);
//     const __m128i isOutOfRange = _mm_or_si128(isGT, isLT);

//     const __m128i offsets = _mm_blendv_epi8(_mm_sub_epi32(symbolsVec, minVec), escapeIdxVec, isOutOfRange);

//     LOG(info) << offsets << store<uint32_t>(offsets);

//     aosToSoa(ArrayView{ret}.template subView<0, 2>(),
//              toSIMDView(unrolledSymbols.frequencies).template subView<0, 1>(),
//              toSIMDView(unrolledSymbols.cumulativeFrequencies).template subView<0, 1>());
//     aosToSoa(ArrayView{ret}.template subView<2, 2>(),
//              toSIMDView(unrolledSymbols.frequencies).template subView<1, 1>(),
//              toSIMDView(unrolledSymbols.cumulativeFrequencies).template subView<1, 1>());
//     //aosToSoa(ret, toSIMDView(unrolledSymbols.frequencies), toSIMDView(unrolledSymbols.cumulativeFrequencies));
//     return {symbolIter - 4, unrolledSymbols};
//   } else {
//     AlignedArray<const Symbol*, simd::SIMDWidth::SSE, 8> ret;
//     ret[7] = lookupSymbol(symbolIter - 1, symbolTable, literals);
//     ret[6] = lookupSymbol(symbolIter - 2, symbolTable, literals);
//     ret[5] = lookupSymbol(symbolIter - 3, symbolTable, literals);
//     ret[4] = lookupSymbol(symbolIter - 4, symbolTable, literals);
//     ret[3] = lookupSymbol(symbolIter - 5, symbolTable, literals);
//     ret[2] = lookupSymbol(symbolIter - 6, symbolTable, literals);
//     ret[1] = lookupSymbol(symbolIter - 7, symbolTable, literals);
//     ret[0] = lookupSymbol(symbolIter - 8, symbolTable, literals);

//     aosToSoa(ArrayView{ret}.template subView<0, 4>(),
//              toSIMDView(unrolledSymbols.frequencies).template subView<0, 1>(),
//              toSIMDView(unrolledSymbols.cumulativeFrequencies).template subView<0, 1>());
//     aosToSoa(ArrayView{ret}.template subView<4, 4>(),
//              toSIMDView(unrolledSymbols.frequencies).template subView<1, 1>(),
//              toSIMDView(unrolledSymbols.cumulativeFrequencies).template subView<1, 1>());
//     return {symbolIter - 8, unrolledSymbols};
//   }
// };

inline constexpr std::array<epi8_t<SIMDWidth::SSE>, 16>
  SSEIncompressibleMapping{{
    {0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0000   0xFFFFu
    {0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0001   0x0FFFu
    {0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0010   0x1FFFu
    {0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0011   0x10FFu
    {0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0100   0x2FFFu
    {0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0101   0x20FFu
    {0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0110   0x21FFu
    {0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b0111   0x210Fu
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b1000   0x3FFFu
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b1001   0x30FFu
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b1010   0x31FFu
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b1011   0x310Fu
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b1100   0x32FFu
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b1101   0x320Fu
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8, 0xFF_u8}, //0b1110   0x321Fu
    {0x0C_u8, 0x0D_u8, 0x0E_u8, 0x0F_u8, 0x08_u8, 0x09_u8, 0x0A_u8, 0x0B_u8, 0x04_u8, 0x05_u8, 0x06_u8, 0x07_u8, 0x00_u8, 0x01_u8, 0x02_u8, 0x03_u8}  //0b1111   0x3210u
  }};

template <SIMDWidth simdWidth_V>
class SymbolMapper;

template <>
class SymbolMapper<SIMDWidth::SSE>
{

 public:
  SymbolMapper(const SymbolTable& symbolTable) : mSymbolTable(&symbolTable),
                                                 mMinVec(_mm_set1_epi32(mSymbolTable->getMinSymbol())),
                                                 mMaxVec(_mm_set1_epi32(mSymbolTable->getMaxSymbol())),
                                                 mEscapeIdxVec(_mm_set1_epi32(mSymbolTable->size() - 1)),
                                                 mIncompressibleCumulatedFrequencies(_mm_set1_epi32(mSymbolTable->getEscapeSymbol().getCumulative())){};

  //   template <typename source_IT>
  //   std::tuple<source_IT, UnrolledSymbols> readSymbols(source_IT symbolIter);
  template <typename source_IT>
  std::tuple<source_IT, UnrolledSymbols> readSymbols(source_IT symbolIter, std::vector<typename std::iterator_traits<source_IT>::value_type>& literals);

 private:
  const SymbolTable* mSymbolTable{};
  __m128i mMinVec;
  __m128i mMaxVec;
  __m128i mEscapeIdxVec;
  __m128i mIncompressibleCumulatedFrequencies;
};

template <typename source_IT>
inline auto SymbolMapper<SIMDWidth::SSE>::readSymbols(source_IT symbolIter, std::vector<typename std::iterator_traits<source_IT>::value_type>& literals) -> std::tuple<source_IT, UnrolledSymbols>
{
  using source_t = typename std::iterator_traits<source_IT>::value_type;

  UnrolledSymbols unrolledSymbols;

  __m128i symbolsVec;

  if constexpr (std::is_pointer_v<source_IT>) {

    symbolsVec = _mm_loadu_si128(reinterpret_cast<const __m128i*>(symbolIter - 4));
    if constexpr (std::is_same_v<source_t, uint8_t>) {
      symbolsVec = _mm_cvtepu8_epi32(symbolsVec);
    } else if constexpr (std::is_same_v<source_t, int8_t>) {
      symbolsVec = _mm_cvtepi8_epi32(symbolsVec);
    } else if constexpr (std::is_same_v<source_t, uint16_t>) {
      symbolsVec = _mm_cvtepu16_epi32(symbolsVec);
    } else if constexpr (std::is_same_v<source_t, int16_t>) {
      symbolsVec = _mm_cvtepi16_epi32(symbolsVec);
    }
    //no conversion needed for int32
  } else {
    AlignedArray<symbol_t, SIMDWidth::SSE, 4> symbols{
      static_cast<symbol_t>(*(symbolIter - 4)),
      static_cast<symbol_t>(*(symbolIter - 3)),
      static_cast<symbol_t>(*(symbolIter - 2)),
      static_cast<symbol_t>(*(symbolIter - 1)),
    };
    symbolsVec = load(toConstSIMDView(symbols));
  }

  // range check
  const __m128i isGT = _mm_cmpgt_epi32(symbolsVec, mMaxVec);
  const __m128i isLT = _mm_cmplt_epi32(symbolsVec, mMinVec);
  const __m128i isOutOfRange = _mm_or_si128(isGT, isLT);
  // make sure we're in the right range
  const __m128i offsetsVec = _mm_blendv_epi8(_mm_sub_epi32(symbolsVec, mMinVec), mEscapeIdxVec, isOutOfRange);

  auto offset = store<uint32_t>(offsetsVec);

  __m128i symbol0 = _mm_load_si128(reinterpret_cast<__m128i const*>(&mSymbolTable->at(offset[0])));
  __m128i symbol1 = _mm_load_si128(reinterpret_cast<__m128i const*>(&mSymbolTable->at(offset[1])));
  __m128i symbol2 = _mm_load_si128(reinterpret_cast<__m128i const*>(&mSymbolTable->at(offset[2])));
  __m128i symbol3 = _mm_load_si128(reinterpret_cast<__m128i const*>(&mSymbolTable->at(offset[3])));

  // Unpack Symbols
  __m128i merged0 = _mm_unpacklo_epi32(symbol0, symbol1);
  __m128i merged1 = _mm_unpacklo_epi32(symbol2, symbol3);
  __m128i frequencies = _mm_unpacklo_epi64(merged0, merged1);
  __m128i cumulativeFrequencies = _mm_unpackhi_epi64(merged0, merged1);

  // find all incompressibleSymbols and pass them to a literals vector
  const __m128i isIncompressible = _mm_cmpeq_epi32(cumulativeFrequencies, mIncompressibleCumulatedFrequencies);

  const uint32_t id = _mm_movemask_ps(_mm_castsi128_ps(isIncompressible));

  const uint32_t nIncompressible = _mm_popcnt_u32(id);
  // if (nIncompressible > 0) {
  //   LOG(info) << "Cumul " << store<uint32_t>(mIncompressibleCumulatedFrequencies) << " vs " << store<uint32_t>(cumulativeFrequencies);
  //   LOG(info) << "Incompressible: " << asHex(store<uint32_t>(isIncompressible)) << "at " << store<uint32_t>(isIncompressible);
  // }

  __m128i shuffleMaskVec = load(toConstSIMDView(SSEIncompressibleMapping[id]));
  symbolsVec = _mm_shuffle_epi8(symbolsVec, shuffleMaskVec);

  UnrolledSymbols unrolledSymbols2;

  store(frequencies, toSIMDView(unrolledSymbols2.frequencies).template subView<0, 1>());
  store(_mm_bsrli_si128(frequencies, 8), toSIMDView(unrolledSymbols2.frequencies).template subView<1, 1>());

  store(cumulativeFrequencies, toSIMDView(unrolledSymbols2.cumulativeFrequencies).template subView<0, 1>());
  store(_mm_bsrli_si128(cumulativeFrequencies, 8), toSIMDView(unrolledSymbols2.cumulativeFrequencies).template subView<1, 1>());

  // store;
  auto incompressibleSymbols = store<uint32_t>(symbolsVec);

  // std::vector<typename std::iterator_traits<source_IT>::value_type> fakeLiterals;
  // fakeLiterals.reserve(4);

  // AlignedArray<const Symbol*, simd::SIMDWidth::SSE, 4>
  //   ret;
  // ret[3] = lookupSymbol(symbolIter - 1, *mSymbolTable, fakeLiterals);
  // ret[2] = lookupSymbol(symbolIter - 2, *mSymbolTable, fakeLiterals);
  // ret[1] = lookupSymbol(symbolIter - 3, *mSymbolTable, fakeLiterals);
  // ret[0] = lookupSymbol(symbolIter - 4, *mSymbolTable, fakeLiterals);

  // aosToSoa(ArrayView{ret}.template subView<0, 2>(),
  //          toSIMDView(unrolledSymbols.frequencies).template subView<0, 1>(),
  //          toSIMDView(unrolledSymbols.cumulativeFrequencies).template subView<0, 1>());
  // aosToSoa(ArrayView{ret}.template subView<2, 2>(),
  //          toSIMDView(unrolledSymbols.frequencies).template subView<1, 1>(),
  //          toSIMDView(unrolledSymbols.cumulativeFrequencies).template subView<1, 1>());
  // //aosToSoa(ret, toSIMDView(unrolledSymbols.frequencies), toSIMDView(unrolledSymbols.cumulativeFrequencies));

  // auto checkEqual = [](auto a, auto b) {
  //   if (a != b) {
  //     LOGP(warning, "{}!={}", a, b);
  //   }
  // };

  // // checks
  // //LOG(info) << "frequency check";
  // checkEqual(unrolledSymbols2.frequencies[0], unrolledSymbols.frequencies[0]);
  // checkEqual(unrolledSymbols2.frequencies[1], unrolledSymbols.frequencies[1]);
  // checkEqual(unrolledSymbols2.frequencies[4], unrolledSymbols.frequencies[4]);
  // checkEqual(unrolledSymbols2.frequencies[5], unrolledSymbols.frequencies[5]);

  // //LOG(info) << "cumulativeFrequencies check";
  // checkEqual(unrolledSymbols2.cumulativeFrequencies[0], unrolledSymbols.cumulativeFrequencies[0]);
  // checkEqual(unrolledSymbols2.cumulativeFrequencies[1], unrolledSymbols.cumulativeFrequencies[1]);
  // checkEqual(unrolledSymbols2.cumulativeFrequencies[4], unrolledSymbols.cumulativeFrequencies[4]);
  // checkEqual(unrolledSymbols2.cumulativeFrequencies[5], unrolledSymbols.cumulativeFrequencies[5]);

  // // LOG(info) << "checking incompressible sizes";
  // checkEqual(_mm_popcnt_u32(id), fakeLiterals.size());

  // // LOG(info) << "incompressible symbols check";
  // for (size_t i = 0; i < fakeLiterals.size(); ++i) {
  //   checkEqual(incompressibleSymbols[i], fakeLiterals[i]);
  // }

  for (size_t i = 0; i < nIncompressible; ++i) {
    literals.push_back(static_cast<typename std::iterator_traits<source_IT>::value_type>(incompressibleSymbols[i]));
  }

  // for (auto& s : fakeLiterals) {
  //   literals.push_back(s);
  // }

  return {symbolIter - 4, unrolledSymbols2};
};

} // namespace simd
} // namespace internal
} // namespace rans
} // namespace o2
#endif /* RANS_INTERNAL_SIMD_SYMBOLMAPPER_H */
