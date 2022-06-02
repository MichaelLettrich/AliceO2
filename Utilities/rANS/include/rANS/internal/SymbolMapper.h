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

/// @file   EncoderBase.h
/// @author michael.lettrich@cern.ch
/// @since  Feb 8, 2021
/// @brief

#ifndef INCLUDE_RANS_SYMBOLMAPPER_H_
#define INCLUDE_RANS_SYMBOLMAPPER_H_

#include <algorithm>
#include <iomanip>
#include <memory>

#include <fairlogger/Logger.h>

#include "rANS/definitions.h"
#include "rANS/internal/helper.h"
#include "rANS/internal/SIMDEncodeCommand.h"

namespace o2
{
namespace rans
{
namespace internal
{

template <typename symbolTable_T, typename coder_T, typename incompressible_IT, typename derived_T>
class SymbolMapperIterface
{
 public:
  using symbolTable_type = symbolTable_T;
  using coder_type = coder_T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using source_type = typename symbolTable_type::source_type;
  using symbol_type = typename symbolTable_type::symbol_type;
  using coderSymbol_type = typename coder_type::symbol_type;
  using incompressible_iterator = incompressible_IT;

  // static_assert(std::is_same_v<typename symbolTable_type::value_type, typename coder_type::symbol_type>);

  template <typename source_IT>
  [[nodiscard]] source_IT unpackSymbols(source_IT sourceIter, coderSymbol_type& symbol, size_type nStreams)
  {
    return static_cast<derived_T*>(this)->unpackSymbols(sourceIter, symbol, nStreams);
  };

  template <typename source_IT>
  [[nodiscard]] source_IT unpackSymbols(source_IT sourceIter, coderSymbol_type& symbol)
  {
    return static_cast<derived_T*>(this)->unpackSymbols(sourceIter, symbol);
  };

  incompressible_iterator getIncompressibleIterator() const { return mIncompressibleIter; };

 protected:
  template <typename source_IT>
  [[nodiscard]] inline const symbol_type& lookupSymbol(source_IT sourceIter)
  {
    LOGP(info, "unpacking {}", fmt::ptr(sourceIter));

    const symbol_type& symbol = (*mSymbolTable)[*sourceIter];

    if constexpr (!std::is_null_pointer_v<incompressible_iterator>) {
      if (mSymbolTable->isEscapeSymbol(symbol)) {
        *mIncompressibleIter++ = *sourceIter;
      }
    }
    return symbol;
  };

  SymbolMapperIterface() = default;
  SymbolMapperIterface(const symbolTable_type& symbolTable,
                       incompressible_IT incompressibleIter = nullptr) : mSymbolTable{&symbolTable},
                                                                         mIncompressibleIter{incompressibleIter} {};

  const symbolTable_type* mSymbolTable{};
  incompressible_iterator mIncompressibleIter{};
}; // namespace SymbolMapperIterface

template <typename symbolTable_T, typename coder_T, typename incompressible_IT = std::nullptr_t>
class SymbolMapper : public SymbolMapperIterface<symbolTable_T,
                                                 coder_T,
                                                 incompressible_IT,
                                                 SymbolMapper<symbolTable_T, coder_T, incompressible_IT>>
{
  using base_type = SymbolMapperIterface<symbolTable_T, coder_T, incompressible_IT, SymbolMapper<symbolTable_T, coder_T, incompressible_IT>>;

 public:
  using symbolTable_type = typename base_type::symbolTable_type;
  using coder_type = typename base_type::coder_type;
  using size_type = typename base_type::size_type;
  using difference_type = typename base_type::difference_type;
  using source_type = typename base_type::source_type;
  using symbol_type = typename base_type::symbol_type;
  using coderSymbol_type = typename base_type::coderSymbol_type;
  using incompressible_iterator = typename base_type::incompressible_iterator;

  static_assert(coder_type::getNstreams() == 1);

  SymbolMapper() = default;

  SymbolMapper(const symbolTable_type& symbolTable, incompressible_IT incompressibleIter = nullptr) : base_type{symbolTable, incompressibleIter} {};

  template <typename source_IT>
  [[nodiscard]] inline source_IT unpackSymbols(source_IT sourceIter, coderSymbol_type& unpacked)
  {
    unpacked = this->lookupSymbol(sourceIter);
    return --sourceIter;
  };

  template <typename source_IT>
  [[nodiscard]] inline source_IT unpackSymbols(source_IT sourceIter, coderSymbol_type& unpacked, size_type nStreams)
  {
    return this->unpackSymbols(sourceIter, unpacked);
  };
};

template <typename symbolTable_T, typename incompressible_IT>
class SymbolMapper<symbolTable_T,
                   SSEEncoderCommand<20>,
                   incompressible_IT> : public SymbolMapperIterface<symbolTable_T,
                                                                    SSEEncoderCommand<20>,
                                                                    incompressible_IT,
                                                                    SymbolMapper<symbolTable_T, incompressible_IT>>
{
  using base_type = SymbolMapperIterface<symbolTable_T, SSEEncoderCommand<20>, incompressible_IT, SymbolMapper<symbolTable_T, incompressible_IT>>;

 public:
  using symbolTable_type = typename base_type::symbolTable_type;
  using coder_type = typename base_type::coder_type;
  using size_type = typename base_type::size_type;
  using difference_type = typename base_type::difference_type;
  using source_type = typename base_type::source_type;
  using symbol_type = typename base_type::symbol_type;
  using coderSymbol_type = typename base_type::coderSymbol_type;
  using incompressible_iterator = typename base_type::incompressible_iterator;

  static_assert(coder_type::getNstreams() == 4);

  SymbolMapper() = default;

  SymbolMapper(const symbolTable_type& symbolTable, incompressible_IT incompressibleIter = nullptr) : base_type{symbolTable, incompressibleIter} {};

  template <typename source_IT>
  [[nodiscard]] inline source_IT unpackSymbols(source_IT sourceIter, coderSymbol_type& unpacked)
  {

    auto unpacker = [&, this](auto srcIter, size_t dstIdx) {
      auto& symbol = this->lookupSymbol(srcIter);
      unpacked.frequencies[dstIdx] = symbol.getFrequency();
      unpacked.cumulativeFrequencies[dstIdx] = symbol.getCumulative();
    };

    using namespace simd;

    unpacker(sourceIter - 0, 5);
    unpacker(sourceIter - 1, 4);
    unpacker(sourceIter - 2, 1);
    unpacker(sourceIter - 3, 0);

    return internal::advanceIter(sourceIter, -coder_type::getNstreams());
  };

  template <typename source_IT>
  [[nodiscard]] inline source_IT unpackSymbols(source_IT sourceIter, coderSymbol_type& unpacked, size_type nActiveStreams)
  {
    const size_type nStreams = coder_type::getNstreams();

    difference_type currentStream = nActiveStreams;

    while (currentStream-- > 0) {
      const auto& symbol = this->lookupSymbol(sourceIter--);
      unpacked.frequencies[currentStream] = symbol.getFrequency();
      unpacked.cumulativeFrequencies[currentStream] = symbol.getCumulative();
    }

    return sourceIter;
  };
}; // namespace internal

template <typename symbolTable_T, typename incompressible_IT>
class SymbolMapper<symbolTable_T,
                   AVXEncoderCommand<20>,
                   incompressible_IT> : public SymbolMapperIterface<symbolTable_T,
                                                                    AVXEncoderCommand<20>,
                                                                    incompressible_IT,
                                                                    SymbolMapper<symbolTable_T, incompressible_IT>>
{
  using base_type = SymbolMapperIterface<symbolTable_T, AVXEncoderCommand<20>, incompressible_IT, SymbolMapper<symbolTable_T, incompressible_IT>>;

 public:
  using symbolTable_type = typename base_type::symbolTable_type;
  using coder_type = typename base_type::coder_type;
  using size_type = typename base_type::size_type;
  using difference_type = typename base_type::difference_type;
  using source_type = typename base_type::source_type;
  using symbol_type = typename base_type::symbol_type;
  using coderSymbol_type = typename base_type::coderSymbol_type;
  using incompressible_iterator = typename base_type::incompressible_iterator;

  static_assert(coder_type::getNstreams() == 8);

  SymbolMapper() = default;

  SymbolMapper(const symbolTable_type& symbolTable, incompressible_IT incompressibleIter = nullptr) : base_type{symbolTable, incompressibleIter} {};

  template <typename source_IT>
  [[nodiscard]] inline source_IT unpackSymbols(source_IT sourceIter, coderSymbol_type& unpacked)
  {

    auto unpacker = [&, this](auto srcIter, size_t dstIdx) {
      auto& symbol = this->lookupSymbol(srcIter);
      unpacked.frequencies[dstIdx] = symbol.getFrequency();
      unpacked.cumulativeFrequencies[dstIdx] = symbol.getCumulative();
    };

    using namespace simd;

    unpacker(sourceIter - 0, 7);
    unpacker(sourceIter - 1, 6);
    unpacker(sourceIter - 2, 5);
    unpacker(sourceIter - 3, 4);
    unpacker(sourceIter - 4, 3);
    unpacker(sourceIter - 5, 2);
    unpacker(sourceIter - 6, 1);
    unpacker(sourceIter - 7, 0);

    return internal::advanceIter(sourceIter, -coder_type::getNstreams());
  };

  template <typename source_IT>
  [[nodiscard]] inline source_IT unpackSymbols(source_IT sourceIter, coderSymbol_type& unpacked, size_type nActiveStreams)
  {
    const size_type nStreams = coder_type::getNstreams();

    difference_type currentStream = nActiveStreams;

    while (currentStream-- > 0) {
      const auto& symbol = this->lookupSymbol(sourceIter--);
      unpacked.frequencies[currentStream] = symbol.getFrequency();
      unpacked.cumulativeFrequencies[currentStream] = symbol.getCumulative();
    }

    return sourceIter;
  };
};

}; // namespace internal
}; // namespace rans
}; // namespace o2

#endif /* INCLUDE_RANS_SYMBOLMAPPER_H_ */