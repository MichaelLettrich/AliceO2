// Copyright 2019-2023 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   DecoderSymbolTable.h
/// @author Michael Lettrich
/// @brief  Maps rANS state information back to source symbol, used for decoding.

#ifndef RANS_INTERNAL_CONTAINERS_DECODERSYMBOLTABLE_H_
#define RANS_INTERNAL_CONTAINERS_DECODERSYMBOLTABLE_H_

#include <vector>
#include <type_traits>
#include <fairlogger/Logger.h>
#include <variant>

#include "rANS/internal/common/utils.h"
#include "rANS/internal/containers/RenormedHistogram.h"
#include "rANS/internal/containers/Symbol.h"
#include "rANS/internal/containers/ReverseSymbolLookupTable.h"
#include "rANS/internal/containers/SymbolTable.h"

namespace o2::rans
{

namespace internal
{

template <typename source_T>
class ReverseLookupDecoderTable
{
 public:
  using source_type = source_T;
  using value_type = internal::DecoderSymbol<source_type>;
  using count_type = uint32_t;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using container_type = std::vector<value_type>;
  using iterator_type = const value_type*;

  // TODO(milettri): fix once ROOT cling respects the standard http://wg21.link/p1286r2
  inline ReverseLookupDecoderTable() noexcept {}; // NOLINT

  explicit ReverseLookupDecoderTable(const RenormedHistogram<source_type>& renormedHistogram);

  [[nodiscard]] inline size_type size() const noexcept { return mLut.size(); };

  [[nodiscard]] inline bool isEscapeSymbol(count_type cumul) const noexcept { return cumul >= this->size(); };

  [[nodiscard]] inline bool hasEscapeSymbol() const noexcept { return this->mEscapeSymbol.getDecoderSymbol().getFrequency() > 0; };

  [[nodiscard]] inline const value_type& getEscapeSymbol() const noexcept { return this->mEscapeSymbol; };

  [[nodiscard]] inline const value_type& operator[](count_type cumul) const noexcept
  {
    assert(cumul < this->size());
    return this->mLut[cumul];
  };

  [[nodiscard]] inline size_type getPrecision() const noexcept { return this->mSymbolTablePrecision; };

  [[nodiscard]] inline const iterator_type begin() const noexcept { return this->mLut.data(); };

  [[nodiscard]] inline const iterator_type end() const noexcept { return this->mLut.data() + this->size(); };

 private:
  size_type mSymbolTablePrecision{};
  value_type mEscapeSymbol{};
  container_type mLut{};
};

template <typename source_T>
ReverseLookupDecoderTable<source_T>::ReverseLookupDecoderTable(const RenormedHistogram<source_type>& renormedHistogram) : mSymbolTablePrecision{renormedHistogram.getRenormingBits()}
{
  if (renormedHistogram.empty()) {
    LOG(warning) << "SymbolStatistics of empty message passed to " << __func__;
  }

  this->mLut.reserve(renormedHistogram.getNumSamples());
  auto histogramView = trim(makeHistogramView(renormedHistogram));

  this->mEscapeSymbol = [&]() -> value_type {
    const count_type symbolFrequency = renormedHistogram.getIncompressibleSymbolFrequency();
    const count_type cumulatedFrequency = renormedHistogram.getNumSamples() - symbolFrequency;
    return {0, symbolFrequency, cumulatedFrequency};
  }();

  source_type sourceSymbol = histogramView.getOffset();
  count_type cumulative = 0;
  for (count_type symbolFrequency : histogramView) {
    this->mLut.insert(mLut.end(), symbolFrequency, {sourceSymbol++, symbolFrequency, cumulative});
    cumulative += symbolFrequency;
  }
};

template <typename source_T>
class CombinedDecoderTable
{
 public:
  using source_type = source_T;
  using value_type = internal::DecoderSymbol<source_type>;
  using count_type = uint32_t;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using container_type = std::vector<value_type>;
  using iterator_type = const value_type*;

  // TODO(milettri): fix once ROOT cling respects the standard http://wg21.link/p1286r2
  inline CombinedDecoderTable() noexcept {}; // NOLINT

  explicit CombinedDecoderTable(const RenormedHistogram<source_type>& renormedHistogram) : mSymbolTable{renormedHistogram}, mRLUT{renormedHistogram}, mIncompressibleSymbol{0, mSymbolTable.getEscapeSymbol()} {};

  [[nodiscard]] inline size_type size() const noexcept { return mRLUT.size(); };

  [[nodiscard]] inline bool isEscapeSymbol(count_type cumul) const noexcept { return mRLUT.isIncompressible(cumul); };

  [[nodiscard]] inline bool hasEscapeSymbol() const noexcept { return mSymbolTable.hasEscapeSymbol(); };

  [[nodiscard]] inline const value_type& getEscapeSymbol() const noexcept { return mIncompressibleSymbol; };

  [[nodiscard]] inline const value_type operator[](count_type cumul) const noexcept
  {
    assert(cumul < this->size());
    source_type symbol = mRLUT[cumul];
    return {symbol, *mSymbolTable.lookupUnsafe(symbol)};
  };

  [[nodiscard]] inline size_type getPrecision() const noexcept { return this->mSymbolTable.getPrecision(); };

  [[nodiscard]] inline const iterator_type begin() const noexcept { return this->mLut.data(); };

  [[nodiscard]] inline const iterator_type end() const noexcept { return this->mLut.data() + this->size(); };

 private:
  SymbolTable<source_type, Symbol> mSymbolTable;
  ReverseSymbolLookupTable<source_type> mRLUT;
  value_type mIncompressibleSymbol{};
};

} // namespace internal

template <typename source_T>
class DecoderSymbolTable
{
 public:
  using source_type = source_T;
  using value_type = internal::DecoderSymbol<source_type>;
  using count_type = uint32_t;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

 private:
  using lut_type = internal::ReverseLookupDecoderTable<source_type>;
  using combined_type = internal::CombinedDecoderTable<source_type>;

  enum class impl_tag { lut,
                        combined };

 public:
  // TODO(milettri): fix once ROOT cling respects the standard http://wg21.link/p1286r2
  inline DecoderSymbolTable() noexcept {}; // NOLINT

  explicit DecoderSymbolTable(const RenormedHistogram<source_type>& renormedHistogram)
  {

    if (renormedHistogram.size() > utils::pow2(renormedHistogram.getRenormingBits())) {
      mImplTag = impl_tag::lut;
      // LOGP(info, "using LUT");
      mImpl = std::make_shared<lut_type>(renormedHistogram);
    } else {
      mImplTag = impl_tag::combined;
      // LOGP(info, "using Combined");
      mImpl = std::make_shared<combined_type>(renormedHistogram);
    }
  };

  [[nodiscard]] inline size_type size() const noexcept
  {
    return mImplTag == impl_tag::lut ? getLUT().size() : getCombined().size();
  };

  [[nodiscard]] inline bool isEscapeSymbol(count_type cumul) const noexcept
  {
    return mImplTag == impl_tag::lut ? getLUT().isEscapeSymbol(cumul) : getCombined().isEscapeSymbol(cumul);
  };

  [[nodiscard]] inline bool hasEscapeSymbol() const noexcept
  {
    return mImplTag == impl_tag::lut ? getLUT().hasEscapeSymbol() : getCombined().hasEscapeSymbol();
  };

  [[nodiscard]] inline const value_type& getEscapeSymbol() const noexcept
  {
    return mImplTag == impl_tag::lut ? getLUT().getEscapeSymbol() : getCombined().getEscapeSymbol();
  };

  [[nodiscard]] inline const value_type operator[](count_type cumul) const noexcept
  {
    return mImplTag == impl_tag::lut ? getLUT()[cumul] : getCombined()[cumul];
  };

  [[nodiscard]] inline size_type getPrecision() const noexcept
  {
    return mImplTag == impl_tag::lut ? getLUT().getPrecision() : getCombined().getPrecision();
  };

 private:
  inline lut_type& getLUT() noexcept { return *reinterpret_cast<lut_type*>(mImpl.get()); };
  inline combined_type& getCombined() noexcept { return *reinterpret_cast<combined_type*>(mImpl.get()); };
  inline const lut_type& getLUT() const noexcept { return *reinterpret_cast<const lut_type*>(mImpl.get()); }
  inline const combined_type& getCombined() const noexcept { return *reinterpret_cast<const combined_type*>(mImpl.get()); }

  impl_tag mImplTag{};
  std::shared_ptr<void> mImpl{};
};

} // namespace o2::rans

#endif /* RANS_INTERNAL_CONTAINERS_DECODERSYMBOLTABLE_H_ */
