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

/// @file   Symbol.h
/// @author Michael Lettrich
/// @since  2019-05-21
/// @brief  Structure containing all relevant information for decoding a rANS encoded symbol

#ifndef RANS_INTERNAL_SYMBOL_H
#define RANS_INTERNAL_SYMBOL_H

#include <cassert>
#include <cstdint>
#include <cstring>

#include <fairlogger/Logger.h>

#include "rANS/definitions.h"
#include "rANS/internal/helper.h"

namespace o2
{
namespace rans
{
namespace internal
{
// Decoder symbols are straightforward.
class Symbol
{
 public:
  using value_type = count_t;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;

  //TODO(milettri): fix once ROOT cling respects the standard http://wg21.link/p1286r2
  constexpr Symbol() noexcept {}; //NOLINT
  constexpr Symbol(value_type frequency, value_type cumulative, size_type symbolTablePrecision)
    : mFrequency(frequency), mCumulative(cumulative)
  {
    (void)symbolTablePrecision; // silence compiler warnings if assert not compiled.
    assert(mCumulative <= pow2(symbolTablePrecision));
    assert(mFrequency <= pow2(symbolTablePrecision) - mCumulative);
  };
  [[nodiscard]] inline constexpr value_type getFrequency() const noexcept { return mFrequency; };
  [[nodiscard]] inline constexpr value_type getCumulative() const noexcept { return mCumulative; };

  [[nodiscard]] inline bool operator==(const Symbol& other) const { return this->mCumulative == other.mCumulative; };

  friend std::ostream& operator<<(std::ostream& os, const Symbol& symbol)
  {
    os << fmt::format("Symbol:{{Frequency: {}, Cumulative: {}}}", symbol.getFrequency(), symbol.getCumulative());
    return os;
  }

 protected:
  value_type mFrequency{};
  value_type mCumulative{};
};

class PrecomputedSymbol : public Symbol
{
  using base_type = Symbol;

 public:
  using value_type = typename base_type::value_type;
  using state_type = uint64_t;
  using size_type = typename base_type::size_type;
  using difference_type = typename base_type::difference_type;

  //TODO(milettri): fix once ROOT cling respects the standard http://wg21.link/p1286r2
  constexpr PrecomputedSymbol() noexcept {}; //NOLINT

  constexpr PrecomputedSymbol(value_type frequency, value_type cumulative, size_t symbolTablePrecision) : base_type{}
  {
    assert(cumulative <= pow2(symbolTablePrecision));
    assert(frequency <= pow2(symbolTablePrecision) - cumulative);

    // Say M := 1 << symbolTablePrecision.
    //
    // The original encoder does:
    //   x_new = (x/frequency)*M + cumulative + (x%frequency)
    //
    // The fast encoder does (schematically):
    //   q     = mul_hi(x, rcp_freq) >> rcp_shift   (division)
    //   r     = x - q*frequency                         (remainder)
    //   x_new = q*M + bias + r                     (new x)
    // plugging in r into x_new yields:
    //   x_new = bias + x + q*(M - frequency)
    //        =: bias + x + q*cmpl_freq             (*)
    //
    // and we can just precompute cmpl_freq. Now we just need to
    // set up our parameters such that the original encoder and
    // the fast encoder agree.

    mFrequency = frequency;
    mFrequencyComplement = static_cast<state_type>((pow2(symbolTablePrecision)) - frequency);
    if (frequency < 2) {
      // frequency=0 symbols are never valid to encode, so it doesn't matter what
      // we set our values to.
      //
      // frequency=1 is tricky, since the reciprocal of 1 is 1; unfortunately,
      // our fixed-point reciprocal approximation can only multiply by values
      // smaller than 1.
      //
      // So we use the "next best thing": rcp_freq=0xffffffff, rcp_shift=0.
      // This gives:
      //   q = mul_hi(x, rcp_freq) >> rcp_shift
      //     = mul_hi(x, (1<<32) - 1)) >> 0
      //     = floor(x - x/(2^32))
      //     = x - 1 if 1 <= x < 2^32
      // and we know that x>0 (x=0 is never in a valid normalization interval).
      //
      // So we now need to choose the other parameters such that
      //   x_new = x*M + cumulative
      // plug it in:
      //     x*M + cumulative                   (desired result)
      //   = bias + x + q*cmpl_freq        (*)
      //   = bias + x + (x - 1)*(M - 1)    (plug in q=x-1, cmpl_freq)
      //   = bias + 1 + (x - 1)*M
      //   = x*M + (bias + 1 - M)
      //
      // so we have cumulative = bias + 1 - M, or equivalently
      //   bias = cumulative + M - 1.
      mReciprocalFrequency = static_cast<state_type>(~0ul);
      mReciprocalShift = 0;
      mCumulative = cumulative + (pow2(symbolTablePrecision)) - 1;
    } else {
      // Alverson, "Integer Division using reciprocals"
      const uint32_t shift = std::ceil(std::log2(frequency));

      // long divide ((uint128) (1 << (shift + 63)) + frequency-1) / frequency
      // by splitting it into two 64:64 bit divides (this works because
      // the dividend has a simple form.)
      uint64_t x0 = frequency - 1;
      const uint64_t x1 = 1ull << (shift + 31);

      const uint64_t t1 = x1 / frequency;
      x0 += (x1 % frequency) << 32;
      const uint64_t t0 = x0 / frequency;

      mReciprocalFrequency = t0 + (t1 << 32);

      mReciprocalShift = shift - 1;

      // With these values, 'q' is the correct quotient, so we
      // have bias=cumulative.
      mCumulative = cumulative;
    }
  };

  inline constexpr state_type getReciprocalFrequency() const noexcept { return mReciprocalFrequency; };
  inline constexpr value_type getFrequencyComplement() const noexcept { return mFrequencyComplement; };
  inline constexpr value_type getReciprocalShift() const noexcept { return mReciprocalShift; };

  friend std::ostream& operator<<(std::ostream& os, const PrecomputedSymbol& symbol)
  {
    os << fmt::format("PrecomputedSymbol{{Frequency: {},Cumulative: {}, ReciprocalFrequency {}, FrequencyComplement {}, mReciprocalShift {}}}",
                      symbol.getFrequency(),
                      symbol.getCumulative(),
                      symbol.getReciprocalFrequency(),
                      symbol.getFrequencyComplement(),
                      symbol.getReciprocalShift());
    return os;
  };

 private:
  state_type mReciprocalFrequency{}; // Fixed-point reciprocal frequency
  value_type mFrequencyComplement{}; // Complement of frequency: (1 << symbolTablePrecision) - frequency
  value_type mReciprocalShift{};     // Reciprocal shift
};

} // namespace internal
} // namespace rans
} // namespace o2

#endif /* RANS_INTERNAL_SYMBOL_H */
