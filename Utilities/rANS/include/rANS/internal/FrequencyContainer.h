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

/// @file   FrequencyContainer.h
/// @author Michael Lettrich
/// @since  2019-05-08
/// @brief Histogram to depict frequencies of source symbols for rANS compression.

#ifndef INCLUDE_RANS_FREQUENCYCONTAINER_H_
#define INCLUDE_RANS_FREQUENCYCONTAINER_H_

#include <cstdint>
#include <string>

namespace o2
{
namespace rans
{

template <class source_T, class index_T, class value_T, class container_T, class derived_T>
class FrequencyContainer
{
 public:
  using source_type = source_T;
  using index_type = index_T;
  using value_type = value_T;
  using container_type = container_T;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;
  using const_iterator = typename container_type::const_iterator;

  // accessors
  [[nodiscard]] inline value_type operator[](source_type sourceSymbol) const { static_cast<derived_T*>(this)->operator[](sourceSymbol); };

  [[nodiscard]] inline const_iterator cbegin() const noexcept { return mContainer.begin(); };

  [[nodiscard]] inline const_iterator cend() const noexcept { return mContainer.end(); };

  [[nodiscard]] inline const_iterator begin() const noexcept { return cbegin(); };

  [[nodiscard]] inline const_iterator end() const noexcept { return cend(); };

  [[nodiscard]] inline size_type size() const noexcept { return static_cast<derived_T*>(this)->size(); };

  [[nodiscard]] inline bool empty() const noexcept { return mNSamples == 0; };

  [[nodiscard]] inline size_type computeNUsedAlphabetSymbols() const noexcept { return static_cast<derived_T*>(this)->computeNUsedAlphabetSymbols(); };

  [[nodiscard]] inline size_type getNumSamples() const noexcept { return mNSamples; };

  [[nodiscard]] inline source_type getOffset() const noexcept { return mOffset; };

  [[nodiscard]] inline container_type release() &&
  {
    container_type t{};
    std::swap(t, mContainer);
    mOffset = 0;
    mNSamples = 0;
    return t;
  };

 protected:
  FrequencyContainer() = default;

  container_type mContainer{};
  source_type mOffset{};
  size_type mNSamples{};
};
} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_FREQUENCYCONTAINER_H_ */
