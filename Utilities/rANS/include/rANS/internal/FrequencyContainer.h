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

#include "rANS/internal/ContainerInterface.h"

namespace o2
{
namespace rans
{
namespace internal
{

template <class source_T, class index_T, class value_T, class container_T, class derived_T>
class FrequencyContainer : public ContainerInterface<source_T, index_T, value_T, container_T,
                                                     FrequencyContainer<source_T, index_T, value_T, container_T, derived_T>>
{
  using base_type = ContainerInterface<source_T, index_T, value_T, container_T,
                                       FrequencyContainer<source_T, index_T, value_T, container_T, derived_T>>;

 public:
  using source_type = typename base_type::source_type;
  using index_type = typename base_type::index_type;
  using value_type = typename base_type::value_type;
  using container_type = typename base_type::container_type;
  using size_type = typename base_type::size_type;
  using difference_type = typename base_type::difference_type;
  using reference = typename base_type::reference;
  using const_reference = typename base_type::const_reference;
  using pointer = typename base_type::pointer;
  using const_pointer = typename base_type::const_pointer;
  using const_iterator = typename base_type::const_iterator;

  // accessors
  [[nodiscard]] inline value_type operator[](source_type sourceSymbol) const { static_cast<const derived_T*>(this)->operator[](sourceSymbol); };

  [[nodiscard]] inline size_type size() const noexcept { return static_cast<const derived_T*>(this)->size(); };

  [[nodiscard]] inline bool empty() const noexcept { return mNSamples == 0; };

  [[nodiscard]] inline size_type computeNUsedAlphabetSymbols() const noexcept { return static_cast<const derived_T*>(this)->computeNUsedAlphabetSymbols(); };

  [[nodiscard]] inline size_type getNumSamples() const noexcept { return mNSamples; };

  friend void swap(FrequencyContainer& a, FrequencyContainer& b) noexcept
  {
    using std::swap;
    swap(a.mNSamples, b.mNSamples);
    swap(static_cast<typename FrequencyContainer::base_type&>(a),
         static_cast<typename FrequencyContainer::base_type&>(b));
  };

 protected:
  FrequencyContainer() = default;

  size_type mNSamples{};
};

} // namespace internal
} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_FREQUENCYCONTAINER_H_ */
