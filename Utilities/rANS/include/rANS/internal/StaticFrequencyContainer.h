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

/// @file   FrequencyTableImpl.h
/// @author Michael Lettrich
/// @since  2019-05-08
/// @brief Histogram to depict frequencies of source symbols for rANS compression.

#ifndef INCLUDE_RANS_STATICFREQUENCYCONTAINER_H_
#define INCLUDE_RANS_STATICFREQUENCYCONTAINER_H_

#include <numeric>
#include <vector>
#include <algorithm>

#include "rANS/definitions.h"
#include "rANS/internal/helper.h"
#include "rANS/internal/FrequencyContainer.h"
#include "rANS/internal/StaticContainerIterator.h"

namespace o2
{
namespace rans
{
namespace internal
{

template <typename source_T>
class StaticFrequencyContainer : public FrequencyContainer<
                                   source_T,
                                   std::make_unsigned_t<source_T>,
                                   count_t,
                                   std::vector<count_t>,
                                   StaticContainerIterator<const StaticFrequencyContainer<source_T>>, StaticFrequencyContainer<source_T>>
{

  using base_type = FrequencyContainer<source_T, std::make_unsigned_t<source_T>, count_t, std::vector<count_t>, StaticContainerIterator<const StaticFrequencyContainer<source_T>>, StaticFrequencyContainer<source_T>>;

 public:
  using source_type = source_T;
  using index_type = std::make_unsigned_t<source_type>;
  using value_type = count_t;
  using container_type = std::vector<value_type>;
  using size_type = typename base_type::size_type;
  using difference_type = typename base_type::difference_type;
  using reference = typename base_type::reference;
  using const_reference = typename base_type::const_reference;
  using pointer = typename base_type::pointer;
  using const_pointer = typename base_type::const_pointer;
  using const_iterator = typename base_type::const_iterator;

  static_assert(sizeof(index_type) <= 2, "This datatype requires a <=16Bit datatype for source_T");

  // accessors
  [[nodiscard]] inline const_iterator cbegin() const noexcept { return StaticContainerIterator{this, this->getOffset()}; };

  [[nodiscard]] inline const_iterator cend() const noexcept { return StaticContainerIterator{this, this->getOffset() + static_cast<int64_t>(this->size())}; };

  [[nodiscard]] inline value_type operator[](source_type sourceSymbol) const { return this->mContainer[static_cast<index_type>(sourceSymbol)]; };

  [[nodiscard]] inline const_pointer data() const noexcept { return this->mContainer.data(); };

  [[nodiscard]] inline constexpr size_type size() const noexcept { return internal::pow2(internal::toBits(sizeof(index_type))); };

  [[nodiscard]] inline size_type computeNUsedAlphabetSymbols() const noexcept
  {
    return std::count_if(this->begin(), this->end(), [](value_type v) { return v != static_cast<value_type>(0); });
  };

  friend void swap(StaticFrequencyContainer& a, StaticFrequencyContainer& b) noexcept
  {
    using std::swap;
    swap(static_cast<typename StaticFrequencyContainer::base_type&>(a),
         static_cast<typename StaticFrequencyContainer::base_type&>(b));
  };

 protected:
  StaticFrequencyContainer()
  {
    this->mContainer.resize(this->size(), 0);
    this->mOffset = std::numeric_limits<source_type>::min();
  };
}; // namespace internal

} // namespace internal
} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_STATICFREQUENCYCONTAINER_H_ */