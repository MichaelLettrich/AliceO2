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

#ifndef INCLUDE_RANS_HASHFREQUENCYCONTAINER_H_
#define INCLUDE_RANS_HASHFREQUENCYCONTAINER_H_

#include <absl/container/flat_hash_map.h>

#include "rANS/definitions.h"
#include "rANS/internal/FrequencyContainer.h"

namespace o2
{
namespace rans
{

template <typename source_T>
class HashFrequencyContainer : public FrequencyContainer<source_T,
                                                         source_T,
                                                         count_t,
                                                         absl::flat_hash_map<source_T, count_t>,
                                                         HashFrequencyContainer<source_T>>
{
  using base_type = FrequencyContainer<source_T, source_T, count_t, absl::flat_hash_map<source_T, count_t>, HashFrequencyContainer<source_T>>;

 public:
  using source_type = source_T;
  using index_type = source_type;
  using value_type = count_t;
  using container_type = absl::flat_hash_map<index_type, value_type>;
  using size_type = typename base_type::size_type;
  using difference_type = typename base_type::difference_type;
  using reference = typename base_type::reference;
  using const_reference = typename base_type::const_reference;
  using pointer = typename base_type::pointer;
  using const_pointer = typename base_type::const_pointer;
  using const_iterator = typename base_type::const_iterator;

  // constructors and operators are implicit.

  // accessors
  [[nodiscard]] inline value_type operator[](source_type sourceSymbol) const
  {
    auto iter = this->mContainer.find(sourceSymbol);
    if (iter == this->mContainer.end()) {
      return 0;
    } else {
      return iter->second;
    }
  };

  [[nodiscard]] inline constexpr size_type size() const noexcept { return this->mContainer.size(); };

  [[nodiscard]] inline size_type computeNUsedAlphabetSymbols() const noexcept { return this->size(); };

 protected:
  HashFrequencyContainer() = default;
};

} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_HASHFREQUENCYCONTAINER_H_ */
