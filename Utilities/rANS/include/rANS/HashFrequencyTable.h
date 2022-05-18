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

#ifndef INCLUDE_RANS_HASHFREQUENCYTABLE_H_
#define INCLUDE_RANS_HASHFREQUENCYTABLE_H_

#include "rANS/definitions.h"
#include "rANS/internal/helper.h"
#include "rANS/internal/FrequencyTableBase.h"
#include "rANS/internal/HashFrequencyContainer.h"

namespace o2
{
namespace rans
{

template <typename source_T>
class HashFrequencyTable : public HashFrequencyContainer<source_T>,
                           public FrequencyTableBase<source_T,
                                                     typename HashFrequencyContainer<source_T>::value_type,
                                                     HashFrequencyTable<source_T>>
{
  using containerBase_type = HashFrequencyContainer<source_T>;
  using frequencyTableBase_type = FrequencyTableBase<source_T, typename HashFrequencyContainer<source_T>::value_type, HashFrequencyTable<source_T>>;

 public:
  using source_type = typename containerBase_type::source_type;
  using index_type = typename containerBase_type::index_type;
  using value_type = typename containerBase_type::value_type;
  using container_type = typename containerBase_type::container_type;
  using size_type = typename containerBase_type::size_type;
  using difference_type = typename containerBase_type::difference_type;
  using reference = typename containerBase_type::reference;
  using const_reference = typename containerBase_type::const_reference;
  using pointer = typename containerBase_type::pointer;
  using const_pointer = typename containerBase_type::const_pointer;
  using const_iterator = typename containerBase_type::const_iterator;

  HashFrequencyTable() = default;

  template <typename freq_IT>
  HashFrequencyTable(freq_IT begin, freq_IT end, source_type offset) : containerBase_type(), frequencyTableBase_type{begin, end, offset} {};

  // operations
  template <typename source_IT>
  HashFrequencyTable& addSamples(source_IT begin, source_IT end);

  using frequencyTableBase_type::addSamples;

  template <typename freq_IT>
  HashFrequencyTable& addFrequencies(freq_IT begin, freq_IT end, source_type offset);

  using frequencyTableBase_type::addFrequencies;
};

template <typename source_T>
template <typename source_IT>
auto HashFrequencyTable<source_T>::addSamples(source_IT begin, source_IT end) -> HashFrequencyTable&
{
  std::for_each(begin, end, [this](const source_type& symbol) { 
      ++this->mNSamples;
      ++this->mContainer[static_cast<index_type>(symbol)]; });
  return *this;
}

template <typename source_T>
template <typename freq_IT>
auto HashFrequencyTable<source_T>::addFrequencies(freq_IT begin, freq_IT end, source_type offset) -> HashFrequencyTable&
{
  source_type sourceSymbol = offset;
  for (auto iter = begin; iter != end; ++iter) {
    auto value = *iter;
    if (value > 0) {
      this->mContainer[sourceSymbol] += value;
      this->mNSamples += value;
    }
    ++sourceSymbol;
  }
  return *this;
}

} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_HASHFREQUENCYTABLE_H_ */
