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

#ifndef INCLUDE_RANS_FREQUENCYTABLEBASE_H_
#define INCLUDE_RANS_FREQUENCYTABLEBASE_H_

#include <gsl/span>

#include "rANS/internal/helper.h"

namespace o2
{
namespace rans
{
template <class source_T, class value_T, class derived_T>
class FrequencyTableBase
{

 public:
  using source_type = source_T;
  using value_type = value_T;

  // operations
  template <typename source_IT>
  inline derived_T& addSamples(source_IT begin, source_IT end)
  {
    static_assert(internal::isCompatibleIter_v<source_type, source_IT>);
    return static_cast<derived_T*>(this)->addSamples(begin, end);
  };

  inline derived_T& addSamples(gsl::span<const source_type> samples)
  {
    return addSamples(samples.data(), samples.data() + samples.size());
  };

  template <typename freq_IT>
  inline derived_T& addFrequencies(freq_IT begin, freq_IT end, source_type offset)
  {
    static_assert(internal::isCompatibleIter_v<value_type, freq_IT>);
    return static_cast<derived_T*>(this)->addFrequencies(begin, end, offset);
  };

  inline derived_T& addFrequencies(gsl::span<const value_type> frequencies, source_type offset)
  {
    return addFrequencies(frequencies.data(), frequencies.data() + frequencies.size(), offset);
  };

  derived_T& operator+(derived_T& other)
  {
    return addFrequencies(other.cbegin(), other.cbegin(), 0);
  };

 protected:
  FrequencyTableBase() = default;

  template <typename freq_IT>
  FrequencyTableBase(freq_IT begin, freq_IT end, source_type offset)
  {
    static_assert(internal::isIntegralIter_v<freq_IT>);
    addFrequencies(begin, end, offset);
  };
};

} // namespace rans
} // namespace o2

#endif /* INCLUDE_RANS_FREQUENCYTABLEBASE_H_ */
