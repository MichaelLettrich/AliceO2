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

/// @file   helpers.h
/// @author Michael Lettrich
/// @brief  common functionality for rANS benchmarks.

#ifndef RANS_BENCHMARKS_UTILS_BUFFERS_H_
#define RANS_BENCHMARKS_UTILS_BUFFERS_H_

#include <algorithm>
#include <vector>
#include <concepts>
#include <cstdint>

#ifdef RANS_PARALLEL_STL
#include <execution>
#endif

namespace o2::rans::benchmarks
{

template <std::integral source_T, std::integral stream_T = uint32_t>
struct EncodeBuffer {

  EncodeBuffer() = default;
  EncodeBuffer(size_t sourceSize) : buffer(2 * sourceSize), literals()
  {
    literals.reserve(sourceSize + 32);
    literalsEnd = literals.data();
    encodeBufferEnd = buffer.data();
  };

  std::vector<stream_T> buffer{};
  std::vector<source_T> literals{};
  stream_T* encodeBufferEnd{buffer.data()};
  source_T* literalsEnd{literals.data()};
};

template <std::integral source_T>
struct DecodeBuffer {

  DecodeBuffer() = default;
  DecodeBuffer(size_t sourceSize) : buffer(sourceSize) {};

  template <typename T>
  bool operator==(const T& correct)
  {
#ifdef RANS_PARALLEL_STL
    return std::equal(std::execution::par_unseq, buffer.begin(), buffer.end(), std::begin(correct), std::end(correct));
#else
    return std::equal(buffer.begin(), buffer.end(), std::begin(correct), std::end(correct));
#endif // RANS_PARALLEL_STL
  }

  std::vector<source_T> buffer{};
};

} // namespace o2::rans::benchmarks

#endif /* RANS_BENCHMARKS_UTILS_BUFFERS_H_ */