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

/// @file   DecoderSymbol.h
/// @author Michael Lettrich
/// @since  2020-04-15
/// @brief  Test rANS encoder/ decoder

#define BOOST_TEST_MODULE Utility test
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <vector>
#include <cstring>
#include <random>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/vector.hpp>

#include "rANS/Packer.h"

using source_type = uint32_t;

std::vector<source_type> generate(source_type minNumber, source_type maxNumber)
{
  std::vector<source_type> ret;

  std::mt19937 mt(0); // same seed we want always the same distrubution of random numbers;
  std::uniform_int_distribution<source_type> dist(minNumber, maxNumber);
  ret.resize(257);
  std::generate(ret.begin(), ret.end(), [&dist, &mt]() { return dist(mt); });
  return ret;
};

BOOST_AUTO_TEST_CASE(test_packUnpack)
{
  constexpr size_t packWidth = 17;

  auto source = generate(0, (1u << packWidth) - 1);
  std::vector<uint8_t> packingBuffer(source.size() * sizeof(source_type), 0);
  std::vector<uint8_t> packingBuffer2(source.size() * sizeof(source_type), 0);

  // use element wise API
  uint64_t* iter = reinterpret_cast<uint64_t*>(packingBuffer.data());
  size_t pos = 0;
  for (auto i : source) {
    o2::rans::pack(iter, pos, i, packWidth);
  };

  // use buffer API
  auto packingBuffer2End = reinterpret_cast<uint8_t*>(o2::rans::packStream(source.begin(),
                                                                           source.end(),
                                                                           reinterpret_cast<uint64_t*>(packingBuffer2.data()),
                                                                           packWidth));

  // check if results are equal;
  BOOST_CHECK_EQUAL_COLLECTIONS(packingBuffer.data(), reinterpret_cast<uint8_t*>(++iter), packingBuffer2.data(), packingBuffer2End);

  //unpack using indexed access
  std::vector<source_type> unpackBuffer;
  iter = reinterpret_cast<uint64_t*>(packingBuffer.data());
  for (size_t i = 0; i < source.size(); ++i) {
    unpackBuffer.push_back(o2::rans::unpackByIndex(iter, i, packWidth));
  };

  std::vector<source_type> unpackBuffer2(source.size(), 0);

  // unpack full buffer
  o2::rans::unpackStream(reinterpret_cast<uint64_t*>(packingBuffer.data()),
                         unpackBuffer2.data(),
                         source.size(),
                         packWidth);

  // compare if both yield the correct result
  BOOST_CHECK_EQUAL_COLLECTIONS(source.begin(), source.end(), unpackBuffer.begin(), unpackBuffer.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(source.begin(), source.end(), unpackBuffer2.begin(), unpackBuffer2.end());
};