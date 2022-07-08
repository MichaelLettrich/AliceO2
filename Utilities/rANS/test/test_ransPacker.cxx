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
#include <fmt/core.h>

#include "rANS/Packer.h"

using source_type = uint64_t;

std::vector<source_type> generate(source_type minNumber, source_type maxNumber)
{
  std::vector<source_type> ret;

  std::mt19937 mt(0); // same seed we want always the same distrubution of random numbers;
  std::uniform_int_distribution<source_type> dist(minNumber, maxNumber);
  ret.resize(257);
  std::generate(ret.begin(), ret.end(), [&dist, &mt]() { return dist(mt); });
  return ret;
};

BOOST_AUTO_TEST_CASE(test_pack)
{

  for (size_t packingWidth = 1; packingWidth <= 63; ++packingWidth) {
    BOOST_TEST_MESSAGE(fmt::format("Checking {} Bit Packing", packingWidth));

    auto source = generate(0, (1ull << packingWidth) - 1);
    std::vector<uint8_t> packingBuffer(source.size() * sizeof(source_type), 0);

    uint64_t* iter = reinterpret_cast<uint64_t*>(packingBuffer.data());
    uint32_t bitOffset = 0;
    for (auto i : source) {
      o2::rans::pack(iter, bitOffset, i, packingWidth);
    };

    std::vector<uint8_t> packingBuffer2(source.size() * sizeof(source_type), 0);

    // use buffer API
    auto packingBuffer2End = reinterpret_cast<uint8_t*>(o2::rans::packStream(source.begin(),
                                                                             source.end(),
                                                                             reinterpret_cast<uint64_t*>(packingBuffer2.data()),
                                                                             packingWidth));

    // check if results are equal;
    BOOST_CHECK_EQUAL_COLLECTIONS(packingBuffer.data(), reinterpret_cast<uint8_t*>(++iter), packingBuffer2.data(), packingBuffer2End);
  }
};

BOOST_AUTO_TEST_CASE(test_packUnpackNext)
{

  for (size_t packingWidth = 1; packingWidth <= 63; ++packingWidth) {
    BOOST_TEST_MESSAGE(fmt::format("Checking {} Bit Packing", packingWidth));
    auto source = generate(0, (1ull << packingWidth) - 1);
    std::vector<uint8_t> packingBuffer(source.size() * sizeof(source_type), 0);

    uint64_t* iter = reinterpret_cast<uint64_t*>(packingBuffer.data());
    uint32_t bitOffset = 0;
    for (auto i : source) {
      o2::rans::pack(iter, bitOffset, i, packingWidth);
    };

    //unpack next
    std::vector<source_type> unpackBuffer;
    const uint64_t* unpackIter = reinterpret_cast<uint64_t*>(packingBuffer.data());
    bitOffset = 0;
    for (size_t i = 0; i < source.size(); ++i) {
      unpackBuffer.push_back(o2::rans::unpackNext(unpackIter, bitOffset, packingWidth));
    };

    // compare if both yield the correct result
    BOOST_CHECK_EQUAL_COLLECTIONS(source.begin(), source.end(), unpackBuffer.begin(), unpackBuffer.end());
  }
};

BOOST_AUTO_TEST_CASE(test_packRUnpackNext)
{

  for (size_t packingWidth = 1; packingWidth <= 63; ++packingWidth) {
    BOOST_TEST_MESSAGE(fmt::format("Checking {} Bit Packing", packingWidth));
    auto source = generate(0, (1ull << packingWidth) - 1);
    std::vector<uint8_t> packingBuffer(source.size() * sizeof(source_type), 0);

    uint64_t* iter = reinterpret_cast<uint64_t*>(packingBuffer.data());
    uint32_t bitOffset = 0;
    for (auto i : source) {
      o2::rans::pack(iter, bitOffset, i, packingWidth);
    };

    //unpack next
    std::vector<source_type> unpackBuffer;
    for (size_t i = 0; i < source.size(); ++i) {
      unpackBuffer.push_back(o2::rans::rUnpackNext(iter, bitOffset, packingWidth));
    };

    // compare if both yield the correct result
    BOOST_CHECK_EQUAL_COLLECTIONS(source.begin(), source.end(), unpackBuffer.rbegin(), unpackBuffer.rend());
  }
};

BOOST_AUTO_TEST_CASE(test_packUnpackByIdx)
{

  for (size_t packingWidth = 1; packingWidth <= 63; ++packingWidth) {
    BOOST_TEST_MESSAGE(fmt::format("Checking {} Bit Packing", packingWidth));
    auto source = generate(0, (1ull << packingWidth) - 1);
    std::vector<uint8_t> packingBuffer(source.size() * sizeof(source_type), 0);

    uint64_t* iter = reinterpret_cast<uint64_t*>(packingBuffer.data());
    uint32_t bitOffset = 0;
    for (auto i : source) {
      o2::rans::pack(iter, bitOffset, i, packingWidth);
    };

    //unpack using indexed access
    std::vector<source_type> unpackBuffer;
    iter = reinterpret_cast<uint64_t*>(packingBuffer.data());
    for (size_t i = 0; i < source.size(); ++i) {
      unpackBuffer.push_back(o2::rans::unpackByIndex(iter, i, packingWidth));
    };

    // compare if both yield the correct result
    BOOST_CHECK_EQUAL_COLLECTIONS(source.begin(), source.end(), unpackBuffer.begin(), unpackBuffer.end());
  }
};

BOOST_AUTO_TEST_CASE(test_packUnpackStream)
{

  for (size_t packingWidth = 1; packingWidth <= 63; ++packingWidth) {
    BOOST_TEST_MESSAGE(fmt::format("Checking {} Bit Packing", packingWidth));
    auto source = generate(0, (1ull << packingWidth) - 1);
    std::vector<uint8_t> packingBuffer(source.size() * sizeof(source_type), 0);

    // use buffer API
    auto packingBufferEnd = reinterpret_cast<uint8_t*>(o2::rans::packStream(source.begin(),
                                                                            source.end(),
                                                                            reinterpret_cast<uint64_t*>(packingBuffer.data()),
                                                                            packingWidth));

    std::vector<source_type> unpackBuffer(source.size(), 0);

    // unpack full buffer
    o2::rans::unpackStream(reinterpret_cast<uint64_t*>(packingBuffer.data()),
                           unpackBuffer.data(),
                           source.size(),
                           packingWidth);

    // compare if both yield the correct result
    BOOST_CHECK_EQUAL_COLLECTIONS(source.begin(), source.end(), unpackBuffer.begin(), unpackBuffer.end());
  }
};

BOOST_AUTO_TEST_CASE(test_packRUnpackEliasDelta)
{

  for (size_t packingWidth = 1; packingWidth <= 32; ++packingWidth) {
    BOOST_TEST_MESSAGE(fmt::format("Checking {} Bit Elias Delta Coder", packingWidth));
    auto source = generate(1ull, (1ull << packingWidth) - 1);
    std::vector<uint8_t> packingBuffer(source.size() * sizeof(source_type), 0);

    uint64_t* iter = reinterpret_cast<uint64_t*>(packingBuffer.data());
    uint32_t bitOffset = 0;
    for (auto i : source) {
      o2::rans::eliasDeltaEncode(iter, bitOffset, i);
    };

    std::vector<uint32_t> unpackBuffer;
    for (size_t i = 0; i < source.size(); ++i) {
      unpackBuffer.push_back(o2::rans::eliasDeltaDecode(iter, bitOffset));
    };

    // compare if both yield the correct result
    BOOST_CHECK_EQUAL_COLLECTIONS(source.begin(), source.end(), unpackBuffer.rbegin(), unpackBuffer.rend());
  }
};