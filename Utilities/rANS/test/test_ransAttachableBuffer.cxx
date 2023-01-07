// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   test_ransAttachableBuffer.h
/// @author Michael Lettrich
/// @brief  Test AttachableBuffer class

#define BOOST_TEST_MODULE Utility test
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#undef NDEBUG
#include <cassert>

#include <vector>
#include <type_traits>
#include <array>

#include <boost/test/unit_test.hpp>

#include <fairlogger/Logger.h>

#include "rANS/internal/containers/AttachableBuffer.h"

using namespace o2::rans::internal;

BOOST_AUTO_TEST_CASE(test_emptyAttachableBuffer)
{
  AttachableBuffer<uint32_t> buffer{};
  BOOST_CHECK_EQUAL(buffer.size(), 0);
  BOOST_CHECK_EQUAL(buffer.data(), nullptr);
  BOOST_CHECK_EQUAL(buffer.empty(), true);
  BOOST_CHECK_EQUAL(buffer.begin(), buffer.end());
  BOOST_CHECK_EQUAL(buffer.cbegin(), buffer.cend());
  BOOST_CHECK((buffer.rbegin() == buffer.rend()));
  BOOST_CHECK((buffer.crbegin() == buffer.crend()));
  auto s = std::move(buffer).release();
  BOOST_CHECK_EQUAL(s.get(), nullptr);

  BOOST_CHECK_EQUAL(buffer.size(), 0);
  BOOST_CHECK_EQUAL(buffer.data(), nullptr);
  BOOST_CHECK_EQUAL(buffer.empty(), true);
}

BOOST_AUTO_TEST_CASE(test_constructAttachableBuffer)
{
  std::array a = {0u, 1u, 2u, 3u, 4u};

  AttachableBuffer<uint32_t> b1{a.data(), a.size()};

  BOOST_CHECK_EQUAL(b1.size(), a.size());
  BOOST_CHECK_EQUAL(b1.data(), a.data());
  BOOST_CHECK_EQUAL(b1.empty(), false);
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.cbegin(), b1.cend(), a.begin(), a.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.begin(), b1.end(), a.begin(), a.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.rbegin(), b1.rend(), a.rbegin(), a.rend());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.crbegin(), b1.crend(), a.crbegin(), a.crend());

  for (size_t i = 0; i < a.size(); i++) {
    BOOST_CHECK_EQUAL(b1[i], a[i]);
  }

  // copy construct
  auto b2{b1};
  BOOST_CHECK_EQUAL(b1.size(), b2.size());
  BOOST_CHECK_EQUAL(b1.data(), b2.data());
  BOOST_CHECK_EQUAL(b2.empty(), false);
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.cbegin(), b1.cend(), b2.begin(), b2.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.begin(), b1.end(), b2.begin(), b2.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.rbegin(), b1.rend(), b2.rbegin(), b2.rend());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.crbegin(), b1.crend(), b2.crbegin(), b2.crend());

  // copy assign
  auto b3 = b1;
  BOOST_CHECK_EQUAL(b1.size(), b3.size());
  BOOST_CHECK_EQUAL(b1.data(), b3.data());
  BOOST_CHECK_EQUAL(b3.empty(), false);
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.cbegin(), b1.cend(), b3.begin(), b3.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.begin(), b1.end(), b3.begin(), b3.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.rbegin(), b1.rend(), b3.rbegin(), b3.rend());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.crbegin(), b1.crend(), b3.crbegin(), b3.crend());

  // move construct
  auto b4{std::move(b2)};
  BOOST_CHECK_EQUAL(b1.size(), b4.size());
  BOOST_CHECK_EQUAL(b1.data(), b4.data());
  BOOST_CHECK_EQUAL(b4.empty(), false);
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.cbegin(), b1.cend(), b4.begin(), b4.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.begin(), b1.end(), b4.begin(), b4.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.rbegin(), b1.rend(), b4.rbegin(), b4.rend());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.crbegin(), b1.crend(), b4.crbegin(), b4.crend());
  BOOST_CHECK_EQUAL(b2.size(), 0);
  BOOST_CHECK_EQUAL(b2.data(), nullptr);
  BOOST_CHECK_EQUAL(b2.empty(), true);

  // move assign
  auto b5 = std::move(b3);
  BOOST_CHECK_EQUAL(b1.size(), b5.size());
  BOOST_CHECK_EQUAL(b1.data(), b5.data());
  BOOST_CHECK_EQUAL(b5.empty(), false);
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.cbegin(), b1.cend(), b5.begin(), b5.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.begin(), b1.end(), b5.begin(), b5.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.rbegin(), b1.rend(), b5.rbegin(), b5.rend());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.crbegin(), b1.crend(), b5.crbegin(), b5.crend());
  BOOST_CHECK_EQUAL(b3.size(), 0);
  BOOST_CHECK_EQUAL(b3.data(), nullptr);
  BOOST_CHECK_EQUAL(b3.empty(), true);

  // release
  auto s = std::move(b1).release();
  BOOST_CHECK_EQUAL(s.get(), nullptr);

  BOOST_CHECK_EQUAL(b1.size(), 0);
  BOOST_CHECK_EQUAL(b1.data(), nullptr);
  BOOST_CHECK_EQUAL(b1.empty(), true);
}

BOOST_AUTO_TEST_CASE(test_constructSelfManagedAttachableBuffer)
{
  std::array a = {0u, 1u, 2u, 3u, 4u};
  AttachableBuffer<uint32_t> b1{a.size()};

  std::copy(a.begin(), a.end(), b1.begin());

  BOOST_CHECK_EQUAL(b1.size(), a.size());
  BOOST_CHECK(b1.data() != nullptr);
  BOOST_CHECK(b1.data() != a.data());
  BOOST_CHECK_EQUAL(b1.empty(), false);
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.cbegin(), b1.cend(), a.begin(), a.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.begin(), b1.end(), a.begin(), a.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.rbegin(), b1.rend(), a.rbegin(), a.rend());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.crbegin(), b1.crend(), a.crbegin(), a.crend());

  for (size_t i = 0; i < a.size(); i++) {
    BOOST_CHECK_EQUAL(b1[i], a[i]);
  }

  // copy construct
  auto b2{b1};
  BOOST_CHECK_EQUAL(b1.size(), b2.size());
  BOOST_CHECK(b1.data() != b2.data());
  BOOST_CHECK_EQUAL(b2.empty(), false);
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.cbegin(), b1.cend(), b2.begin(), b2.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.begin(), b1.end(), b2.begin(), b2.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.rbegin(), b1.rend(), b2.rbegin(), b2.rend());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.crbegin(), b1.crend(), b2.crbegin(), b2.crend());

  // copy assign
  auto b3 = b1;
  BOOST_CHECK_EQUAL(b1.size(), b3.size());
  BOOST_CHECK(b1.data() != b3.data());
  BOOST_CHECK_EQUAL(b3.empty(), false);
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.cbegin(), b1.cend(), b3.begin(), b3.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.begin(), b1.end(), b3.begin(), b3.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.rbegin(), b1.rend(), b3.rbegin(), b3.rend());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.crbegin(), b1.crend(), b3.crbegin(), b3.crend());

  // move construct
  auto d = b2.data();
  auto b4{std::move(b2)};
  BOOST_CHECK_EQUAL(b1.size(), b4.size());
  BOOST_CHECK_EQUAL(b4.data(), d);
  BOOST_CHECK_EQUAL(b4.empty(), false);
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.cbegin(), b1.cend(), b4.begin(), b4.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.begin(), b1.end(), b4.begin(), b4.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.rbegin(), b1.rend(), b4.rbegin(), b4.rend());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.crbegin(), b1.crend(), b4.crbegin(), b4.crend());
  BOOST_CHECK_EQUAL(b2.size(), 0);
  BOOST_CHECK_EQUAL(b2.data(), nullptr);
  BOOST_CHECK_EQUAL(b2.empty(), true);

  // move assign
  d = b3.data();
  auto b5 = std::move(b3);
  BOOST_CHECK_EQUAL(b1.size(), b5.size());
  BOOST_CHECK_EQUAL(b5.data(), d);
  BOOST_CHECK_EQUAL(b5.empty(), false);
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.cbegin(), b1.cend(), b5.begin(), b5.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.begin(), b1.end(), b5.begin(), b5.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.rbegin(), b1.rend(), b5.rbegin(), b5.rend());
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.crbegin(), b1.crend(), b5.crbegin(), b5.crend());
  BOOST_CHECK_EQUAL(b3.size(), 0);
  BOOST_CHECK_EQUAL(b3.data(), nullptr);
  BOOST_CHECK_EQUAL(b3.empty(), true);

  // release
  d = b1.data();
  auto s = std::move(b1).release();
  BOOST_CHECK_EQUAL(s.get(), d);

  BOOST_CHECK_EQUAL(b1.size(), 0);
  BOOST_CHECK_EQUAL(b1.data(), nullptr);
  BOOST_CHECK_EQUAL(b1.empty(), true);
}

BOOST_AUTO_TEST_CASE(test_operationsAttachableBuffer)
{
  std::array a = {0u, 1u, 2u, 3u, 4u};
  AttachableBuffer<uint32_t> b1{a.size()};

  std::copy(a.begin(), a.end(), b1.begin());
  BOOST_CHECK_EQUAL(b1.size(), a.size());
  BOOST_CHECK(b1.data() != nullptr);
  BOOST_CHECK(b1.data() != a.data());
  BOOST_CHECK_EQUAL(b1.empty(), false);
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.begin(), b1.end(), a.begin(), a.end());

  const size_t newSize = 10;
  auto d = b1.data();
  b1.resize(newSize);
  BOOST_CHECK_EQUAL(b1.size(), newSize);
  BOOST_CHECK(b1.data() != nullptr);
  BOOST_CHECK(b1.data() != a.data());
  BOOST_CHECK(b1.data() != d);
  BOOST_CHECK_EQUAL(b1.empty(), false);
  BOOST_CHECK_EQUAL_COLLECTIONS(b1.cbegin(), b1.cbegin() + a.size(), a.begin(), a.end());

  AttachableBuffer<uint32_t> b2{a.data(), a.size()};
  BOOST_CHECK_EQUAL(b2.size(), a.size());
  BOOST_CHECK_EQUAL(b2.data(), a.data());
  BOOST_CHECK_EQUAL(b2.empty(), false);
  BOOST_CHECK_EQUAL_COLLECTIONS(b2.begin(), b2.end(), a.begin(), a.end());

  d = b2.data();
  b2.resize(newSize);
  BOOST_CHECK_EQUAL(b2.size(), newSize);
  BOOST_CHECK(b2.data() != nullptr);
  BOOST_CHECK(b2.data() != a.data());
  BOOST_CHECK(b2.data() != d);
  BOOST_CHECK_EQUAL(b2.empty(), false);
  BOOST_CHECK_EQUAL_COLLECTIONS(b2.cbegin(), b2.cbegin() + a.size(), a.begin(), a.end());

  b2.assign(5);
  for (const auto& elem : b2) {
    BOOST_CHECK_EQUAL(elem, 5);
  }

  b2.attach(a.data(), a.size());
  BOOST_CHECK_EQUAL(b2.size(), a.size());
  BOOST_CHECK_EQUAL(b2.data(), a.data());
  BOOST_CHECK_EQUAL(b2.empty(), false);
  BOOST_CHECK_EQUAL_COLLECTIONS(b2.begin(), b2.end(), a.begin(), a.end());
}
