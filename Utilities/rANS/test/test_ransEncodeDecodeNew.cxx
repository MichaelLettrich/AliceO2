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

#include <boost/test/unit_test.hpp>
#include <boost/mp11.hpp>

#include <gsl/span>

#include "rANS/typetraits.h"
#include "rANS/renorm.h"
#include "rANS/rans.h"

using namespace o2::rans;

struct Empty {
  inline static const std::string Data{};

  using iterator_type = decltype(Data.begin());
};

struct Full {

  inline static const std::string Data = R"(Sed ut perspiciatis, unde omnis iste natus error sit voluptatem accusantium
doloremque laudantium, totam rem aperiam eaque ipsa, quae ab illo inventore veritatis
et quasi architecto beatae vitae dicta sunt, explicabo. nemo enim ipsam voluptatem,
quia voluptas sit, aspernatur aut odit aut fugit, sed quia consequuntur magni dolores
eos, qui ratione voluptatem sequi nesciunt, neque porro quisquam est, qui dolorem ipsum,
quia dolor sit, amet, consectetur, adipisci velit, sed quia non numquam eius modi tempora
incidunt, ut labore et dolore magnam aliquam quaerat voluptatem. ut enim ad minima veniam,
quis nostrum exercitationem ullam corporis suscipit laboriosam, nisi ut aliquid ex ea
commodi consequatur? quis autem vel eum iure reprehenderit, qui in ea voluptate velit
esse, quam nihil molestiae consequatur, vel illum, qui dolorem eum fugiat,
quo voluptas nulla pariatur?)";

  using iterator_type = decltype(Data.begin());
};

// template <typename T>
// inline constexpr bool isIncompressibleOnlyTestCase = std::is_same_v<TestString<Empty, Full>, T>;

using string_types = boost::mp11::mp_list<Empty, Full>;

using container_types = boost::mp11::mp_list<std::integral_constant<ContainerTag, ContainerTag::Dynamic>,
                                             std::integral_constant<ContainerTag, ContainerTag::Static>>;

using coder_types = boost::mp11::mp_list<std::integral_constant<CoderTag, CoderTag::Compat>,
                                         std::integral_constant<CoderTag, CoderTag::SingleStream>,
                                         std::integral_constant<CoderTag, CoderTag::SSE>,
                                         std::integral_constant<CoderTag, CoderTag::AVX2>>;

using testCase_types = boost::mp11::mp_product<boost::mp11::mp_list, container_types, coder_types, string_types, string_types>;

using ransSource_type = char;
using ransCoder_type = uint64_t;
using ransStream_type = uint32_t;

inline constexpr size_t NRansStreams = o2::rans::NStreams;
inline constexpr size_t RansRenormingPrecision = 16;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_encodeDecode, test_types, testCase_types)
{
  using container_type = boost::mp11::mp_at_c<test_types, 0>;
  using coder_type = boost::mp11::mp_at_c<test_types, 1>;
  using dictString_type = boost::mp11::mp_at_c<test_types, 2>;
  using encodeString_type = boost::mp11::mp_at_c<test_types, 3>;

  constexpr ContainerTag containerTag = container_type::value;
  constexpr CoderTag coderTag = coder_type::value;
  const std::string& dictString = dictString_type::Data;
  const std::string& encodeString = encodeString_type::Data;

  //TODO(milettri): renorming is not satisfactory.
  size_t precision = dictString.size() == 0 ? 0 : RansRenormingPrecision;

  auto encoder = makeEncoder<coderTag>::template fromSamples<typename dictString_type::iterator_type, containerTag>(dictString.begin(), dictString.end(), precision);

  auto decoder = makeDecoder<>::fromSamples<typename dictString_type::iterator_type, containerTag>(dictString.begin(), dictString.end(), precision);

  if (dictString == encodeString) {
    std::vector<ransStream_type> encodeBuffer(encodeString.size());
    auto encodeBufferEnd = encoder.process(encodeString.begin(), encodeString.end(), encodeBuffer.begin());
    std::vector<ransStream_type> encodeBuffer2(encodeString.size());
    auto encodeBuffer2End = encoder.process(gsl::span<const ransSource_type>(encodeString), gsl::make_span(encodeBuffer2));

    BOOST_CHECK_EQUAL_COLLECTIONS(encodeBuffer.begin(), encodeBufferEnd, encodeBuffer2.data(), encodeBuffer2End);

    std::vector<ransSource_type> decodeBuffer(encodeString.size());
    decoder.process(encodeBufferEnd, decodeBuffer.begin(), encodeString.size(), encoder.getNStreams());

    BOOST_CHECK_EQUAL_COLLECTIONS(decodeBuffer.begin(), decodeBuffer.end(), encodeString.begin(), encodeString.end());
  }

  std::vector<ransSource_type> literals(encodeString.size());
  std::vector<ransStream_type> encodeBuffer(encodeString.size());
  auto [encodeBufferEnd, literalBufferEnd] = encoder.process(encodeString.begin(), encodeString.end(), encodeBuffer.begin(), literals.begin());
  std::vector<ransStream_type> encodeBuffer2(encodeString.size());
  std::vector<ransSource_type> literals2(encodeString.size());
  auto [encodeBuffer2End, literalBuffer2End] = encoder.process(gsl::span<const ransSource_type>(encodeString), gsl::make_span(encodeBuffer2), literals2.begin());

  BOOST_CHECK_EQUAL_COLLECTIONS(encodeBuffer.begin(), encodeBufferEnd, encodeBuffer2.data(), encodeBuffer2End);
  BOOST_CHECK_EQUAL_COLLECTIONS(literals.begin(), literalBufferEnd, literals2.begin(), literalBuffer2End);

  std::vector<ransSource_type> decodeBuffer(encodeString.size());
  decoder.process(encodeBufferEnd, decodeBuffer.begin(), encodeString.size(), encoder.getNStreams(), literalBufferEnd);

  BOOST_CHECK_EQUAL_COLLECTIONS(decodeBuffer.begin(), decodeBuffer.end(), encodeString.begin(), encodeString.end());
};
