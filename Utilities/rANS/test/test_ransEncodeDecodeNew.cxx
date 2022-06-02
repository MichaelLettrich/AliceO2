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
#include <boost/mpl/vector.hpp>
#include <boost/mp11.hpp>

#include <gsl/span>

#include "rANS/EncoderFacade.h"
#include "rANS/StaticFrequencyTable.h"
#include "rANS/DynamicFrequencyTable.h"
#include "rANS/HashFrequencyTable.h"
#include "rANS/RenormedFrequencies.h"
#include "rANS/StaticSymbolTable.h"
#include "rANS/DynamicSymbolTable.h"
#include "rANS/HashSymbolTable.h"
#include "rANS/renorm.h"
#include "rANS/internal/SingleStreamEncodeCommand.h"
#include "rANS/internal/SIMDEncodeCommand.h"

#include "rANS/SIMDDecoder.h"
#include "rANS/LiteralSIMDDecoder.h"

#include "rANS/FrequencyTable.h"
#include "rANS/rans.h"

using namespace o2::rans;

struct EmptyTestString {
  std::string data{};
};

struct FullTestString : public EmptyTestString {
  FullTestString()
  {
    data = R"(Sed ut perspiciatis, unde omnis iste natus error sit voluptatem accusantium
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
  }
};

// template <typename source_T,
//           template <typename> class frequencyTable_T,
//           template <typename> class renormedFrequencyTable_T,
//           template <typename, typename> class symbolTable_T>
// struct FrequencyTableBuilder {

//   using source_type = source_T;
//   using frequencyTable_type = frequencyTable_T<source_type>;
//   using renormedFrequencyTable_type = renormedFrequencyTable_T<source_type>;
//   template <typename symbol_T>
//   using symbolTable_type = symbolTable_T<source_type, symbol_T>;

//   inline static constexpr size_t renormingPrecision = 16;

//   template <typename source_IT>
//   static constexpr renormedFrequencyTable_type operator(source_IT begin, source_IT end)
//   {
//     frequencyTable_type fTable{};
//     fTable.addSamples(begin, end);
//     return renormCutoffIncompressible<>(fTable, renormingPrecision);
//   };
// };
//
// using staticFrequencyBuilder_type = FrequencyTableBuilder<char,
//                                                           StaticFrequencyTable,
//                                                           RenormedStaticFrequencyTable,
//                                                           StaticSymbolTable>;

// using dynamicFrequencyBuilder_type = FrequencyTableBuilder<char,
//                                                            DynamicFrequencyTable,
//                                                            RenormedDynamicFrequencyTable,
//                                                            DynamicSymbolTable>;

// using hashFrequencyBuilder_type = FrequencyTableBuilder<char,
//                                                         HashFrequencyTable,
//                                                         RenormedHashFrequencyTable,
//                                                         HashSymbolTable>;

// template <template <typename...> typename... F>
// using mp_list_q = boost::mp11::mp_list<boost::mp11::mp_quote<F>...>;

// using testString_types = boost::mp11::mp_list<EmptyTestString, FullTestString>;

// using encoder_types = boost::mp11::mp_list<
//   o2::rans::internal::CompatEncoderCommand<31>,
//   o2::rans::internal::SingleStreamEncoderCommand<31>>;

// using frequencyTableBuilders = boost::mp11::mp_list<staticFrequencyBuilder_type,
//                                                     dynamicFrequencyBuilder_type,
//                                                     hashFrequencyBuilder_type>;

// template <
//   class encoderCommand_T,
//   class frequencyBuilder_T,
//   class dictString_T, class testString_T>
// struct EncodeDecodeBase {
//  public:
//   static constexpr size_t NStreams = 2;

//   using encoderCommand_type = encoderCommand_T;
//   using encoderSymbol_type = typename encoderCommand_type::encoder_symbol;
//   using symbolTable_type = typename frequencyBuilder_T::symbolTable_type<encoderSymbol_type>;
//   using encoder_type = EncoderFacade<encoderCommand_type, symbolTable_type, NStreams>;

//   EncodeDecodeBase()
//   {

//     dictString_T source;
//     std::string& s = source.data;

//     auto renormedFrequencyTable = frequencyBuilder_T::(std::begin(s), std::end(s));

//     encoder = encoder_type{renormedFrequencyTable};
//     decoder = decoder_type{renormedFrequencyTable};
//   }

//   virtual void encode() = 0;
//   virtual void decode() = 0;

//   void check()
//   {
//     testString_T testString;
//     BOOST_CHECK_EQUAL_COLLECTIONS(testString.data.begin(), testString.data.end(), decodeBuffer.begin(), decodeBuffer.end());
//   }

//   testString_T source;
//   encoder_T<typename params_t::coder_t, typename params_t::stream_t, typename params_t::source_t> encoder{};
//   decoder_T<typename params_t::coder_t, typename params_t::stream_t, typename params_t::source_t> decoder{};
//   std::vector<typename Params<coder_T>::stream_t> encodeBuffer{};
//   std::vector<typename Params<coder_T>::source_t> decodeBuffer{};
// };

// using outer_types = mp_list_q<std::vector, std::set>;
// using inner_types = boost::mp11::mp_list<int, double>;

// using TypeList = boost::mp11::mp_product<boost::mp11::mp_invoke_q,
//                                          outer_types,
//                                          inner_types>;

// using encoderCommand_type = internal::SingleStreamEncoderCommand<31>;
// using symbolTable_type = StaticSymbolTable<char, typename encoderCommand_type::symbol_type>;

// using fastRANSEncoder_type = EncoderFacade<encoderCommand_type, symbolTable_type, 2>;

BOOST_FIXTURE_TEST_CASE(test_ModularEncoders, FullTestString)
{

  using namespace o2::rans;

  StaticFrequencyTable<char> frequencyTable;
  std::string empty;
  frequencyTable.addSamples(gsl::make_span(empty));

  auto renormed = renormCutoffIncompressible<>(frequencyTable);

  using encoderCommand_type = internal::SingleStreamEncoderCommand<31>;
  using symbolTable_type = StaticSymbolTable<char, typename encoderCommand_type::symbol_type>;

  using fastRANSEncoder_type = EncoderFacade<encoderCommand_type, symbolTable_type, 2>;

  fastRANSEncoder_type encoder{renormed};

  std::vector<uint32_t> encodeBuffer(1 << 16, 0);
  std::vector<char> incompressible{};
  const auto& sourceData = data;

  auto encodeEnd = encoder.process(data.data(), data.data() + data.size(), encodeBuffer.begin(), std::back_inserter(incompressible));

  auto fTable = renormCutoffIncompressible<>(makeFrequencyTableFromSamples(empty.begin(), empty.end()), renormed.getRenormingBits());
  std::vector<char> incompressible2{};
  LiteralEncoder64<char> encoder2{fTable};
  LiteralDecoder64<char> decoder{fTable};

  std::vector<uint32_t> encodeBuffer2(1 << 16, 0);
  auto encodeEnd2 = encoder2.process(data.data(), data.data() + data.size(), encodeBuffer2.begin(), incompressible2);

  BOOST_CHECK_EQUAL_COLLECTIONS(encodeBuffer.begin(), encodeEnd, encodeBuffer2.begin(), encodeEnd2);

  std::vector<char> decodeBuffer(data.size(), 0);
  decoder.process(encodeEnd, decodeBuffer.begin(), data.size(), incompressible);

  BOOST_CHECK_EQUAL_COLLECTIONS(decodeBuffer.begin(), decodeBuffer.end(), data.begin(), data.end());
};

BOOST_FIXTURE_TEST_CASE(test_ModularEncodersLiterals, FullTestString)
{

  using namespace o2::rans;

  StaticFrequencyTable<char> frequencyTable;
  frequencyTable.addSamples(gsl::make_span(data));

  auto renormed = renormCutoffIncompressible<>(frequencyTable);

  using encoderCommand_type = internal::SingleStreamEncoderCommand<31>;
  using symbolTable_type = StaticSymbolTable<char, typename encoderCommand_type::symbol_type>;

  using fastRANSEncoder_type = EncoderFacade<encoderCommand_type, symbolTable_type, 2>;

  fastRANSEncoder_type encoder{renormed};

  std::vector<uint32_t> encodeBuffer(1 << 16, 0);

  const auto& sourceData = data;

  auto encodeEnd = encoder.process(data.data(), data.data() + data.size(), encodeBuffer.begin());

  auto fTable = renormCutoffIncompressible<>(makeFrequencyTableFromSamples(data.begin(), data.end()), renormed.getRenormingBits());

  Encoder64<char> encoder2{fTable};
  Decoder64<char> decoder{fTable};

  std::vector<uint32_t> encodeBuffer2(1 << 16, 0);
  auto encodeEnd2 = encoder2.process(data.data(), data.data() + data.size(), encodeBuffer2.begin());

  BOOST_CHECK_EQUAL_COLLECTIONS(encodeBuffer.begin(), encodeEnd, encodeBuffer2.begin(), encodeEnd2);

  std::vector<char> decodeBuffer(data.size(), 0);
  decoder.process(encodeEnd, decodeBuffer.begin(), data.size());

  BOOST_CHECK_EQUAL_COLLECTIONS(decodeBuffer.begin(), decodeBuffer.end(), data.begin(), data.end());
};

BOOST_FIXTURE_TEST_CASE(test_SSEEncoder, FullTestString)
{

  using namespace o2::rans;

  StaticFrequencyTable<char> frequencyTable;
  frequencyTable.addSamples(gsl::make_span(data));

  auto renormed = renormCutoffIncompressible<>(frequencyTable);

  using encoderCommand_type = internal::SSEEncoderCommand<20>;
  using symbolTable_type = StaticSymbolTable<char, internal::simd::Symbol>;

  using encoder_type = EncoderFacade<encoderCommand_type, symbolTable_type, 4>;

  encoder_type encoder{renormed};

  std::vector<uint32_t> encodeBuffer(1 << 16, 0);

  const auto& sourceData = data;

  auto encodeEnd = encoder.process(data.data(), data.data() + data.size(), encodeBuffer.begin());

  auto fTable = renormCutoffIncompressible<>(makeFrequencyTableFromSamples(data.begin(), data.end()), renormed.getRenormingBits());

  SIMDDecoder<uint64_t, uint32_t, char, 4, 2> decoder{fTable};
  std::vector<char> decodeBuffer(data.size(), 0);
  decoder.process(encodeEnd, decodeBuffer.begin(), data.size());

  BOOST_CHECK_EQUAL_COLLECTIONS(decodeBuffer.begin(), decodeBuffer.end(), data.begin(), data.end());
};

BOOST_FIXTURE_TEST_CASE(test_SSEEncoderLiterals, FullTestString)
{

  using namespace o2::rans;

  StaticFrequencyTable<char> frequencyTable;
  frequencyTable.addSamples(gsl::make_span(data));

  auto renormed = renormCutoffIncompressible<>(frequencyTable);

  using encoderCommand_type = internal::SSEEncoderCommand<20>;
  using symbolTable_type = StaticSymbolTable<char, internal::simd::Symbol>;

  using encoder_type = EncoderFacade<encoderCommand_type, symbolTable_type, 4>;

  encoder_type encoder{renormed};

  std::vector<uint32_t> encodeBuffer(1 << 16, 0);
  std::vector<char> incompressible{};

  const auto& sourceData = data;

  auto encodeEnd = encoder.process(data.data(), data.data() + data.size(), encodeBuffer.begin(), std::back_inserter(incompressible));

  auto fTable = renormCutoffIncompressible<>(makeFrequencyTableFromSamples(data.begin(), data.end()), renormed.getRenormingBits());

  LiteralSIMDDecoder<uint64_t, uint32_t, char, 4, 2> decoder{fTable};
  std::vector<char> decodeBuffer(data.size(), 0);
  decoder.process(encodeEnd, decodeBuffer.begin(), data.size(), incompressible);

  BOOST_CHECK_EQUAL_COLLECTIONS(decodeBuffer.begin(), decodeBuffer.end(), data.begin(), data.end());
};

BOOST_FIXTURE_TEST_CASE(test_AVXEncoder, FullTestString)
{

  using namespace o2::rans;

  StaticFrequencyTable<char> frequencyTable;
  frequencyTable.addSamples(gsl::make_span(data));

  auto renormed = renormCutoffIncompressible<>(frequencyTable);

  using encoderCommand_type = internal::AVXEncoderCommand<20>;
  using symbolTable_type = StaticSymbolTable<char, internal::simd::Symbol>;

  using encoder_type = EncoderFacade<encoderCommand_type, symbolTable_type, 8>;

  encoder_type encoder{renormed};

  std::vector<uint32_t> encodeBuffer(1 << 16, 0);

  const auto& sourceData = data;

  auto encodeEnd = encoder.process(data.data(), data.data() + data.size(), encodeBuffer.begin());

  auto fTable = renormCutoffIncompressible<>(makeFrequencyTableFromSamples(data.begin(), data.end()), renormed.getRenormingBits());

  SIMDDecoder<uint64_t, uint32_t, char, 8, 4> decoder{fTable};
  std::vector<char> decodeBuffer(data.size(), 0);
  decoder.process(encodeEnd, decodeBuffer.begin(), data.size());

  BOOST_CHECK_EQUAL_COLLECTIONS(decodeBuffer.begin(), decodeBuffer.end(), data.begin(), data.end());
};

BOOST_FIXTURE_TEST_CASE(test_AVXEncoderLiterals, FullTestString)
{

  using namespace o2::rans;

  StaticFrequencyTable<char> frequencyTable;
  frequencyTable.addSamples(gsl::make_span(data));

  auto renormed = renormCutoffIncompressible<>(frequencyTable);

  using encoderCommand_type = internal::AVXEncoderCommand<20>;
  using symbolTable_type = StaticSymbolTable<char, internal::simd::Symbol>;

  using encoder_type = EncoderFacade<encoderCommand_type, symbolTable_type, 8>;

  encoder_type encoder{renormed};

  std::vector<uint32_t> encodeBuffer(1 << 16, 0);
  std::vector<char> incompressible{};

  const auto& sourceData = data;

  auto encodeEnd = encoder.process(data.data(), data.data() + data.size(), encodeBuffer.begin(), std::back_inserter(incompressible));

  auto fTable = renormCutoffIncompressible<>(makeFrequencyTableFromSamples(data.begin(), data.end()), renormed.getRenormingBits());

  LiteralSIMDDecoder<uint64_t, uint32_t, char, 8, 4> decoder{fTable};
  std::vector<char> decodeBuffer(data.size(), 0);
  decoder.process(encodeEnd, decodeBuffer.begin(), data.size(), incompressible);

  BOOST_CHECK_EQUAL_COLLECTIONS(decodeBuffer.begin(), decodeBuffer.end(), data.begin(), data.end());
};

// struct Params {
//   using stream_type = uint32_t;
//   using source_type = char;
//   static constexpr size_t symbolTablePrecision = 16;
// };

// template <
//   template <typename, typename, typename> class encoder_T,
//   template <typename, typename, typename> class decoder_T,
//   typename coder_T, class dictString_T, class testString_T>
// struct EncodeDecodeBase {
//  public:
//   using params_t = Params<coder_T>;

//   EncodeDecodeBase()
//   {
//     dictString_T source;
//     std::string& s = source.data;
//     o2::rans::RenormedFrequencyTable frequencyTable = o2::rans::renorm(o2::rans::makeFrequencyTableFromSamples(std::begin(s), std::end(s)), params_t::symbolTablePrecision);

//     encoder = decltype(encoder)(frequencyTable);
//     decoder = decltype(decoder)(frequencyTable);

//     const auto [min, max] = [&s]() {
//       const auto [minIter, maxIter] = std::minmax_element(s.begin(), s.end());
//       const char min = minIter == s.end() ? 0 : *minIter;
//       const char max = maxIter == s.end() ? 0 : *maxIter + 1;
//       return std::make_tuple(min, max);
//     }();

//     const size_t alphabetRangeBits = o2::rans::internal::numBitsForNSymbols(max - min + 1 + 1);

//     BOOST_CHECK_EQUAL(encoder.getSymbolTablePrecision(), params_t::symbolTablePrecision);
//     BOOST_CHECK_EQUAL(encoder.getAlphabetRangeBits(), alphabetRangeBits);
//     BOOST_CHECK_EQUAL(encoder.getMinSymbol(), min);
//     BOOST_CHECK_EQUAL(encoder.getMaxSymbol(), max);

//     BOOST_CHECK_EQUAL(decoder.getSymbolTablePrecision(), params_t::symbolTablePrecision);
//     BOOST_CHECK_EQUAL(decoder.getAlphabetRangeBits(), alphabetRangeBits);
//     BOOST_CHECK_EQUAL(decoder.getMinSymbol(), min);
//     BOOST_CHECK_EQUAL(decoder.getMaxSymbol(), max);
//   }

//   virtual void encode() = 0;
//   virtual void decode() = 0;

//   void check()
//   {
//     testString_T testString;
//     BOOST_CHECK_EQUAL_COLLECTIONS(testString.data.begin(), testString.data.end(), decodeBuffer.begin(), decodeBuffer.end());
//   }

//   testString_T source;
//   encoder_T<typename params_t::coder_t, typename params_t::stream_t, typename params_t::source_t> encoder{};
//   decoder_T<typename params_t::coder_t, typename params_t::stream_t, typename params_t::source_t> decoder{};
//   std::vector<typename Params<coder_T>::stream_t> encodeBuffer{};
//   std::vector<typename Params<coder_T>::source_t> decodeBuffer{};
// };

// template <typename coder_T, class dictString_T, class testString_T>
// struct EncodeDecode : public EncodeDecodeBase<o2::rans::Encoder, o2::rans::Decoder, coder_T, dictString_T, testString_T> {
//   void encode() override
//   {
//     BOOST_CHECK_NO_THROW(this->encoder.process(std::begin(this->source.data), std::end(this->source.data), std::back_inserter(this->encodeBuffer)));
//   };
//   void decode() override
//   {
//     BOOST_CHECK_NO_THROW(this->decoder.process(this->encodeBuffer.end(), std::back_inserter(this->decodeBuffer), this->source.data.size()));
//   };
// };

// template <typename coder_T, typename stream_T, typename source_V>
// using simdEncoderSSE_t = o2::rans::SIMDEncoder<coder_T, stream_T, source_V, 4, 2>;

// template <typename coder_T, typename stream_T, typename source_V>
// using simdDecoderSSE_t = o2::rans::SIMDDecoder<coder_T, stream_T, source_V, 4, 2>;

// template <typename coder_T, class dictString_T, class testString_T>
// struct EncodeDecodeSSE : public EncodeDecodeBase<simdEncoderSSE_t, simdDecoderSSE_t, coder_T, dictString_T, testString_T> {
//   void encode() override
//   {
//     BOOST_CHECK_NO_THROW(this->encoder.process(std::begin(this->source.data), std::end(this->source.data), std::back_inserter(this->encodeBuffer)));
//   };
//   void decode() override
//   {
//     BOOST_CHECK_NO_THROW(this->decoder.process(this->encodeBuffer.end(), std::back_inserter(this->decodeBuffer), this->source.data.size()));
//   };
// };

// template <typename coder_T, typename stream_T, typename source_V>
// using simdEncoderAVX_t = o2::rans::SIMDEncoder<coder_T, stream_T, source_V, 8, 4>;

// template <typename coder_T, typename stream_T, typename source_V>
// using simdDecoderAVX_t = o2::rans::SIMDDecoder<coder_T, stream_T, source_V, 8, 4>;

// template <typename coder_T, class dictString_T, class testString_T>
// struct EncodeDecodeAVX : public EncodeDecodeBase<simdEncoderAVX_t, simdDecoderAVX_t, coder_T, dictString_T, testString_T> {
//   void encode() override
//   {
//     BOOST_CHECK_NO_THROW(this->encoder.process(std::begin(this->source.data), std::end(this->source.data), std::back_inserter(this->encodeBuffer)));
//   };
//   void decode() override
//   {
//     BOOST_CHECK_NO_THROW(this->decoder.process(this->encodeBuffer.end(), std::back_inserter(this->decodeBuffer), this->source.data.size()));
//   };
// };

// template <typename coder_T, class dictString_T, class testString_T>
// struct EncodeDecodeLiteral : public EncodeDecodeBase<o2::rans::LiteralEncoder, o2::rans::LiteralDecoder, coder_T, dictString_T, testString_T> {
//   void encode() override
//   {
//     BOOST_CHECK_NO_THROW(this->encoder.process(std::begin(this->source.data), std::end(this->source.data), std::back_inserter(this->encodeBuffer), literals));
//   };
//   void decode() override
//   {
//     BOOST_CHECK_NO_THROW(this->decoder.process(this->encodeBuffer.end(), std::back_inserter(this->decodeBuffer), this->source.data.size(), literals));
//     BOOST_CHECK(literals.empty());
//   };

//   std::vector<typename Params<coder_T>::source_t> literals;
// };

// template <typename coder_T, typename stream_T, typename source_V>
// using literalSimdEncoderSSE_t = o2::rans::LiteralSIMDEncoder<coder_T, stream_T, source_V, 4, 2>;

// template <typename coder_T, typename stream_T, typename source_V>
// using literalSimdDecoderSSE_t = o2::rans::LiteralSIMDDecoder<coder_T, stream_T, source_V, 4, 2>;

// template <typename coder_T, class dictString_T, class testString_T>
// struct EncodeDecodeLiteralSSE : public EncodeDecodeBase<literalSimdEncoderSSE_t, literalSimdDecoderSSE_t, coder_T, dictString_T, testString_T> {
//   void encode() override
//   {
//     BOOST_CHECK_NO_THROW(this->encoder.process(std::begin(this->source.data), std::end(this->source.data), std::back_inserter(this->encodeBuffer), literals));
//   };
//   void decode() override
//   {
//     BOOST_CHECK_NO_THROW(this->decoder.process(this->encodeBuffer.end(), std::back_inserter(this->decodeBuffer), this->source.data.size(), literals));
//     BOOST_CHECK(literals.empty());
//   };

//   std::vector<typename Params<coder_T>::source_t> literals;
// };

// template <typename coder_T, typename stream_T, typename source_V>
// using literalSimdEncoderAVX_t = o2::rans::LiteralSIMDEncoder<coder_T, stream_T, source_V, 8, 4>;

// template <typename coder_T, typename stream_T, typename source_V>
// using literalSimdDecoderAVX_t = o2::rans::LiteralSIMDDecoder<coder_T, stream_T, source_V, 8, 4>;

// template <typename coder_T, class dictString_T, class testString_T>
// struct EncodeDecodeLiteralAVX : public EncodeDecodeBase<literalSimdEncoderAVX_t, literalSimdDecoderAVX_t, coder_T, dictString_T, testString_T> {
//   void encode() override
//   {
//     BOOST_CHECK_NO_THROW(this->encoder.process(std::begin(this->source.data), std::end(this->source.data), std::back_inserter(this->encodeBuffer), literals));
//   };
//   void decode() override
//   {
//     BOOST_CHECK_NO_THROW(this->decoder.process(this->encodeBuffer.end(), std::back_inserter(this->decodeBuffer), this->source.data.size(), literals));
//     BOOST_CHECK(literals.empty());
//   };

//   std::vector<typename Params<coder_T>::source_t> literals;
// };

// template <typename coder_T, class dictString_T, class testString_T>
// struct EncodeDecodeDedup : public EncodeDecodeBase<o2::rans::DedupEncoder, o2::rans::DedupDecoder, coder_T, dictString_T, testString_T> {
//   void encode() override
//   {
//     BOOST_CHECK_NO_THROW(this->encoder.process(std::begin(this->source.data), std::end(this->source.data), std::back_inserter(this->encodeBuffer), duplicates));
//   };
//   void decode() override
//   {
//     BOOST_CHECK_NO_THROW(this->decoder.process(this->encodeBuffer.end(), std::back_inserter(this->decodeBuffer), this->source.data.size(), duplicates));
//   };

//   using params_t = Params<coder_T>;
//   typename o2::rans::DedupEncoder<typename params_t::coder_t,
//                                   typename params_t::stream_t,
//                                   typename params_t::source_t>::duplicatesMap_t duplicates;
// };

// using testCase_t = boost::mpl::vector<EncodeDecode<uint32_t, EmptyTestString, EmptyTestString>,
//                                       EncodeDecode<uint64_t, EmptyTestString, EmptyTestString>,
//                                       EncodeDecode<uint32_t, FullTestString, FullTestString>,
//                                       EncodeDecode<uint64_t, FullTestString, FullTestString>,
//                                       EncodeDecodeSSE<uint64_t, EmptyTestString, EmptyTestString>,
//                                       EncodeDecodeSSE<uint64_t, FullTestString, FullTestString>,
//                                       EncodeDecodeAVX<uint64_t, EmptyTestString, EmptyTestString>,
//                                       EncodeDecodeAVX<uint64_t, FullTestString, FullTestString>,
//                                       EncodeDecodeLiteral<uint32_t, EmptyTestString, EmptyTestString>,
//                                       EncodeDecodeLiteral<uint64_t, EmptyTestString, EmptyTestString>,
//                                       EncodeDecodeLiteral<uint32_t, FullTestString, FullTestString>,
//                                       EncodeDecodeLiteral<uint64_t, FullTestString, FullTestString>,
//                                       EncodeDecodeLiteral<uint32_t, EmptyTestString, FullTestString>,
//                                       EncodeDecodeLiteral<uint64_t, EmptyTestString, FullTestString>,
//                                       EncodeDecodeLiteralSSE<uint64_t, EmptyTestString, EmptyTestString>,
//                                       EncodeDecodeLiteralSSE<uint64_t, FullTestString, FullTestString>,
//                                       EncodeDecodeLiteralSSE<uint64_t, EmptyTestString, FullTestString>,
//                                       EncodeDecodeLiteralAVX<uint64_t, EmptyTestString, EmptyTestString>,
//                                       EncodeDecodeLiteralAVX<uint64_t, FullTestString, FullTestString>,
//                                       EncodeDecodeLiteralAVX<uint64_t, EmptyTestString, FullTestString>>;
// // EncodeDecodeDedup<uint32_t, EmptyTestString, EmptyTestString>,
// // EncodeDecodeDedup<uint64_t, EmptyTestString, EmptyTestString>,
// // EncodeDecodeDedup<uint32_t, FullTestString, FullTestString>,
// // EncodeDecodeDedup<uint64_t, FullTestString, FullTestString>>;

// BOOST_AUTO_TEST_CASE_TEMPLATE(test_encodeDecode, testCase_T, testCase_t)
// {
//   testCase_T testCase;
//   testCase.encode();
//   testCase.decode();
//   testCase.check();
// };
