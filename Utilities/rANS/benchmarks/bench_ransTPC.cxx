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

/// @file   bench_ransTPC.h
/// @author Michael Lettrich
/// @brief  main benchmark executable that performs compression/decompression of TPC data while taking detailed metrics.

#include <vector>
#include <cstring>
#include <execution>
#include <iterator>
#include <iostream>

#include <boost/program_options.hpp>
#include <boost/mp11.hpp>
#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/ostreamwrapper.h>
#include <fairlogger/Logger.h>

#include "rANS/factory.h"
#include "rANS/histogram.h"
#include "rANS/serialize.h"

#include "utils/TPCData.h"
#include "utils/Buffers.h"

#ifdef ENABLE_VTUNE_PROFILER
#include <ittnotify.h>
#endif

namespace bpo = boost::program_options;
namespace bmp = boost::mp11;
namespace rbm = o2::rans::benchmarks;

using json_writer_type = rapidjson::Writer<rapidjson::OStreamWrapper>;

static json_writer_type writer;

using namespace o2::rans;

using coder_types = bmp::mp_list<std::integral_constant<CoderTag, CoderTag::Compat>
#ifdef RANS_SINGLE_STREAM
                                 ,
                                 std::integral_constant<CoderTag, CoderTag::SingleStream>
#endif /* RANS_SINGLE_STREAM */
#ifdef RANS_SSE
                                 ,
                                 std::integral_constant<CoderTag, CoderTag::SSE>
#endif /* RANS_SSE */
#ifdef RANS_AVX2
                                 ,
                                 std::integral_constant<CoderTag, CoderTag::AVX2>
#endif /* RANS_AVX2 */
                                 >;

std::string toString(CoderTag tag)
{
  switch (tag) {
    case CoderTag::Compat:
      return {"Compat"};
      break;
    case CoderTag::SingleStream:
      return {"SingleStream"};
      break;
    case CoderTag::SSE:
      return {"SSE"};
      break;
    case CoderTag::AVX2:
      return {"AVX2"};
      break;
    default:
      throw Exception("Invalid");
      break;
  };
};

class TimingDecorator
{
 public:
  using jsonWriter_type = rapidjson::Writer<rapidjson::OStreamWrapper>;

  TimingDecorator() = default;
  TimingDecorator(jsonWriter_type& writer) : mWriter{&writer} {};

  template <typename F, std::enable_if_t<!std::is_void_v<std::invoke_result_t<F>>, bool> = true>
  decltype(auto) timeAndLog(const std::string& keyName, const std::string& logMessage, F functor)
  {
    mTimer.start();
    decltype(auto) ret = functor();
    mTimer.stop();
    mWriter->Key(keyName.c_str());
    const double_t msDuration = mTimer.getDurationMS();
    mWriter->Double(msDuration);
    LOGP(info, "{} in {} ms", logMessage, msDuration);

    return ret;
  };

  template <typename F, std::enable_if_t<std::is_void_v<std::invoke_result_t<F>>, bool> = true>
  void timeAndLog(const std::string& keyName, const std::string& logMessage, F functor)
  {
    mTimer.start();
    functor();
    mTimer.stop();
    mWriter->Key(keyName.c_str());
    const double_t msDuration = mTimer.getDurationMS();
    mWriter->Double(msDuration);
    LOGP(info, "{} in {} ms", logMessage, msDuration);
  };

 private:
  jsonWriter_type* mWriter{};
  utils::RANSTimer mTimer{};
};

template <rbm::TPCField field_V, CoderTag coderTag_V>
void ransEncodeDecode(const rbm::TPCData& dataIn, rbm::TPCData& dataOut)
{
  static constexpr auto fieldTrait = rbm::traits::FieldTrait<field_V>{};
  using source_type = typename decltype(fieldTrait)::value_type;

  auto [inputBegin, inputEnd] = rbm::getIter<field_V>(dataIn);
  const size_t inputSize = std::distance(inputBegin, inputEnd);

  if (inputSize == 0) {
    LOG(warning) << "Skipping 0 Byte field " << fieldTrait.name;
    return;
  }

  utils::RANSTimer timer{};
  TimingDecorator t{writer};

  writer.Key(fieldTrait.name.data());
  writer.StartObject();

  rbm::EncodeBuffer<source_type> encodeBuffer{inputSize};
  encodeBuffer.literals.resize(inputSize, 0);
  encodeBuffer.literalsEnd = encodeBuffer.literals.data();

  auto [outputBegin, outputEnd] = rbm::getIter<field_V>(dataOut);

  writer.Key("Timing");
  writer.StartObject();

  LOGP(info, "processing: {} (nItems: {}, size: {} MiB)", fieldTrait.name, inputSize, inputSize * sizeof(source_type) / 1024.0 / 1024.0);
  auto histogram = t.timeAndLog(
    "FrequencyTable", "Built Frequency Table", [&]() { return makeDenseHistogram::fromSamples(inputBegin, inputEnd); });

  auto tmpHist = histogram;
  Metrics<source_type> metrics{};
  RenormedDenseHistogram<source_type> renormedHistogram{};
  t.timeAndLog("Renorming", "Renormed Frequency Table", [&]() mutable {
    metrics = Metrics<source_type>{histogram};
    renormedHistogram = renorm(std::move(tmpHist), metrics);
  });

  auto encoder = t.timeAndLog("Encoder", "Built Encoder", [&]() { return makeDenseEncoder<coderTag_V>::fromRenormed(renormedHistogram); });

  t.timeAndLog("Encoding", "Encoded", [&]() mutable {
#ifdef ENABLE_VTUNE_PROFILER
    __itt_resume();
#endif
    if (renormedHistogram.hasIncompressibleSymbol()) {
      std::tie(encodeBuffer.encodeBufferEnd, encodeBuffer.literalsEnd) = encoder.process(inputBegin, inputEnd, encodeBuffer.buffer.data(), encodeBuffer.literalsEnd);
    } else {
      encodeBuffer.encodeBufferEnd = encoder.process(inputBegin, inputEnd, encodeBuffer.buffer.data());
    }
#ifdef ENABLE_VTUNE_PROFILER
    __itt_pause();
#endif
  });
  LOGP(info, "Encoded {} Bytes", inputSize * sizeof(source_type));

  std::vector<uint8_t> dict(histogram.size() * sizeof(uint64_t), 0);
  auto dictEnd = t.timeAndLog("WriteDict", "Serialized Dict", [&]() { return compressRenormedDictionary(encoder.getSymbolTable(), dict.data()); });
  LOGP(info, "Serialized Dict of {} Bytes", std::distance(dict.data(), dictEnd));
  auto recoveredHistogram = t.timeAndLog("ReadDict", "Read Dict", [&]() {
    const source_type min =  encoder.getSymbolTable().getOffset();
    const source_type max =  min + std::max<source_type>(static_cast<int64_t>(encoder.getSymbolTable().size())-1,0);
    return readRenormedDictionary(dict.data(), dictEnd,min,max,renormedHistogram.getRenormingBits()); });
  auto decoder = makeDecoder<>::fromRenormed(renormedHistogram);
  auto recoveredDecoder = t.timeAndLog("MakeDecoder", "Built Decoder", [&]() { return makeDecoder<>::fromRenormed(recoveredHistogram); });

  t.timeAndLog("Decoding", "Decoded", [&]() mutable {
    if (encodeBuffer.literalsEnd == encodeBuffer.literals.data()) {
      decoder.process(encodeBuffer.encodeBufferEnd, outputBegin, inputSize, encoder.getNStreams());
    } else {
      decoder.process(encodeBuffer.encodeBufferEnd, outputBegin, inputSize, encoder.getNStreams(), encodeBuffer.literalsEnd);
    }
  });

  // Verify correctness
  auto [checkBegin, checkEnd] = rbm::getIter<field_V>(static_cast<const rbm::TPCData&>(dataOut));

  if (!(std::equal(inputBegin, inputEnd, checkBegin))) {
    LOGP(warning, "Missmatch between original and decoded Message");
  }
  LOG(info) << "finished: " << fieldTrait.name;

  writer.EndObject(); // Timing

  const auto& datasetProperties = metrics.getDatasetProperties();

  // Frequency Table
  // ##########################
  writer.Key("FrequencyTable");
  writer.StartObject();
  writer.Key("nSamples");
  writer.Uint64(histogram.getNumSamples());
  writer.Key("Min");
  writer.Int(datasetProperties.min);
  writer.Key("Max");
  writer.Int(datasetProperties.max);
  writer.Key("alphabetRangeBits");
  writer.Int(datasetProperties.alphabetRangeBits);
  writer.Key("nUsedAlphabetSymbols");
  writer.Uint(datasetProperties.nUsedAlphabetSymbols);
  writer.Key("IncompressibleFrequency");
  writer.Uint(0);
  writer.EndObject(); // FrequencyTable

  // RescaledFrequencies
  // ##########################
  const Metrics<source_type> renormedMetrics{histogram};
  const auto& renormedDatasetProperties = renormedMetrics.getDatasetProperties();

  writer.Key("RescaledFrequencies");
  writer.StartObject();
  writer.Key("nSamples");
  writer.Uint64(renormedHistogram.getNumSamples());
  writer.Key("Min");
  writer.Int(renormedDatasetProperties.min);
  writer.Key("Max");
  writer.Int(renormedDatasetProperties.max);
  writer.Key("alphabetRangeBits");
  writer.Int(renormedDatasetProperties.alphabetRangeBits);
  writer.Key("nUsedAlphabetSymbols");
  writer.Uint(renormedDatasetProperties.nUsedAlphabetSymbols);
  writer.Key("IncompressibleFrequency");
  writer.Uint(renormedHistogram.getIncompressibleSymbolFrequency());
  writer.Key("RenormingBits");
  writer.Uint(renormedHistogram.getRenormingBits());
  writer.EndObject(); // RescaledFrequencies

  // Message Properties
  // ##########################
  writer.Key("Message");
  writer.StartObject();
  writer.Key("Size");
  writer.Uint64(inputSize);
  writer.Key("SymbolSize");
  writer.Uint(sizeof(source_type));
  writer.Key("Entropy");
  writer.Double(datasetProperties.entropy);
  writer.Key("ExpectedCodewordLength");
  writer.Double(computeExpectedCodewordLength<>(histogram, renormedHistogram));
  writer.EndObject(); // Message

  // Compression Properties
  // ##########################
  writer.Key("Compression");
  writer.StartObject();
  writer.Key("EncodeBufferSize");
  writer.Uint64(std::distance(encodeBuffer.buffer.data(), encodeBuffer.encodeBufferEnd) * sizeof(uint32_t));
  writer.Key("LiteralSize");
  writer.Uint64(std::distance(encodeBuffer.literals.data(), encodeBuffer.literalsEnd) * sizeof(source_type));
  writer.Key("DictSize");
  writer.Uint64(std::distance(dict.data(), dictEnd));
  writer.EndObject(); // Compression

  writer.EndObject(); // Encode/Decode Run
};

template <CoderTag coderTag_V, rbm::TPCField... fields_V>
void tpcEncodeDecode(const std::string& name, const rbm::TPCData& dataIn, rbm::TPCData& dataOut, bmp::mp_list_v<fields_V...>)
{
  writer.Key(name.c_str());
  writer.StartObject();

  (ransEncodeDecode<fields_V, coderTag_V>(dataIn, dataOut), ...);

  writer.EndObject();
};

using encoder_types = bmp::mp_product<bmp::mp_list, coder_types>;

int main(int argc, char* argv[])
{
  bpo::options_description options("Allowed options");
  // clang-format off
  options.add_options()
    ("help,h", "print usage message")
    ("in,i",bpo::value<std::string>(), "file to process")
    ("out,o",bpo::value<std::string>(), "json output file")
    ("log_severity,l",bpo::value<std::string>(), "severity of FairLogger");
  // clang-format on

  bpo::variables_map vm;
  bpo::store(bpo::parse_command_line(argc, argv, options), vm);
  bpo::notify(vm);

  if (vm.count("help")) {
    std::cout << options << "\n";
    return 0;
  }

  const std::string inFile = [&]() {
    if (vm.count("in")) {
      return vm["in"].as<std::string>();
    } else {
      LOG(error) << "missing path to input file";
      exit(1);
    }
  }();

  const std::string outFile = [&]() {
    if (vm.count("out")) {
      return vm["out"].as<std::string>();
    } else {
      return std::string("out.json");
    }
  }();

  if (vm.count("log_severity")) {
    fair::Logger::SetConsoleSeverity(vm["log_severity"].as<std::string>().c_str());
  }

  std::ofstream of{outFile};
  if (!of) {
    std::runtime_error(fmt::format("could not open output file at path {}", inFile));
  }

  rapidjson::OStreamWrapper stream{of};
  writer.Reset(stream);
  writer.StartObject();

  rbm::TPCData tpcDataIn = rbm::loadTPCData(inFile);
  rbm::TPCData tpcDataOut = rbm::TPCData{tpcDataIn.nTracks,
                                         tpcDataIn.nAttachedClusters,
                                         tpcDataIn.nUnattachedClusters,
                                         tpcDataIn.nAttachedClustersReduced,
                                         tpcDataIn.nTriggers};

  LOGP(info, "loaded Compressed Clusters from {}", inFile);
  LOG(info) << "######################################################";
  bmp::mp_for_each<encoder_types>([&](auto L) {
    using coder_type = bmp::mp_at_c<decltype(L), 0>;
    constexpr CoderTag coderTag = coder_type::value;
    const std::string encoderTitle = toString(coderTag);

    LOGP(info, "start rANS {}/Decode", encoderTitle);

    tpcEncodeDecode<coderTag>(encoderTitle, tpcDataIn, tpcDataOut, rbm::traits::TPCFieldsAll{});

    LOG(info) << "######################################################";
  });
  writer.EndObject();
  stream.Flush();
  of.close();
};
