// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   bench_ransTPCCompressionMetrics.cxx
/// @author Michael Lettrich
/// @brief  Computes compression metrics for TPC data for rANS entropy coding

#include <string>
#include <vector>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <regex>
#include <iostream>
#include <tuple>
#include <ranges>

#include <boost/program_options.hpp>
#include <boost/mp11.hpp>
#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/ostreamwrapper.h>
#include <fairlogger/Logger.h>

#include "utils/TPCData.h"
#include "utils/CompressedField.h"
#include "utils/Buffers.h"

#include "rANS/factory.h"
#include "rANS/metrics.h"
#include "rANS/encode.h"
#include "rANS/decode.h"
#include "rANS/utils.h"
#include "rANS/serialize.h"
#include "rANS/histogram.h"
#include "rANS/pack.h"

namespace bpo = boost::program_options;
namespace bmp = boost::mp11;
namespace rbm = o2::rans::benchmarks;
namespace fs = std::filesystem;

using json_writer_type = rapidjson::Writer<rapidjson::OStreamWrapper>;

static json_writer_type writer;

template <typename histogram_T, typename renormed_histogram_T>
void writeMessageProperties(const histogram_T& histogram, const renormed_histogram_T& renormedHistogram, const Metrics<typename histogram_T::source_type>& metrics)
{
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

  auto [renormedMin, renormedMax] = getMinMax(renormedHistogram);
  // RescaledFrequencies
  // ##########################
  writer.Key("RescaledFrequencies");
  writer.StartObject();
  writer.Key("nSamples");
  writer.Uint64(renormedHistogram.getNumSamples());
  writer.Key("Min");
  writer.Int(renormedMin);
  writer.Key("Max");
  writer.Int(renormedMax);
  writer.Key("alphabetRangeBits");
  writer.Int(utils::getRangeBits(renormedMin, renormedMax));
  writer.Key("nUsedAlphabetSymbols");
  writer.Uint(countNUsedAlphabetSymbols(renormedHistogram));
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
  writer.Uint64(datasetProperties.nSamples);
  writer.Key("SymbolSize");
  writer.Uint(sizeof(typename histogram_T::source_type));
  writer.Key("Entropy");
  writer.Double(datasetProperties.entropy);
  writer.Key("ExpectedCodewordLength");
  writer.Double(computeExpectedCodewordLength<>(histogram, renormedHistogram));
  writer.EndObject(); // Message
}

template <rbm::TPCField field_V>
void ransEncodeDecode(const rbm::TPCData& dataIn, rbm::TPCData& dataOut)
{
  using namespace o2::rans;

  const auto fieldTrait = rbm::traits::FieldTrait<field_V>{};
  using source_type = typename decltype(fieldTrait)::value_type;

  auto [inputBegin, inputEnd] = rbm::getIter<field_V>(dataIn);
  const size_t inputSize = std::distance(inputBegin, inputEnd);

  if (inputSize == 0) {
    LOGP(warning, "Skipping 0 Byte field {}", fieldTrait.name);
    return;
  }

  std::apply([&](auto... memberPtrs) {
    (std::fill((dataOut.*memberPtrs).begin(), (dataOut.*memberPtrs).end(), source_type{0}), ...);
  },
             fieldTrait.members);

  LOGP(info, "Processing field: {}", fieldTrait.name);
  writer.Key(fieldTrait.name.data());
  writer.StartObject();

  rbm::EncodeBuffer<source_type> encodeBuffer{inputSize};
  encodeBuffer.literals.resize(inputSize, 0);
  encodeBuffer.literalsEnd = encodeBuffer.literals.data();

  auto [outputBegin, outputEnd] = rbm::getIter<field_V>(dataOut);

  auto histogram = [&]() { 
    if constexpr (sizeof(source_type) < sizeof(uint32_t)) {
      return makeDenseHistogram::fromSamples(inputBegin, inputEnd);
    } else{
      return makeSparseHistogram::fromSamples(inputBegin, inputEnd);
    }
    return makeDenseHistogram::fromSamples(inputBegin, inputEnd); }();
  Metrics<source_type> metrics{histogram};
  writeFrequencyTableStats(histogram, metrics);

  auto renormedHistogram = renorm(histogram, metrics);
  writeRescaledFrequencyStats(renormedHistogram);
  writeMessageStats(histogram, renormedHistogram, metrics);

  histogram = {}; // free memory

  auto encoder = makeDenseEncoder<>::fromRenormed(renormedHistogram);
  std::tie(encodeBuffer.encodeBufferEnd, encodeBuffer.literalsEnd) = encoder.process(inputBegin, inputEnd, encodeBuffer.buffer.data(), encodeBuffer.literalsEnd);
  encoder = {}; // free memory

  std::vector<uint8_t> dict(512 * 1024 * 1024, 0); // 512 MB should be enough for everyone
  auto dictEnd = compressRenormedDictionary(renormedHistogram, dict.data());

  auto decoder = makeDecoder<>::fromRenormed(renormedHistogram);
  renormedHistogram = {}; // free memory

  decoder.process(encodeBuffer.encodeBufferEnd, outputBegin, inputSize, encoder.getNStreams(), encodeBuffer.literalsEnd);

  // Verify correctness
  auto [checkBegin, checkEnd] = rbm::getIter<field_V>(static_cast<const rbm::TPCData&>(dataOut));

  if (!(std::equal(inputBegin, inputEnd, checkBegin))) {
    LOGP(error, "Missmatch between original and decoded Message");
  }
  LOGP(info, "finished: {}", fieldTrait.name);

  decoder = {}; // free memory

  // Compression Properties
  // ##########################
  writer.Key("Compression");
  writer.StartObject();
  writer.Key("EncodeBufferSize");
  writer.Uint64(std::distance(encodeBuffer.buffer.data(), encodeBuffer.encodeBufferEnd) * sizeof(uint32_t));
  writer.Key("LiteralSize");
  writer.Uint64(
    [&]() -> uint64_t {
      const size_t extent = std::distance(encodeBuffer.literals.data(), encodeBuffer.literalsEnd);
      if (extent == 0) {
        return 0;
      }
      const auto [minIt, maxIt] = std::minmax_element(encodeBuffer.literals.data(), encodeBuffer.literalsEnd);
      const size_t symbolRangeBits = utils::getRangeBits(*minIt, *maxIt);
      auto size = computePackingBufferSize<>(extent, symbolRangeBits);
      return size;
    }());
  writer.Key("DictSize");
  writer.Uint64(std::distance(dict.data(), dictEnd));
  writer.EndObject(); // Compression

  writer.EndObject(); // field
  writer.Flush();
};

template <rbm::TPCField... fields_V>
void run(const rbm::TPCData& dataIn, rbm::TPCData& dataOut, bmp::mp_list_v<fields_V...>)
{
  (ransEncodeDecode<fields_V>(dataIn, dataOut), ...);
};

static void writeCTFTimingInfo(const fs::path& path)
{
  auto filename = path.filename().string();
  // three groups of digits separated by underscores ctf_<RunNumber>_<FirstTForbit>_<TFCounter>_TPC.root
  std::regex pattern(R"(ctf_(\d+)_(\d+)_(\d+)_TPC\.root)");
  std::smatch matches;

  if (std::regex_search(filename, matches, pattern)) {
    uint64_t runNumber = std::stoull(matches[1]);
    uint64_t firstTForbit = std::stoull(matches[2]);
    uint64_t tfCounter = std::stoull(matches[3]);

    writer.Key("CTFTimingInfo");
    writer.StartObject();
    writer.Key("RunNumber");
    writer.Uint64(runNumber);
    writer.Key("FirstTForbit");
    writer.Uint64(firstTForbit);
    writer.Key("TFCounter");
    writer.Uint64(tfCounter);
    writer.EndObject();
  } else {
    LOGP(warning, "Could not extract CTF timing information from filename {} ", filename);
  }
};

int main(int argc, char* argv[])
{
  static constexpr const std::string HELP = "help";
  static constexpr const std::string IN = "in";
  static constexpr const std::string OUT = "out";
  static constexpr const std::string LOG_SEVERITY = "log-severity";
  static constexpr const std::string STRICT = "strict";
  static constexpr const std::string FORCE = "force";

  bool strictMode = {false};
  bool forceOverwrite{false};
  fs::path outputDir{"."};
  std::string logSeverity{"info"};
  std::vector<fs::path> inputFiles{};

  bpo::options_description options("Allowed options");
  // clang-format off
    options.add_options()
        ((HELP + ",h").c_str(), "print usage message")
        ((OUT + ",o").c_str(), bpo::value<fs::path>(&outputDir), "output directory for json files, default: working directory")
        ((LOG_SEVERITY + ",l").c_str(), bpo::value<std::string>(&logSeverity), "severity of FairLogger")
        ((STRICT + ",s").c_str(), bpo::bool_switch(&strictMode), "abort if input files are missing (default: skip missing)")
        ((FORCE + ",f").c_str(), bpo::bool_switch(&forceOverwrite), "force overwrite of existing output files");
  // clang-format on

  bpo::options_description positional("Positional arguments");
  // clang-format off
  positional.add_options()
      (IN.c_str(), bpo::value<std::vector<fs::path>>(&inputFiles), "input files");
  // clang-format on

  bpo::options_description all_options{};
  all_options.add(options).add(positional);
  bpo::positional_options_description p;
  p.add(IN.c_str(), -1);

  bpo::variables_map vm;
  try {
    bpo::store(bpo::command_line_parser(argc, argv).options(all_options).positional(p).run(), vm);
    bpo::notify(vm);
  } catch (const bpo::error& ex) {
    std::cerr << fmt::format("Command Line Error: {}\n", ex.what());
    return 1;
  }

  auto print_help = [&]() {
    std::cout << "Usage: o2-bench-rans-TPCCompressionMetrics [-o <out_dir>] [-l <log>] [-f] [-s] [-h] <file1> [file2...]\n\n";
    std::cout << options << "\n";
  };

  if (vm.count(HELP)) {
    print_help();
    return 0;
  }

  if (inputFiles.empty()) {
    std::cerr << "No input files specified.\n";
    print_help();
    return 1;
  }

  if (!fs::exists(outputDir)) {
    LOGP(warning, "Output directory {} does not exist. Creating it.", outputDir.c_str());
    fs::create_directories(outputDir);
  }

  if (vm.count(LOG_SEVERITY)) {
    fair::Logger::SetConsoleSeverity(vm[LOG_SEVERITY].as<std::string>().c_str());
  }

  for (const auto& inFilePath : inputFiles) {
    if (!fs::exists(inFilePath)) {
      if (strictMode) {
        LOGP(error, "Input file {} does not exist. Aborting due to strict mode.", inFilePath.c_str());
        return 1;
      } else {
        LOGP(warning, "Input file {} does not exist. Skipping.", inFilePath.c_str());
        continue;
      }
    }

    fs::path outFile = outputDir / inFilePath.filename();
    outFile.replace_extension(".json");
    if (fs::exists(outFile) && !forceOverwrite) {
      LOGP(warning, "Output file {} already exists. Use -f to force overwrite. Skipping.", outFile.c_str());
      continue;
    }

    std::ofstream of{outFile};
    if (!of) {
      std::runtime_error(fmt::format("could not open output file at path {}", outFile.c_str()));
    }

    rapidjson::OStreamWrapper stream{of};
    writer.Reset(stream);
    writer.StartObject();

    LOGP(info, "Processing input file: {}", inFilePath.c_str());

    rbm::TPCData tpcDataIn = rbm::loadTPCData(inFilePath);
    rbm::TPCData dataOut{
      tpcDataIn.nTracks,
      tpcDataIn.nAttachedClusters,
      tpcDataIn.nUnattachedClusters,
      tpcDataIn.nAttachedClustersReduced,
      tpcDataIn.nTriggers};

    LOGP(info, "loaded Compressed Clusters from {}", inFilePath.c_str());
    writeCTFTimingInfo(inFilePath);
    LOG(info) << "Starting rANS TPC Encode Decode";
    writer.Key("Blocks");
    writer.StartObject();
    run(tpcDataIn, dataOut, rbm::traits::TPCFieldsAll{});
    writer.EndObject();
    LOG(info) << "Finished rANS TPC Encode Decode";
    writer.EndObject();
    stream.Flush();
    of.close();
  }
};
