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

#include "rANS/rans.h"
#include "rANS/LiteralSIMDEncoder.h"
#include "rANS/LiteralSIMDDecoder.h"
#include "rANS/typetraits.h"
#include "rANS/renorm.h"

#include "helpers.h"

#ifdef ENABLE_VTUNE_PROFILER
#include <ittnotify.h>
#endif

namespace bpo = boost::program_options;
using namespace o2::rans;

using ransCoder_type = uint64_t;
using ransStream_type = uint32_t;

template <typename source_T>
using decoder_type = LiteralSIMDDecoder<ransCoder_type, ransStream_type, source_T, NStreams, 1>;

// using container_types = boost::mp11::mp_list<std::integral_constant<ContainerTag, ContainerTag::Dynamic>,
//                                              std::integral_constant<ContainerTag, ContainerTag::Static>,
//                                              std::integral_constant<ContainerTag, ContainerTag::Hash>>;

using container_types = boost::mp11::mp_list<std::integral_constant<ContainerTag, ContainerTag::Dynamic>>;

// using coder_types = boost::mp11::mp_list<std::integral_constant<CoderTag, CoderTag::Compat>,
//                                          std::integral_constant<CoderTag, CoderTag::SingleStream>,
//                                          std::integral_constant<CoderTag, CoderTag::SSE>,
//                                          std::integral_constant<CoderTag, CoderTag::AVX2>>;

using coder_types = boost::mp11::mp_list<std::integral_constant<CoderTag, CoderTag::SingleStream>>;

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
      throw std::runtime_error("Invalid");
      break;
  };
};

std::string toString(ContainerTag tag)
{
  switch (tag) {
    case ContainerTag::Dynamic:
      return {"Dynamic"};
      break;
    case ContainerTag::Static:
      return {"Static"};
      break;
    case ContainerTag::Hash:
      return {"Hash"};
      break;
    default:
      throw std::runtime_error("Invalid");
      break;
  };
};

std::string makeEncoderTitle(ContainerTag containerTag, CoderTag coderTag)
{
  return toString(containerTag) + toString(coderTag) + "Encoder";
};

template <typename source_T, ContainerTag containerTag_V, CoderTag coderTag_V>
void ransEncodeDecode(const std::string& name, const std::vector<source_T>& inputData, rapidjson::Writer<rapidjson::OStreamWrapper>& writer)
{
  using source_type = source_T;
  constexpr ContainerTag containerTag = [&]() {
    if (containerTag_V != ContainerTag::Static) {
      return containerTag_V;
    } else {
      if (sizeof(source_T) < 4) {
        return ContainerTag::Static;
      } else {
        return ContainerTag::Dynamic;
      }
    }
  }();
  internal::RANSTimer timer{};
  writer.Key(name.c_str());
  writer.StartObject();

  EncodeBuffer<source_type> encodeBuffer{inputData.size()};
  encodeBuffer.literals.resize(inputData.size(), 0);
  encodeBuffer.literalsEnd = encodeBuffer.literals.data();
  DecodeBuffer<source_type> decodeBuffer{inputData.size()};

  writer.Key("Timing");
  writer.StartObject();

  LOGP(info, "processing: {} (nItems: {}, size: {} MiB)", name, inputData.size(), inputData.size() * sizeof(source_type) / 1024.0 / 1024.0);
  timer.start();
  auto frequencyTable = makeFrequencyTable<containerTag>::fromSamples(gsl::span<const source_type>(inputData));

  timer.stop();
  writer.Key("FrequencyTable");
  writer.Double(timer.getDurationMS());
  LOGP(info, "Built Frequency Table in {} ms", timer.getDurationMS());

  timer.start();
  auto renormedFrequencyTable = renormCutoffIncompressible<>(frequencyTable);

  timer.stop();
  writer.Key("Renorming");
  writer.Double(timer.getDurationMS());
  LOGP(info, "Renormed Frequency Table in {} ms", timer.getDurationMS());
  timer.start();
  auto encoder = makeEncoder<coderTag_V>::fromRenormed(renormedFrequencyTable);
  timer.stop();
  writer.Key("Encoder");
  writer.Double(timer.getDurationMS());
  LOGP(info, "Built Encoder in {} ms", timer.getDurationMS());
  timer.start();
#ifdef ENABLE_VTUNE_PROFILER
  __itt_resume();
#endif
  std::tie(encodeBuffer.encodeBufferEnd, encodeBuffer.literalsEnd) = encoder.process(inputData.data(), inputData.data() + inputData.size(), encodeBuffer.buffer.data(), encodeBuffer.literalsEnd);
#ifdef ENABLE_VTUNE_PROFILER
  __itt_pause();
#endif
  timer.stop();
  writer.Key("Encoding");
  writer.Double(timer.getDurationMS());
  LOGP(info, "Encoded {} Bytes in {} ms", inputData.size() * sizeof(source_type), timer.getDurationMS());

  auto decoderFrequencyTable = makeFrequencyTableFromSamples(inputData.begin(), inputData.end());
  auto decoderRenormed = renormCutoffIncompressible<>(decoderFrequencyTable, renormedFrequencyTable.getRenormingBits());
  decoder_type<source_type> decoder{decoderRenormed};
  encodeBuffer.literals.resize(std::distance(encodeBuffer.literals.data(), encodeBuffer.literalsEnd));
  decoder.process(encodeBuffer.encodeBufferEnd, decodeBuffer.buffer.data(), inputData.size(), encodeBuffer.literals);
  if (!(decodeBuffer == inputData)) {
    LOGP(warning, "Missmatch between original and decoded Message");
  }
  LOG(info) << "finished: " << name;

  writer.EndObject(); // Timing

  // Frequency Table
  // ##########################
  writer.Key("FrequencyTable");
  writer.StartObject();
  writer.Key("nSamples");
  writer.Uint64(decoderFrequencyTable.getNumSamples());
  writer.Key("Min");
  writer.Int(decoderFrequencyTable.getMinSymbol());
  writer.Key("Max");
  writer.Int(decoderFrequencyTable.getMaxSymbol());
  writer.Key("alphabetRangeBits");
  writer.Int(decoderFrequencyTable.getAlphabetRangeBits());
  writer.Key("nUsedAlphabetSymbols");
  writer.Uint(decoderFrequencyTable.getNUsedAlphabetSymbols());
  writer.Key("IncompressibleFrequency");
  writer.Uint(decoderFrequencyTable.getIncompressibleSymbolFrequency());
  writer.EndObject(); // FrequencyTable

  // RescaledFrequencies
  //##########################
  writer.Key("RescaledFrequencies");
  writer.StartObject();
  writer.Key("nSamples");
  writer.Uint64(decoderRenormed.getNumSamples());
  writer.Key("Min");
  writer.Int(decoderRenormed.getMinSymbol());
  writer.Key("Max");
  writer.Int(decoderRenormed.getMaxSymbol());
  writer.Key("alphabetRangeBits");
  writer.Int(decoderRenormed.getAlphabetRangeBits());
  writer.Key("nUsedAlphabetSymbols");
  writer.Uint(decoderRenormed.getNUsedAlphabetSymbols());
  writer.Key("IncompressibleFrequency");
  writer.Uint(decoderRenormed.getIncompressibleSymbolFrequency());
  writer.Key("RenormingBits");
  writer.Uint(decoderRenormed.getRenormingBits());
  writer.EndObject(); // RescaledFrequencies

  // Message Properties
  //##########################
  writer.Key("Message");
  writer.StartObject();
  writer.Key("Size");
  writer.Uint64(inputData.size());
  writer.Key("SymbolSize");
  writer.Uint(sizeof(source_type));
  writer.Key("Entropy");
  writer.Double(computeEntropy(decoderFrequencyTable));
  writer.Key("ExpectedCodewordLength");
  writer.Double(computeExpectedCodewordLength<source_type>(decoderFrequencyTable, decoderRenormed));
  writer.EndObject(); // Message

  // Message Properties
  //##########################
  writer.Key("Compression");
  writer.StartObject();
  writer.Key("EncodeBufferSize");
  writer.Uint64(std::distance(encodeBuffer.buffer.data(), encodeBuffer.encodeBufferEnd) * sizeof(uint32_t));
  writer.Key("LiteralSize");
  writer.Uint64(encodeBuffer.literals.size() * sizeof(source_type));
  writer.EndObject(); // Compression

  writer.EndObject(); // Encode/Decode Run
};

template <ContainerTag containerTag_V, CoderTag coderTag_V>
void encodeTPC(const std::string& name, const TPCCompressedClusters& compressedClusters, bool mergeColumns, rapidjson::Writer<rapidjson::OStreamWrapper>& writer)
{
  writer.Key(name.c_str());
  writer.StartObject();
  ransEncodeDecode<uint16_t, containerTag_V, coderTag_V>("qTotA", compressedClusters.qTotA, writer);
  ransEncodeDecode<uint16_t, containerTag_V, coderTag_V>("qMaxA", compressedClusters.qMaxA, writer);
  ransEncodeDecode<uint8_t, containerTag_V, coderTag_V>("flagsA", compressedClusters.flagsA, writer);
  ransEncodeDecode<uint8_t, containerTag_V, coderTag_V>("rowDiffA", compressedClusters.rowDiffA, writer);
  ransEncodeDecode<uint8_t, containerTag_V, coderTag_V>("sliceLegDiffA", compressedClusters.sliceLegDiffA, writer);
  ransEncodeDecode<uint16_t, containerTag_V, coderTag_V>("padResA", compressedClusters.padResA, writer);
  ransEncodeDecode<uint32_t, containerTag_V, coderTag_V>("timeResA", compressedClusters.timeResA, writer);
  ransEncodeDecode<uint8_t, containerTag_V, coderTag_V>("sigmaPadA", compressedClusters.sigmaPadA, writer);
  ransEncodeDecode<uint8_t, containerTag_V, coderTag_V>("sigmaTimeA", compressedClusters.sigmaTimeA, writer);
  ransEncodeDecode<uint8_t, containerTag_V, coderTag_V>("qPtA", compressedClusters.qPtA, writer);
  ransEncodeDecode<uint8_t, containerTag_V, coderTag_V>("rowA", compressedClusters.rowA, writer);
  ransEncodeDecode<uint8_t, containerTag_V, coderTag_V>("sliceA", compressedClusters.sliceA, writer);
  // ransEncodeDecode<uint32_t, containerTag_V,coderTag_V,>("timeA", compressedClusters.timeA, writer);
  ransEncodeDecode<uint16_t, containerTag_V, coderTag_V>("padA", compressedClusters.padA, writer);
  ransEncodeDecode<uint16_t, containerTag_V, coderTag_V>("qTotU", compressedClusters.qTotU, writer);
  ransEncodeDecode<uint16_t, containerTag_V, coderTag_V>("qMaxU", compressedClusters.qMaxU, writer);
  ransEncodeDecode<uint8_t, containerTag_V, coderTag_V>("flagsU", compressedClusters.flagsU, writer);
  ransEncodeDecode<uint16_t, containerTag_V, coderTag_V>("padDiffU", compressedClusters.padDiffU, writer);
  ransEncodeDecode<uint32_t, containerTag_V, coderTag_V>("timeDiffU", compressedClusters.timeDiffU, writer);
  ransEncodeDecode<uint8_t, containerTag_V, coderTag_V>("sigmaPadU", compressedClusters.sigmaPadU, writer);
  ransEncodeDecode<uint8_t, containerTag_V, coderTag_V>("sigmaTimeU", compressedClusters.sigmaTimeU, writer);
  ransEncodeDecode<uint16_t, containerTag_V, coderTag_V>("nTrackClusters", compressedClusters.nTrackClusters, writer);
  ransEncodeDecode<uint32_t, containerTag_V, coderTag_V>("nSliceRowClusters", compressedClusters.nSliceRowClusters, writer);

  writer.EndObject();
};

using encoder_types = boost::mp11::mp_product<boost::mp11::mp_list, container_types, coder_types>;

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
  rapidjson::Writer<rapidjson::OStreamWrapper> writer{stream};
  writer.StartObject();

  TPCCompressedClusters compressedClusters = readFile(inFile);
  LOG(info) << "loaded Compressed Clusters from file";
  LOG(info) << "######################################################";
  boost::mp11::mp_for_each<encoder_types>([&](auto L) {
    using container_type = boost::mp11::mp_at_c<decltype(L), 0>;
    using coder_type = boost::mp11::mp_at_c<decltype(L), 1>;
    constexpr ContainerTag containerTag = container_type::value;
    constexpr CoderTag coderTag = coder_type::value;
    const std::string encoderTitle = makeEncoderTitle(containerTag, coderTag);

    LOGP(info, "start rANS {}/Decode", encoderTitle);
    encodeTPC<containerTag, coderTag>(encoderTitle, compressedClusters, false, writer);
    LOG(info) << "######################################################";
  });
  writer.EndObject();
  writer.Flush();
  of.close();
};