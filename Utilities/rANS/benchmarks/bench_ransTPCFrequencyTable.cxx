#include <vector>
#include <cstring>
#include <execution>
#include <iterator>
#include <iostream>

#ifdef ENABLE_VTUNE_PROFILER
#include <ittnotify.h>
#endif

#include <boost/program_options.hpp>
#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/ostreamwrapper.h>
#include <fairlogger/Logger.h>

#include "rANS/rans.h"
#include "rANS/LiteralSIMDEncoder.h"
#include "rANS/LiteralSIMDDecoder.h"
#include "rANS/StaticFrequencyTable.h"
#include "rANS/HashFrequencyTable.h"

namespace bpo = boost::program_options;

template <typename source_T>
double_t computeExpectedCodewordLength(const o2::rans::FrequencyTable& frequencies, const o2::rans::RenormedFrequencyTable& rescaled)
{

  using symbol_t = o2::rans::symbol_t;
  using count_t = o2::rans::count_t;
  double_t expectedCodewordLength = 0;
  count_t trueIncompressibleFrequency = frequencies.getIncompressibleSymbolFrequency();

  auto getRescaledFrequency = [&rescaled](symbol_t sourceSymbol) {
    if (sourceSymbol >= rescaled.getMinSymbol() && sourceSymbol <= rescaled.getMaxSymbol()) {
      return rescaled[sourceSymbol];
    } else {
      return static_cast<count_t>(0);
    }
  };

  // all "normal symbols"
  for (symbol_t sourceSymbol = frequencies.getMinSymbol(); sourceSymbol <= frequencies.getMaxSymbol(); ++sourceSymbol) {

    const count_t frequency = frequencies[sourceSymbol];
    if (frequency) {
      const count_t rescaledFrequency = getRescaledFrequency(sourceSymbol);

      const double_t trueProbability = static_cast<double_t>(frequency) / frequencies.getNumSamples();

      if (rescaledFrequency) {
        const double_t rescaledProbability = static_cast<double_t>(rescaledFrequency) / rescaled.getNumSamples();
        expectedCodewordLength -= trueProbability * std::log2(rescaledProbability);
      } else {
        trueIncompressibleFrequency += frequency;
      }
    }
  }
  // incompressibleSymbol:
  const double_t trueProbability = static_cast<double_t>(trueIncompressibleFrequency) / frequencies.getNumSamples();
  const double_t rescaledProbability = static_cast<double_t>(rescaled.getIncompressibleSymbolFrequency()) / rescaled.getNumSamples();

  expectedCodewordLength -= trueProbability * std::log2(rescaledProbability);
  expectedCodewordLength += trueProbability * std::log2(o2::rans::internal::toBits(sizeof(source_T)));

  return expectedCodewordLength;
};

struct TPCCompressedClusters {

  TPCCompressedClusters() = default;
  TPCCompressedClusters(size_t nTracks,
                        size_t nAttachedClusters,
                        size_t nUnattachedClusters,
                        size_t nAttachedClustersReduced) : nTracks(nTracks),
                                                           nAttachedClusters(nAttachedClusters),
                                                           nUnattachedClusters(nUnattachedClusters),
                                                           nAttachedClustersReduced(nAttachedClustersReduced)
  {
    qTotA.resize(this->nAttachedClusters);
    qMaxA.resize(this->nAttachedClusters);
    flagsA.resize(this->nAttachedClusters);
    rowDiffA.resize(this->nAttachedClustersReduced);
    sliceLegDiffA.resize(this->nAttachedClustersReduced);
    padResA.resize(this->nAttachedClustersReduced);
    timeResA.resize(this->nAttachedClustersReduced);
    sigmaPadA.resize(this->nAttachedClusters);
    sigmaTimeA.resize(this->nAttachedClusters);

    qPtA.resize(this->nTracks);
    rowA.resize(this->nTracks);
    sliceA.resize(this->nTracks);
    timeA.resize(this->nTracks);
    padA.resize(this->nTracks);

    qTotU.resize(this->nUnattachedClusters);
    qMaxU.resize(this->nUnattachedClusters);
    flagsU.resize(this->nUnattachedClusters);
    padDiffU.resize(this->nUnattachedClusters);
    timeDiffU.resize(this->nUnattachedClusters);
    sigmaPadU.resize(this->nUnattachedClusters);
    sigmaTimeU.resize(this->nUnattachedClusters);

    nTrackClusters.resize(this->nTracks);
    nSliceRowClusters.resize(this->nSliceRows);
  };

  size_t nTracks = 0;
  size_t nAttachedClusters = 0;
  size_t nUnattachedClusters = 0;
  size_t nAttachedClustersReduced = 0;
  size_t nSliceRows = 36 * 152;

  std::vector<uint16_t> qTotA{};
  std::vector<uint16_t> qMaxA{};
  std::vector<uint8_t> flagsA{};
  std::vector<uint8_t> rowDiffA{};
  std::vector<uint8_t> sliceLegDiffA{};
  std::vector<uint16_t> padResA{};
  std::vector<uint32_t> timeResA{};
  std::vector<uint8_t> sigmaPadA{};
  std::vector<uint8_t> sigmaTimeA{};

  std::vector<uint8_t> qPtA{};
  std::vector<uint8_t> rowA{};
  std::vector<uint8_t> sliceA{};
  std::vector<int32_t> timeA{};
  std::vector<uint16_t> padA{};

  std::vector<uint16_t> qTotU{};
  std::vector<uint16_t> qMaxU{};
  std::vector<uint8_t> flagsU{};
  std::vector<uint16_t> padDiffU{};
  std::vector<uint32_t> timeDiffU{};
  std::vector<uint8_t> sigmaPadU{};
  std::vector<uint8_t> sigmaTimeU{};

  std::vector<uint16_t> nTrackClusters{};
  std::vector<uint32_t> nSliceRowClusters{};
};

class TPCJsonHandler : public rapidjson::BaseReaderHandler<rapidjson::UTF8<>, TPCJsonHandler>
{
 private:
  class CurrentVector
  {
   public:
    CurrentVector() = default;

    template <typename T>
    CurrentVector(std::vector<T>& vec) : mVectorConcept{std::make_unique<VectorWrapper<T>>(vec)} {};

    inline void push_back(unsigned i) { mVectorConcept->push_back(i); };

    struct VectorConcept {
      virtual ~VectorConcept() = default;
      virtual void push_back(unsigned i) = 0;
    };

    template <typename T>
    struct VectorWrapper : VectorConcept {
      VectorWrapper(std::vector<T>& vector) : mVector(vector){};

      void push_back(unsigned i) override { mVector.push_back(static_cast<T>(i)); };

     private:
      std::vector<T>& mVector{};
    };

    std::unique_ptr<VectorConcept> mVectorConcept{};
  };

 public:
  bool Null() { return true; };
  bool Bool(bool b) { return true; };
  bool Int(int i) { return this->Uint(static_cast<unsigned>(i)); };
  bool Uint(unsigned i)
  {
    mCurrentVector.push_back(i);
    return true;
  };
  bool Int64(int64_t i) { return this->Uint(static_cast<unsigned>(i)); };
  bool Uint64(uint64_t i) { return this->Uint(static_cast<unsigned>(i)); };
  bool Double(double d) { return this->Uint(static_cast<unsigned>(d)); };
  bool RawNumber(const Ch* str, rapidjson::SizeType length, bool copy) { return true; };
  bool String(const Ch* str, rapidjson::SizeType length, bool copy) { return true; };
  bool StartObject() { return true; };
  bool Key(const Ch* str, rapidjson::SizeType length, bool copy)
  {
    if (str == std::string{"qTotA"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.qTotA};
    } else if (str == std::string{"qMaxA"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.qMaxA};
    } else if (str == std::string{"flagsA"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.flagsA};
    } else if (str == std::string{"rowDiffA"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.rowDiffA};
    } else if (str == std::string{"sliceLegDiffA"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.sliceLegDiffA};
    } else if (str == std::string{"padResA"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.padResA};
    } else if (str == std::string{"timeResA"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.timeResA};
    } else if (str == std::string{"sigmaPadA"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.sigmaPadA};
    } else if (str == std::string{"sigmaTimeA"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.sigmaTimeA};
    } else if (str == std::string{"qPtA"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.qPtA};
    } else if (str == std::string{"rowA"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.rowA};
    } else if (str == std::string{"sliceA"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.sliceA};
    } else if (str == std::string{"timeA"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.timeA};
    } else if (str == std::string{"padA"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.padA};
    } else if (str == std::string{"qTotU"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.qTotU};
    } else if (str == std::string{"qMaxU"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.qMaxU};
    } else if (str == std::string{"flagsU"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.flagsU};
    } else if (str == std::string{"padDiffU"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.padDiffU};
    } else if (str == std::string{"timeDiffU"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.timeDiffU};
    } else if (str == std::string{"sigmaPadU"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.sigmaPadU};
    } else if (str == std::string{"sigmaTimeU"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.sigmaTimeU};
    } else if (str == std::string{"nTrackClusters"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.nTrackClusters};
    } else if (str == std::string{"nSliceRowClusters"}) {
      LOGP(info, "parsing {}", str);
      mCurrentVector = CurrentVector{mCompressedClusters.nSliceRowClusters};
    } else {
      throw std::runtime_error(fmt::format("invalid key: {}", str));
      return false;
    }
    return true;
  };
  bool EndObject(rapidjson::SizeType memberCount)
  {
    LOGP(info, "parsed {} objects", memberCount);
    return true;
  };
  bool StartArray()
  {
    return true;
  };
  bool EndArray(rapidjson::SizeType elementCount)
  {
    LOGP(info, "parsed {} elements", elementCount);
    return true;
  };

  TPCCompressedClusters release() && { return std::move(mCompressedClusters); };

 private:
  TPCCompressedClusters mCompressedClusters;
  CurrentVector mCurrentVector;
};

TPCCompressedClusters readFile(const std::string& filename)
{
  TPCCompressedClusters compressedClusters;
  std::ifstream is(filename, std::ios_base::in);
  if (is) {
    rapidjson::IStreamWrapper isWrapper{is};
    // rapidjson::Document document;
    TPCJsonHandler handler;
    rapidjson::Reader reader;
    reader.Parse(isWrapper, handler);

    compressedClusters = std::move(handler).release();

    compressedClusters.nAttachedClusters = compressedClusters.qTotA.size();
    compressedClusters.nAttachedClustersReduced = compressedClusters.rowDiffA.size();
    compressedClusters.nUnattachedClusters = compressedClusters.qTotA.size();
    compressedClusters.nTracks = compressedClusters.nTrackClusters.size();
    is.close();
  } else {
    throw std::runtime_error(fmt::format("Could not open file {}", filename));
  }
  return compressedClusters;
};

template <typename source_T>
double_t buildDynamicFrequencyTable(const std::vector<source_T>& inputData)
{
  using namespace o2;
  rans::internal::RANSTimer timer{};
  timer.start();
  auto frequencyTable = rans::makeFrequencyTableFromSamples(inputData.begin(), inputData.end());
  timer.stop();
  return timer.getDurationMS();
};

template <typename source_T>
double_t buildBoundedDynamicFrequencyTable(const std::vector<source_T>& inputData)
{
  using namespace o2;
  rans::internal::RANSTimer timer{};

  timer.start();
  auto frequencyTable = rans::makeFrequencyTableFromSamples(inputData.begin(), inputData.end());
  timer.stop();
  return timer.getDurationMS();
};

template <typename source_T>
double_t buildStaticFrequencyTable(const std::vector<source_T>& inputData)
{
  using namespace o2;
  rans::internal::RANSTimer timer{};

  timer.start();
  rans::StaticFrequencyTable<source_T> frequencyTable{};
  frequencyTable.addSamples(gsl::make_span(inputData));
  timer.stop();
  return timer.getDurationMS();
};

template <typename source_T>
double_t buildHashFrequencyTable(const std::vector<source_T>& inputData)
{
  using namespace o2;
  rans::internal::RANSTimer timer{};

  timer.start();
  rans::HashFrequencyTable<source_T> frequencyTable{};
  frequencyTable.addSamples(gsl::make_span(inputData));
  timer.stop();
  return timer.getDurationMS();
};

template <typename source_T,
          bool dynamicEnable = false,
          bool boundedDynamicEnable = false,
          bool staticEnable = false,
          bool hashEnable = false>
void processColumn(const std::string& name, const std::vector<source_T>& inputData, rapidjson::Writer<rapidjson::OStreamWrapper>& writer)
{
  using namespace o2;

  writer.Key(name.c_str());
  writer.StartObject();

  writer.Key("Timing");
  writer.StartObject();

  LOGP(info, "processing: {} (nItems: {}, size: {} MiB)", name, inputData.size(), inputData.size() * sizeof(source_T) / 1024.0 / 1024.0);

  auto computeBandwidth = [&](double_t timerMS) {
    return (static_cast<double_t>(inputData.size() * sizeof(source_T)) / 1024 / 1024) / (timerMS / 1000);
  };

  writer.Key("DynamicFrequencyTable");
  if constexpr (dynamicEnable) {
    const double_t timeMs = buildDynamicFrequencyTable<source_T>(inputData);
    writer.Double(timeMs);
    LOGP(info, "Built DynamicFrequencyTable in {} ms ( {} MiB/s)", timeMs, computeBandwidth(timeMs));
  } else {
    writer.String("NaN");
  };

  writer.Key("DynamicBoundedFrequencyTable");
  if constexpr (boundedDynamicEnable) {
    const double_t timeMs = buildDynamicFrequencyTable<source_T>(inputData);
    writer.Double(timeMs);
    LOGP(info, "Built DynamicBoundedFrequencyTable in {} ms ( {} MiB/s)", timeMs, computeBandwidth(timeMs));
  } else {
    writer.String("NaN");
  };

  writer.Key("StaticFrequencyTable");
  if constexpr (staticEnable) {
    const double_t timeMs = buildStaticFrequencyTable<source_T>(inputData);
    writer.Double(timeMs);
    LOGP(info, "Built StaticFrequencyTable in {} ms ( {} MiB/s)", timeMs, computeBandwidth(timeMs));
  } else {
    writer.String("NaN");
  };

  writer.Key("HashFrequencyTable");
  if constexpr (hashEnable) {
    const double_t timeMs = buildHashFrequencyTable<source_T>(inputData);
    writer.Double(timeMs);
    LOGP(info, "Built HashFrequencyTable in {} ms ( {} MiB/s)", timeMs, computeBandwidth(timeMs));
  } else {
    writer.String("NaN");
  };

  writer.EndObject(); // Timing

  // Message Properties
  //##########################
  writer.Key("Message");
  writer.StartObject();
  writer.Key("Size");
  writer.Uint64(inputData.size());
  writer.Key("SymbolSize");
  writer.Uint(sizeof(source_T));
  writer.EndObject(); // Message

  writer.EndObject(); // Dataset
};

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

  processColumn<uint16_t, true, true, true, true>("qTotA", compressedClusters.qTotA, writer);
  processColumn<uint16_t, true, true, true, true>("qMaxA", compressedClusters.qMaxA, writer);
  processColumn<uint8_t, true, true, true, true>("flagsA", compressedClusters.flagsA, writer);
  processColumn<uint8_t, true, true, true, true>("rowDiffA", compressedClusters.rowDiffA, writer);
  processColumn<uint8_t, true, true, true, true>("sliceLegDiffA", compressedClusters.sliceLegDiffA, writer);
  processColumn<uint16_t, true, true, true, true>("padResA", compressedClusters.padResA, writer);
  processColumn<uint32_t, true, true, false, true>("timeResA", compressedClusters.timeResA, writer);
  processColumn<uint8_t, true, true, true, true>("sigmaPadA", compressedClusters.sigmaPadA, writer);
  processColumn<uint8_t, true, true, true, true>("sigmaTimeA", compressedClusters.sigmaTimeA, writer);
  processColumn<uint8_t, true, true, true, true>("qPtA", compressedClusters.qPtA, writer);
  processColumn<uint8_t, true, true, true, true>("rowA", compressedClusters.rowA, writer);
  processColumn<uint8_t, true, true, true, true>("sliceA", compressedClusters.sliceA, writer);
  processColumn<int32_t, false, false, false, true>("timeA", compressedClusters.timeA, writer);
  processColumn<uint16_t, true, true, true, true>("padA", compressedClusters.padA, writer);
  processColumn<uint16_t, true, true, true, true>("qTotU", compressedClusters.qTotU, writer);
  processColumn<uint16_t, true, true, true, true>("qMaxU", compressedClusters.qMaxU, writer);
  processColumn<uint8_t, true, true, true, true>("flagsU", compressedClusters.flagsU, writer);
  processColumn<uint16_t, true, true, true, true>("padDiffU", compressedClusters.padDiffU, writer);
  processColumn<uint32_t, true, true, false, true>("timeDiffU", compressedClusters.timeDiffU, writer);
  processColumn<uint8_t, true, true, true, true>("sigmaPadU", compressedClusters.sigmaPadU, writer);
  processColumn<uint8_t, true, true, true, true>("sigmaTimeU", compressedClusters.sigmaTimeU, writer);
  processColumn<uint16_t, true, true, true, true>("nTrackClusters", compressedClusters.nTrackClusters, writer);
  processColumn<uint32_t, true, true, false, true>("nSliceRowClusters", compressedClusters.nSliceRowClusters, writer);

  writer.EndObject();
  writer.Flush();
  of.close();
};