#include <vector>
#include <cstring>
#include <random>
#include <algorithm>
#include <execution>
#include <iterator>

#include <benchmark/benchmark.h>

#include "rANS/rans.h"
#include "rANS/SIMDDecoder.h"
#include "rANS/typetraits.h"
#include "rANS/renorm.h"

#ifdef ENABLE_VTUNE_PROFILER
#include <ittnotify.h>
#endif

#include "helpers.h"

using namespace o2::rans;

using ransCoder_type = uint64_t;
using ransStream_type = uint32_t;

template <typename source_T>
using decoder_type = SIMDDecoder<ransCoder_type, ransStream_type, source_T, NStreams, 1>;

inline constexpr size_t MessageSize = 1ull << 22;

template <typename source_T>
class SourceMessageProxy
{
 public:
  SourceMessageProxy(size_t messageSize)
  {
    if (mSourceMessage.empty()) {
      std::mt19937 mt(0); // same seed we want always the same distrubution of random numbers;
      const size_t draws = std::min(1ul << 20, static_cast<size_t>(std::numeric_limits<source_T>::max()));
      const double probability = 0.5;
      std::binomial_distribution<source_T> dist(draws, probability);
      const size_t sourceSize = messageSize / sizeof(source_T) + 1;
      mSourceMessage.resize(sourceSize);
      std::generate(std::execution::par_unseq, mSourceMessage.begin(), mSourceMessage.end(), [&dist, &mt]() { return dist(mt); });
    }
  }

  const auto& get() const { return mSourceMessage; };

 private:
  std::vector<source_T> mSourceMessage{};
};

inline const SourceMessageProxy<uint8_t> sourceMessage8{MessageSize};
inline const SourceMessageProxy<uint16_t> sourceMessage16{MessageSize};
inline const SourceMessageProxy<uint32_t> sourceMessage32{MessageSize};

template <typename T>
const auto& getMessage()
{
  if constexpr (std::is_same_v<uint8_t, T>) {
    return sourceMessage8.get();
  } else if constexpr (std::is_same_v<uint16_t, T>) {
    return sourceMessage16.get();
  } else {
    return sourceMessage32.get();
  }
};

template <typename source_T, ContainerTag containerTag_V, CoderTag coderTag_V>
void ransCompressionBenchmark(benchmark::State& st)
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

  const auto& inputData = getMessage<source_type>();
  EncodeBuffer<source_type> encodeBuffer{inputData.size()};
  DecodeBuffer<source_type> decodeBuffer{inputData.size()};

  auto frequencyTable = makeFrequencyTable<containerTag>::fromSamples(gsl::span<const source_type>(inputData));
  auto renormedFrequencyTable = renormCutoffIncompressible<>(frequencyTable, 20, 10);
  auto encoder = makeEncoder<coderTag_V>::fromRenormed(renormedFrequencyTable);
  auto encodeBufferEnd = encodeBuffer.encodeIter;

#ifdef ENABLE_VTUNE_PROFILER
  __itt_resume();
#endif
  for (auto _ : st) {
    benchmark::DoNotOptimize(encodeBufferEnd = encoder.process(inputData.data(), inputData.data() + inputData.size(), encodeBuffer.encodeIter));
  }
#ifdef ENABLE_VTUNE_PROFILER
  __itt_pause();
#endif

  auto decoderFrequencyTable = makeFrequencyTableFromSamples(inputData.begin(), inputData.end());
  auto decoderRenormed = renormCutoffIncompressible<>(decoderFrequencyTable, renormedFrequencyTable.getRenormingBits(), 10);
  decoder_type<source_type> decoder{decoderRenormed};
  decoder.process(encodeBufferEnd, decodeBuffer.buffer.data(), inputData.size());
  if (!(decodeBuffer == inputData)) {
    st.SkipWithError("Missmatch between encoded and decoded Message");
  }

  st.SetItemsProcessed(static_cast<int64_t>(inputData.size()) * static_cast<int64_t>(st.iterations()));
  st.SetBytesProcessed(static_cast<int64_t>(inputData.size()) * sizeof(source_type) * static_cast<int64_t>(st.iterations()));
};

BENCHMARK(ransCompressionBenchmark<uint8_t, ContainerTag::Dynamic, CoderTag::Compat>);
BENCHMARK(ransCompressionBenchmark<uint16_t, ContainerTag::Dynamic, CoderTag::Compat>);
BENCHMARK(ransCompressionBenchmark<uint32_t, ContainerTag::Dynamic, CoderTag::Compat>);

BENCHMARK(ransCompressionBenchmark<uint8_t, ContainerTag::Static, CoderTag::Compat>);
BENCHMARK(ransCompressionBenchmark<uint16_t, ContainerTag::Static, CoderTag::Compat>);

BENCHMARK(ransCompressionBenchmark<uint8_t, ContainerTag::Hash, CoderTag::Compat>);
BENCHMARK(ransCompressionBenchmark<uint16_t, ContainerTag::Hash, CoderTag::Compat>);
BENCHMARK(ransCompressionBenchmark<uint32_t, ContainerTag::Hash, CoderTag::Compat>);

//########################################################################################

BENCHMARK(ransCompressionBenchmark<uint8_t, ContainerTag::Dynamic, CoderTag::SingleStream>);
BENCHMARK(ransCompressionBenchmark<uint16_t, ContainerTag::Dynamic, CoderTag::SingleStream>);
BENCHMARK(ransCompressionBenchmark<uint32_t, ContainerTag::Dynamic, CoderTag::SingleStream>);

BENCHMARK(ransCompressionBenchmark<uint8_t, ContainerTag::Static, CoderTag::SingleStream>);
BENCHMARK(ransCompressionBenchmark<uint16_t, ContainerTag::Static, CoderTag::SingleStream>);

BENCHMARK(ransCompressionBenchmark<uint8_t, ContainerTag::Hash, CoderTag::SingleStream>);
BENCHMARK(ransCompressionBenchmark<uint16_t, ContainerTag::Hash, CoderTag::SingleStream>);
BENCHMARK(ransCompressionBenchmark<uint32_t, ContainerTag::Hash, CoderTag::SingleStream>);

//########################################################################################

BENCHMARK(ransCompressionBenchmark<uint8_t, ContainerTag::Dynamic, CoderTag::SSE>);
BENCHMARK(ransCompressionBenchmark<uint16_t, ContainerTag::Dynamic, CoderTag::SSE>);
BENCHMARK(ransCompressionBenchmark<uint32_t, ContainerTag::Dynamic, CoderTag::SSE>);

BENCHMARK(ransCompressionBenchmark<uint8_t, ContainerTag::Static, CoderTag::SSE>);
BENCHMARK(ransCompressionBenchmark<uint16_t, ContainerTag::Static, CoderTag::SSE>);

BENCHMARK(ransCompressionBenchmark<uint8_t, ContainerTag::Hash, CoderTag::SSE>);
BENCHMARK(ransCompressionBenchmark<uint16_t, ContainerTag::Hash, CoderTag::SSE>);
BENCHMARK(ransCompressionBenchmark<uint32_t, ContainerTag::Hash, CoderTag::SSE>);

//########################################################################################

BENCHMARK(ransCompressionBenchmark<uint8_t, ContainerTag::Dynamic, CoderTag::AVX2>);
BENCHMARK(ransCompressionBenchmark<uint16_t, ContainerTag::Dynamic, CoderTag::AVX2>);
BENCHMARK(ransCompressionBenchmark<uint32_t, ContainerTag::Dynamic, CoderTag::AVX2>);

BENCHMARK(ransCompressionBenchmark<uint8_t, ContainerTag::Static, CoderTag::AVX2>);
BENCHMARK(ransCompressionBenchmark<uint16_t, ContainerTag::Static, CoderTag::AVX2>);

BENCHMARK(ransCompressionBenchmark<uint8_t, ContainerTag::Hash, CoderTag::AVX2>);
BENCHMARK(ransCompressionBenchmark<uint16_t, ContainerTag::Hash, CoderTag::AVX2>);
BENCHMARK(ransCompressionBenchmark<uint32_t, ContainerTag::Hash, CoderTag::AVX2>);

BENCHMARK_MAIN();