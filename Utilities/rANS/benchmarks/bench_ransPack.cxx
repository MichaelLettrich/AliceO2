#include <vector>
#include <cstring>
#include <random>
#include <algorithm>
#include <execution>
#include <iterator>

#include <benchmark/benchmark.h>

#include "rANS/internal/pack/pack.h"

#ifdef ENABLE_VTUNE_PROFILER
#include <ittnotify.h>
#endif

#include "helpers.h"

using namespace o2::rans;

inline constexpr size_t MessageSize = 1ull << 22;

using source_type = uint32_t;

template <typename source_T>
std::vector<source_T> makeRandomUniformVector(size_t nelems, source_T min = std::numeric_limits<source_T>::max(), source_T max = std::numeric_limits<source_T>::max())
{
  std::vector<source_T> result(nelems, 0);
  std::mt19937 mt(0); // same seed we want always the same distrubution of random numbers;
  std::uniform_int_distribution<source_T> dist(min, max);

  std::generate(std::execution::par_unseq, result.begin(), result.end(), [&dist, &mt]() { return dist(mt); });
  return result;
};

static void copyBenchmark(benchmark::State& state)
{
  std::vector<source_type> src = makeRandomUniformVector<source_type>(MessageSize);
  std::vector<source_type> dst(MessageSize, 0);
  for (auto _ : state) {
    std::copy(src.begin(), src.end(), dst.begin());
  };
  state.SetItemsProcessed(src.size() * state.iterations());
  state.SetBytesProcessed(src.size() * sizeof(source_type) * state.iterations());
};

static void packingBenchmark(benchmark::State& state)
{
  size_t packingBits = state.range(0);

  std::vector<source_type> src = makeRandomUniformVector<source_type>(MessageSize, 0, internal::pow2(packingBits) - 1);
  std::vector<uint32_t> dst(MessageSize, 0);
  for (auto _ : state) {
    internal::BitPtr iter{dst.data()};
    for (auto i : src) {
      iter = internal::pack(iter, i, packingBits);
    }
  };
  state.SetItemsProcessed(src.size() * state.iterations());
  state.SetBytesProcessed(src.size() * sizeof(uint32_t) * state.iterations());

  std::vector<uint32_t> unpacked(MessageSize, 0);

  internal::BitPtr iter{dst.data()};
  for (size_t i = 0; i < src.size(); ++i) {
    unpacked[i] = internal::unpack<uint32_t>(iter, packingBits);
    iter += packingBits;
  }
  if (!std::equal(unpacked.begin(), unpacked.end(), src.begin())) {
    state.SkipWithError("error in packing");
  }
};

static void fastPackBenchmark(benchmark::State& state)
{
  size_t packingBits = state.range(0);

  std::vector<source_type> src = makeRandomUniformVector<source_type>(MessageSize, 0, internal::pow2(packingBits) - 1);
  std::vector<uint32_t> dst(MessageSize, 0);
#ifdef ENABLE_VTUNE_PROFILER
  __itt_resume();
#endif
  for (auto _ : state) {
    pack(src.data(), src.size(), dst.data(), packingBits, 0u);
  };
#ifdef ENABLE_VTUNE_PROFILER
  __itt_pause();
#endif
  state.SetItemsProcessed(src.size() * state.iterations());
  state.SetBytesProcessed(src.size() * sizeof(uint32_t) * state.iterations());

  std::vector<uint32_t> unpacked(MessageSize, 0);

  internal::BitPtr iter{dst.data()};
  for (size_t i = 0; i < src.size(); ++i) {
    unpacked[i] = internal::unpack<uint32_t>(iter, packingBits);
    // LOGP(info, "[{}]{:0" + std::to_string(packingBits) + "b}", i, unpacked[i]);
    iter += packingBits;
  }
  size_t i = 0;
  if (!std::equal(unpacked.begin(), unpacked.end(), src.begin(), [&i](auto a, auto b) -> bool {
        if (a != b) {
          LOGP(info, "[{}]{:0x}!={:0x}", i++, a, b);
        }
        return a == b;
      })) {
    state.SkipWithError("error in packing");
  }
};

BENCHMARK(copyBenchmark);
BENCHMARK(packingBenchmark)->DenseRange(1, 32, 1);
BENCHMARK(fastPackBenchmark)->DenseRange(1, 32, 1);
BENCHMARK_MAIN();
