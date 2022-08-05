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

/// \file CombinedIterator.h
/// \brief
/// \author michael.lettrich@cern.ch

#include <type_traits>
#include <cstdint>
#include <stdexcept>

#include <rapidjson/writer.h>
#include "rANS/definitions.h"
#include "rANS/typetraits.h"
#include "rANS/utils/HistogramView.h"

#include "rANS/Packer.h"

namespace o2
{
namespace rans
{
namespace utils
{

enum class SourceType : uint8_t { Char,
                                  Uint8,
                                  Int8,
                                  Uint16,
                                  Int16,
                                  Uint32,
                                  Int32 };

template <typename T>
struct toSourceType;

template <>
struct toSourceType<char> : public std::integral_constant<SourceType, SourceType::Char> {
};
template <>
struct toSourceType<uint8_t> : public std::integral_constant<SourceType, SourceType::Uint8> {
};
template <>
struct toSourceType<int8_t> : public std::integral_constant<SourceType, SourceType::Int8> {
};
template <>
struct toSourceType<uint16_t> : public std::integral_constant<SourceType, SourceType::Uint16> {
};
template <>
struct toSourceType<int16_t> : public std::integral_constant<SourceType, SourceType::Int16> {
};
template <>
struct toSourceType<uint32_t> : public std::integral_constant<SourceType, SourceType::Uint32> {
};
template <>
struct toSourceType<int32_t> : public std::integral_constant<SourceType, SourceType::Int32> {
};

template <typename T>
inline constexpr SourceType toSourceType_v = toSourceType<T>::value;

template <SourceType sourceType_V>
struct getSourceType;

template <>
struct getSourceType<SourceType::Char> : public std::integral_constant<char, 0> {
};
template <>
struct getSourceType<SourceType::Uint8> : public std::integral_constant<uint8_t, 0> {
};
template <>
struct getSourceType<SourceType::Int8> : public std::integral_constant<int8_t, 0> {
};
template <>
struct getSourceType<SourceType::Uint16> : public std::integral_constant<uint16_t, 0> {
};
template <>
struct getSourceType<SourceType::Int16> : public std::integral_constant<int16_t, 0> {
};
template <>
struct getSourceType<SourceType::Uint32> : public std::integral_constant<uint32_t, 0> {
};
template <>
struct getSourceType<SourceType::Int32> : public std::integral_constant<int32_t, 0> {
};

template <SourceType sourceType_V>
using getSourceType_t = typename getSourceType<sourceType_V>::value_type;

template <typename source_T>
struct FlatContainer {
  std::vector<source_T> index{};
  std::vector<count_t> count{};
  count_t incompressibleCount{};
};

template <typename container_T>
count_t getFrequency(const container_T& container, typename container_T::source_type sourceSymbol)
{
  const auto& ret = container[sourceSymbol];
  if constexpr (isSymbolTable_v<container_T>) {
    return container.isEscapeSymbol(ret) ? 0 : ret.getFrequency();
  } else {
    return ret;
  }
};

template <typename container_T>
count_t getFrequency(const container_T& container, typename container_T::const_iterator iter)
{
  const auto& ret = *iter;
  if constexpr (isSymbolTable_v<container_T>) {
    return container.isEscapeSymbol(ret) ? 0 : ret.getFrequency();
  } else {
    return ret;
  }
};

template <typename container_T>
count_t getIncompressibleFrequency(const container_T& container)
{
  if constexpr (isSymbolTable_v<container_T>) {
    container.getEscapeSymbol().getFrequency();
  } else {
    return 0;
  }
};

template <typename container_T>
auto flattenContainer(const container_T& container) -> FlatContainer<typename container_T::source_type>
{
  if (container.empty()) {
    return {};
  }

  FlatContainer<typename container_T::source_type> flattened;

  flattened.index.reserve(container.size());
  flattened.count.reserve(container.size());

  uint32_t index = container.getOffset();
  for (auto iter = container.begin(); iter != container.end(); ++iter) {
    const count_t count = getFrequency(container, iter);
    if (count > 0) {
      flattened.index.push_back(index);
      flattened.count.push_back(count);
    }
    ++index;
  };

  return flattened;
};

// template <typename container_T>
// class makeFrequencyContainer
// {
//   using source_type = typename container_T::source_type;
//   using index_type = typename container_T::index_type;
//   using container_type = typename container_T::container_type;

//  public:
//   static container_T fromFlatContainer(FlatContainer<source_type>& flatContainer)
//   {
//     constexpr ContainerTag tag = getContainerTag_v<container_T>;
//     const size_t size = flatContainer.index.back() - index.front() + 1;
//     std::vector<source_type> container(size, 0);

//     container_type container();
//   }

//  private:
//   static container_type makeBaseContainer(

//   )
// };

// template <typename source_T, ContainerTag tag_V, std::enable_if_t<tag_V == ContainerTag::Static, bool> = true>
// decltype(auto) makeContainerImpl(FlatContainer<source_T>)
// {
//   using source_type = source_T;
//   using frequencyTable_type = typename ContainerTraits<tag_V>::frequencyTable_type<source_T>;
//   using container_type = typename frequencyTable_type::container_type;
//   using index_type = typename frequencyTable_type::index_type;

//   container_type container;

// }

template <typename container_T, typename jsonBuffer_T>
void toJSON(const container_T& container, rapidjson::Writer<jsonBuffer_T>& writer)
{
  using source_type = typename container_T::source_type;

  auto writeArray = [&](const auto& array, const std::string& name) {
    writer.Key(name.c_str());
    writer.StartArray();
    for (auto& i : array) {
      writer.Int64(i);
    }
    writer.EndArray();
  };

  FlatContainer<source_type> flattened = flattenContainer(container);
  flattened.incompressibleCount = getIncompressibleFrequency(container);
  const uint32_t extent = flattened.index.size() > 0 ? flattened.index.back() - flattened.index.front() + 1 : 0;

  writer.StartObject();
  writer.Key("SourceType");
  writer.Uint(static_cast<uint32_t>(toSourceType_v<source_type>));
  writer.Key("Extent");
  writer.Uint(extent);
  writeArray(flattened.index, "Index");
  writeArray(flattened.count, "Value");
  writer.EndObject();
};

template <typename container_T, typename dest_IT>
dest_IT toCompressedBinary(const container_T& container, dest_IT dstBufferBegin)
{
  using source_type = typename container_T::source_type;
  FlatContainer<source_type> flattened = flattenContainer(container);
  flattened.incompressibleCount = getIncompressibleFrequency(container);
  const uint32_t extent = flattened.index.size() > 0 ? flattened.index.back() - flattened.index.front() + 1 : 0;

  uint32_t bitOffset = 0;
  uint64_t* iter = reinterpret_cast<uint64_t*>(dstBufferBegin);
  //iterate backwards, store all but 0;
  for (size_t i = flattened.index.size(); i-- > 1ull;) {
    source_type indexDelta = flattened.index[i] - flattened.index[i - 1];
    assert(flattened.count[i] > 0);
    assert(indexDelta > 0);
    eliasDeltaEncode(iter, bitOffset, flattened.count[i]);
    eliasDeltaEncode(iter, bitOffset, indexDelta);
  }
  eliasDeltaEncode(iter, bitOffset, flattened.count[0]);
  eliasDeltaEncode(iter, bitOffset, static_cast<uint32_t>(flattened.index[0])); //TODO(milettri): zigzag encode
  pack(iter, bitOffset, extent, 32);
  eliasDeltaEncode(iter, bitOffset, (flattened.incompressibleCount + 1));
  eliasDeltaEncode(iter, bitOffset, 1);

  return reinterpret_cast<dest_IT>(++iter);
}

// template <typename source_T>
// dest_IT toCompressedBinary(const container_T&, dest_IT dstBufferBegin)
// {
//   using source_type = typename container_T::source_type;
//   FlatContainer<source_type> flattened = flattenContainer(container);
//   flattened.incompressibleCount = getIncompressibleFrequency(container);
//   const uint32_t extent = flattened.index.size() > 0 ? flattened.index.back() - flattened.index.front() + 1 : 0;

//   uint32_t bitOffset = 0;
//   uint64_t* iter = reinterpret_cast<uint64_t*>(dstBufferBegin);
//   //iterate backwards, store all but 0;
//   for (size_t i = flattened.index.size(); flattened-- > 1;) {
//     source_type indexDelta = flattened.index[i] - flattened.index[i - 1];
//     assert(flattened.count[i] > 0);
//     assert(indexDelta > 0);
//     eliasDeltaEncode(iter, bitOffset, flattened.count[i]);
//     eliasDeltaEncode(iter, bitOffset, indexDelta);
//   }
//   eliasDeltaEncode(iter, bitOffset, flattened.count[0]);
//   eliasDeltaEncode(iter, bitOffset, static_cast<uint32_t>(flattened.index[0])); //TODO(milettri): zigzag encode
//   pack(iter, bitOffset, extent, 32);
//   eliasDeltaEncode(iter, bitOffset, (flattened.incompressibleCount + 1));
//   eliasDeltaEncode(iter, bitOffset, static_cast<uint32_t>(toSourceType_v<source_type>));
//   eliasDeltaEncode(iter, bitOffset, 1);

//   return ++iter;
// }

} // namespace utils
} // namespace rans
} // namespace o2