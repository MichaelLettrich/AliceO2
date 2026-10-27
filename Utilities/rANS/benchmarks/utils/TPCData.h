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

/// @file   helpers.h
/// @author Michael Lettrich
/// @brief  common functionality for rANS benchmarks.

#ifndef RANS_BENCHMARKS_UTILS_TPCDATA_H_
#define RANS_BENCHMARKS_UTILS_TPCDATA_H_

#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <tuple>
#include <utility>
#include <type_traits>
#include <filesystem>

#include <boost/mp11.hpp>

#include "rANS/iterator.h"

class TTree; // forward declaration to avoid including ROOT headers here

#define DEFINE_FIELD_TRAIT_SIMPLE(ENUM_VAL, NAME_STR, MEMBER_PTR)          \
  template <>                                                              \
  struct FieldTrait<ENUM_VAL> {                                            \
    using strategy_type = detail::IdentityFunctor;                         \
    static constexpr auto members = std::make_tuple(MEMBER_PTR);           \
    using value_type = detail::extract_inner_type_t<decltype(MEMBER_PTR)>; \
    static constexpr std::string_view name = NAME_STR;                     \
  }

#define DEFINE_FIELD_TRAIT_PACKED(ENUM_VAL, NAME_STR, MEMBER_PTR_A, MEMBER_PTR_B, BITS_A, BITS_B) \
  template <>                                                                                     \
  struct FieldTrait<ENUM_VAL> {                                                                   \
    using strategy_type = detail::PackedFieldsFunctor<BITS_A, BITS_B>;                            \
    static constexpr auto members = std::make_tuple(MEMBER_PTR_A, MEMBER_PTR_B);                  \
    using value_type = strategy_type::value_type;                                                 \
    static constexpr std::string_view name = NAME_STR;                                            \
  }

namespace o2::rans::benchmarks
{

/**
 * @brief Struct holding TPC data for rANS benchmarks
 */
struct TPCData {
  TPCData() = default;
  TPCData(size_t nTracks,
          size_t nAttachedClusters,
          size_t nUnattachedClusters,
          size_t nAttachedClustersReduced,
          size_t nTriggers);

  void resizeTracks(size_t n);
  void resizeAttachedClusters(size_t n);
  void resizeUnattachedClusters(size_t n);
  void resizeAttachedClustersReduced(size_t n);
  void resizeTriggers(size_t n);

  size_t nTracks = 0;
  size_t nAttachedClusters = 0;
  size_t nUnattachedClusters = 0;
  size_t nAttachedClustersReduced = 0;
  size_t nSliceRows = 36 * 152;
  size_t nTriggers = 0;

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
  std::vector<uint32_t> timeA{};
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

  std::vector<int32_t> deltaOrbit{};
  std::vector<int16_t> deltaBC{};
  std::vector<uint8_t> triggerType{};
};

/**
 * @brief load TPC data from a ROOT file
 *
 * @param fileName  name of the ROOT file
 * @return TPCData populated TPCData struct
 */
TPCData loadTPCData(const std::filesystem::path& fileName);

/**
 * @brief Enumeration of TPC fields
 *
 * Holds the names of all data fields produced by TPC online reconstruction.
 * For better compression, some fields can be packed together, these are listed
 * under "collapsed fields".
 */
enum class TPCField {
  QTotA,
  QMaxA,
  FlagsA,
  RowDiffA,
  SliceLegDiffA,
  PadResA,
  TimeResA,
  SigmaPadA,
  SigmaTimeA,
  QPtA,
  RowA,
  SliceA,
  TimeA,
  PadA,
  QTotU,
  QMaxU,
  FlagsU,
  PadDiffU,
  TimeDiffU,
  SigmaPadU,
  SigmaTimeU,
  NTrackClusters,
  NSliceRowClusters,
  DeltaOrbit,
  DeltaBC,
  TriggerType,
  // collapsed fields
  QTotQMaxA,
  RowDiffSliceLegDiffA,
  SigmaPadSigmaTimeA,
  QTotQMaxU,
  SigmaPadSigmaTimeU,
};

namespace detail
{
struct IdentityFunctor {
};

inline static constexpr int NBitsQTot = 16;
inline static constexpr int NBitsQMax = 10;
inline static constexpr int NBitsSigmaPad = 8;
inline static constexpr int NBitsSigmaTime = 8;
inline static constexpr int NBitsRowDiff = 8;
inline static constexpr int NBitsSliceLegDiff = 7;

template <size_t N>
using fit_uint_t = std::conditional_t<(N > 16), uint32_t,
                                      std::conditional_t<(N > 8), uint16_t, uint8_t>>;

template <size_t A, size_t B>
using combinedType_t = fit_uint_t<A + B>;

template <typename T>
struct extract_inner_type;

template <typename MemberContainer, typename ClassType>
struct extract_inner_type<MemberContainer ClassType::*> {
  using type = typename MemberContainer::value_type;
};

template <typename T>
using extract_inner_type_t = typename extract_inner_type<T>::type;

template <size_t highBits_V, size_t lowBits_V>
class PackedFieldsFunctor
{
 private:
  static constexpr size_t shift = lowBits_V;

 public:
  using value_type = combinedType_t<highBits_V, lowBits_V>;

  template <typename iterA_T, typename iterB_T>
  inline value_type operator()(iterA_T iterA, iterB_T iterB) const
  {
    return *iterB + (static_cast<value_type>(*iterA) << shift);
  };

  template <typename iterA_T, typename iterB_T>
  inline void operator()(iterA_T iterA, iterB_T iterB, value_type value) const
  {
    *iterA = value >> shift;
    *iterB = value & ((0x1 << shift) - 0x1);
  };
};

template <typename valueA_T, typename valueB_T, class F>
auto makeInputIterators(const std::vector<valueA_T>& a, const std::vector<valueB_T>& b, F functor)
{
  using namespace o2::rans;

  auto advanceIter = [](auto iter, size_t nElements) {
    auto tmp = iter;
    std::advance(tmp, nElements);
    return tmp;
  };

  const auto iterA = a.data();
  const auto iterB = b.data();
  const size_t nElements = std::min(a.size(), b.size());

  assert(iterA != nullptr && iterB != nullptr);
  assert(a.size() == b.size() && "Input vectors must have the same size");

  auto begin = CombinedInputIterator{iterA, iterB, functor};
  auto end = [&]() {
    auto tmpA = advanceIter(iterA, nElements);
    auto tmpB = advanceIter(iterB, nElements);
    return CombinedInputIterator{tmpA, tmpB, functor};
  }();

  return std::make_tuple(begin, end);
};

template <typename T>
auto makeInputIterators(const std::vector<T>& v)
{
  return std::make_tuple(v.data(), v.data() + v.size());
}

template <class Functor_T, typename T1, typename T2>
auto makeOutputIterators(Functor_T, std::vector<T1>& a, std::vector<T2>& b)
{
  using namespace o2::rans;

  auto advanceIter = [](auto iter, size_t nElements) {
    auto tmp = iter;
    std::advance(tmp, nElements);
    return tmp;
  };

  auto iterA = a.data();
  auto iterB = b.data();
  const size_t nElements = std::min(a.size(), b.size());

  assert(iterA != nullptr && iterB != nullptr);
  assert(a.size() == b.size() && "Input vectors must have the same size");

  auto begin = CombinedOutputIteratorFactory<typename Functor_T::value_type>::makeIter(iterA, iterB, Functor_T{});
  auto end = [&]() {
    auto tmpA = advanceIter(iterA, nElements);
    auto tmpB = advanceIter(iterB, nElements);
    return CombinedOutputIteratorFactory<typename Functor_T::value_type>::makeIter(tmpA, tmpB, Functor_T{});
  }();

  return std::make_tuple(begin, end);
};

template <typename T>
auto makeOutputIterators(std::vector<T>& v)
{
  return std::make_tuple(v.data(), v.data() + v.size());
}

template <typename T>
auto makeIteratorsDispatch(IdentityFunctor, const std::vector<T>& vec)
{
  return makeInputIterators(vec);
}

template <typename Functor_T, typename T1, typename T2>
auto makeIteratorsDispatch(Functor_T, const std::vector<T1>& v1, const std::vector<T2>& v2)
{
  return makeInputIterators(v1, v2, Functor_T{});
}

template <typename T>
auto makeIteratorsDispatch(IdentityFunctor, std::vector<T>& vec)
{
  return makeOutputIterators(vec);
}

template <typename Functor_T, typename T1, typename T2>
auto makeIteratorsDispatch(Functor_T, std::vector<T1>& v1, std::vector<T2>& v2)
{
  return makeOutputIterators(Functor_T{}, v1, v2);
}

} // namespace detail

namespace traits
{
template <TPCField field_V>
struct FieldTrait;

DEFINE_FIELD_TRAIT_SIMPLE(TPCField::QTotA, "qTotA", &TPCData::qTotA);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::QMaxA, "qMaxA", &TPCData::qMaxA);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::FlagsA, "flagsA", &TPCData::flagsA);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::RowDiffA, "rowDiffA", &TPCData::rowDiffA);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::SliceLegDiffA, "sliceLegDiffA", &TPCData::sliceLegDiffA);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::PadResA, "padResA", &TPCData::padResA);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::TimeResA, "timeResA", &TPCData::timeResA);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::SigmaPadA, "sigmaPadA", &TPCData::sigmaPadA);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::SigmaTimeA, "sigmaTimeA", &TPCData::sigmaTimeA);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::QPtA, "qPtA", &TPCData::qPtA);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::RowA, "rowA", &TPCData::rowA);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::SliceA, "sliceA", &TPCData::sliceA);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::TimeA, "timeA", &TPCData::timeA);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::PadA, "padA", &TPCData::padA);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::QTotU, "qTotU", &TPCData::qTotU);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::QMaxU, "qMaxU", &TPCData::qMaxU);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::FlagsU, "flagsU", &TPCData::flagsU);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::PadDiffU, "padDiffU", &TPCData::padDiffU);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::TimeDiffU, "timeDiffU", &TPCData::timeDiffU);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::SigmaPadU, "sigmaPadU", &TPCData::sigmaPadU);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::SigmaTimeU, "sigmaTimeU", &TPCData::sigmaTimeU);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::NTrackClusters, "nTrackClusters", &TPCData::nTrackClusters);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::NSliceRowClusters, "nSliceRowClusters", &TPCData::nSliceRowClusters);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::DeltaOrbit, "TrigOrbitInc", &TPCData::deltaOrbit);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::DeltaBC, "TrigBCInc", &TPCData::deltaBC);
DEFINE_FIELD_TRAIT_SIMPLE(TPCField::TriggerType, "TrigType", &TPCData::triggerType);

// Packed Fields
DEFINE_FIELD_TRAIT_PACKED(TPCField::QTotQMaxA,
                          "QTotQMaxA",
                          &TPCData::qTotA,
                          &TPCData::qMaxA,
                          detail::NBitsQTot,
                          detail::NBitsQMax);
DEFINE_FIELD_TRAIT_PACKED(TPCField::RowDiffSliceLegDiffA,
                          "RowDiffSliceLegDiffA",
                          &TPCData::rowDiffA,
                          &TPCData::sliceLegDiffA,
                          detail::NBitsRowDiff,
                          detail::NBitsSliceLegDiff);
DEFINE_FIELD_TRAIT_PACKED(TPCField::SigmaPadSigmaTimeA,
                          "SigmaPadSigmaTimeA",
                          &TPCData::sigmaPadA,
                          &TPCData::sigmaTimeA,
                          detail::NBitsSigmaPad,
                          detail::NBitsSigmaTime);
DEFINE_FIELD_TRAIT_PACKED(TPCField::QTotQMaxU,
                          "QTotQMaxU",
                          &TPCData::qTotU,
                          &TPCData::qMaxU,
                          detail::NBitsQTot,
                          detail::NBitsQMax);
DEFINE_FIELD_TRAIT_PACKED(TPCField::SigmaPadSigmaTimeU,
                          "SigmaPadSigmaTimeU",
                          &TPCData::sigmaPadU,
                          &TPCData::sigmaTimeU,
                          detail::NBitsSigmaPad,
                          detail::NBitsSigmaTime);

// shorthands:

template <TPCField field_V>
using tpcField_value_t = typename FieldTrait<field_V>::value_type;

template <TPCField field_V>
inline static constexpr std::string_view tpcField_name = FieldTrait<field_V>::name;

using TPCFields = boost::mp11::mp_list_v<
  TPCField::QTotA,
  TPCField::QMaxA,
  TPCField::FlagsA,
  TPCField::RowDiffA,
  TPCField::SliceLegDiffA,
  TPCField::PadResA,
  TPCField::TimeResA,
  TPCField::SigmaPadA,
  TPCField::SigmaTimeA,
  TPCField::QPtA,
  TPCField::RowA,
  TPCField::SliceA,
  TPCField::TimeA,
  TPCField::PadA,
  TPCField::QTotU,
  TPCField::QMaxU,
  TPCField::FlagsU,
  TPCField::PadDiffU,
  TPCField::TimeDiffU,
  TPCField::SigmaPadU,
  TPCField::SigmaTimeU,
  TPCField::NTrackClusters,
  TPCField::NSliceRowClusters,
  TPCField::DeltaOrbit,
  TPCField::DeltaBC,
  TPCField::TriggerType>;

using TPCFieldsCollapsedOnly = boost::mp11::mp_list_v<
  TPCField::QTotQMaxA,
  TPCField::RowDiffSliceLegDiffA,
  TPCField::SigmaPadSigmaTimeA,
  TPCField::QTotQMaxU,
  TPCField::SigmaPadSigmaTimeU>;

using TPCFieldsAll = boost::mp11::mp_append<TPCFields, TPCFieldsCollapsedOnly>;

using TPCFieldsCollapsed = boost::mp11::mp_list_v<
  TPCField::QTotQMaxA,
  TPCField::FlagsA,
  TPCField::RowDiffSliceLegDiffA,
  TPCField::PadResA,
  TPCField::TimeResA,
  TPCField::SigmaPadSigmaTimeA,
  TPCField::QPtA,
  TPCField::RowA,
  TPCField::SliceA,
  TPCField::TimeA,
  TPCField::PadA,
  TPCField::QTotQMaxU,
  TPCField::FlagsU,
  TPCField::PadDiffU,
  TPCField::TimeDiffU,
  TPCField::SigmaPadSigmaTimeU,
  TPCField::NTrackClusters,
  TPCField::NSliceRowClusters,
  TPCField::DeltaOrbit,
  TPCField::DeltaBC,
  TPCField::TriggerType>;

} // namespace traits

template <TPCField field_V>
decltype(auto) getIter(const TPCData& data)
{
  using trait_t = traits::FieldTrait<field_V>;

  auto unpacker = [&](auto... memberPtrs) {
    return detail::makeIteratorsDispatch(
      typename trait_t::strategy_type{},
      (data.*memberPtrs)...);
  };

  return std::apply(unpacker, trait_t::members);
};

template <TPCField field_V>
decltype(auto) getIter(TPCData& data)
{
  using trait_t = traits::FieldTrait<field_V>;
  auto unpacker = [&](auto... memberPtrs) {
    return detail::makeIteratorsDispatch(
      typename trait_t::strategy_type{},
      (data.*memberPtrs)...);
  };
  return std::apply(unpacker, trait_t::members);
};

bool operator==(const TPCData& lhs, const TPCData& rhs);

} // namespace o2::rans::benchmarks

#endif /* RANS_BENCHMARKS_HELPERS_TPCDATA_H_ */