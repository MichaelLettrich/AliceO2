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

/// \file   CTFCoder.h
/// \author ruben.shahoyan@cern.ch
/// \brief class for entropy encoding/decoding of TPC compressed clusters data

#ifndef O2_TPC_CTFCODER_H
#define O2_TPC_CTFCODER_H

#include <algorithm>
#include <iterator>
#include <string>
#include <cassert>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <vector>
#include "DataFormatsTPC/CTF.h"
#include "DataFormatsTPC/CompressedClusters.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsBase/CTFCoderBase.h"
#include "rANS/iterator.h"

#include <ittnotify.h>

inline __itt_domain* tpcDomain = __itt_domain_create("o2.tpc.ctfcoder");

inline __itt_string_handle* encodeTaskQTotA = __itt_string_handle_create("encodeQTotA");
inline __itt_string_handle* encodeTaskQMaxA = __itt_string_handle_create("encodeQMaxA");
inline __itt_string_handle* encodeTaskFlagsA = __itt_string_handle_create("encodeFlagsA");
inline __itt_string_handle* encodeTaskRowDiffA = __itt_string_handle_create("encodeRowDiffA");
inline __itt_string_handle* encodeTaskSliceLegDiffA = __itt_string_handle_create("encodeSliceLegDiffA");
inline __itt_string_handle* encodeTaskPadResA = __itt_string_handle_create("encodePadResA");
inline __itt_string_handle* encodeTaskTimeResA = __itt_string_handle_create("encodeTimeResA");
inline __itt_string_handle* encodeTaskSigmaPadA = __itt_string_handle_create("encodeSigmaPadA");
inline __itt_string_handle* encodeTaskSigmaTimeA = __itt_string_handle_create("encodeSigmaTimeA");
inline __itt_string_handle* encodeTaskQPtA = __itt_string_handle_create("encodeQPtA");
inline __itt_string_handle* encodeTaskRowA = __itt_string_handle_create("encodeRowA");
inline __itt_string_handle* encodeTaskSliceA = __itt_string_handle_create("encodeSliceA");
inline __itt_string_handle* encodeTaskTimeA = __itt_string_handle_create("encodeTimeA");
inline __itt_string_handle* encodeTaskPadA = __itt_string_handle_create("encodePadA");
inline __itt_string_handle* encodeTaskQTotU = __itt_string_handle_create("encodeQTotU");
inline __itt_string_handle* encodeTaskQMaxU = __itt_string_handle_create("encodeQMaxU");
inline __itt_string_handle* encodeTaskFlagsU = __itt_string_handle_create("encodeFlagsU");
inline __itt_string_handle* encodeTaskPadDiffU = __itt_string_handle_create("encodePadDiffU");
inline __itt_string_handle* encodeTaskTimeDiffU = __itt_string_handle_create("encodeTimeDiffU");
inline __itt_string_handle* encodeTaskSigmaPadU = __itt_string_handle_create("encodeSigmaPadU");
inline __itt_string_handle* encodeTaskSigmaTimeU = __itt_string_handle_create("encodeSigmaTimeU");
inline __itt_string_handle* encodeTaskNTrackClusters = __itt_string_handle_create("encodeNTrackClusters");
inline __itt_string_handle* encodeTaskNSliceRowClusters = __itt_string_handle_create("encodeNSliceRowClusters");

inline __itt_string_handle* decodeTaskQTotA = __itt_string_handle_create("decodeQTotA");
inline __itt_string_handle* decodeTaskQMaxA = __itt_string_handle_create("decodeQMaxA");
inline __itt_string_handle* decodeTaskFlagsA = __itt_string_handle_create("decodeFlagsA");
inline __itt_string_handle* decodeTaskRowDiffA = __itt_string_handle_create("decodeRowDiffA");
inline __itt_string_handle* decodeTaskSliceLegDiffA = __itt_string_handle_create("decodeSliceLegDiffA");
inline __itt_string_handle* decodeTaskPadResA = __itt_string_handle_create("decodePadResA");
inline __itt_string_handle* decodeTaskTimeResA = __itt_string_handle_create("decodeTimeResA");
inline __itt_string_handle* decodeTaskSigmaPadA = __itt_string_handle_create("decodeSigmaPadA");
inline __itt_string_handle* decodeTaskSigmaTimeA = __itt_string_handle_create("decodeSigmaTimeA");
inline __itt_string_handle* decodeTaskQPtA = __itt_string_handle_create("decodeQPtA");
inline __itt_string_handle* decodeTaskRowA = __itt_string_handle_create("decodeRowA");
inline __itt_string_handle* decodeTaskSliceA = __itt_string_handle_create("decodeSliceA");
inline __itt_string_handle* decodeTaskTimeA = __itt_string_handle_create("decodeTimeA");
inline __itt_string_handle* decodeTaskPadA = __itt_string_handle_create("decodePadA");
inline __itt_string_handle* decodeTaskQTotU = __itt_string_handle_create("decodeQTotU");
inline __itt_string_handle* decodeTaskQMaxU = __itt_string_handle_create("decodeQMaxU");
inline __itt_string_handle* decodeTaskFlagsU = __itt_string_handle_create("decodeFlagsU");
inline __itt_string_handle* decodeTaskPadDiffU = __itt_string_handle_create("decodePadDiffU");
inline __itt_string_handle* decodeTaskTimeDiffU = __itt_string_handle_create("decodeTimeDiffU");
inline __itt_string_handle* decodeTaskSigmaPadU = __itt_string_handle_create("decodeSigmaPadU");
inline __itt_string_handle* decodeTaskSigmaTimeU = __itt_string_handle_create("decodeSigmaTimeU");
inline __itt_string_handle* decodeTaskNTrackClusters = __itt_string_handle_create("decodeNTrackClusters");
inline __itt_string_handle* decodeTaskNSliceRowClusters = __itt_string_handle_create("decodeNSliceRowClusters");

class TTree;

namespace o2
{
namespace tpc
{

namespace detail
{

template <int A, int B>
struct combinedType {
  using type = std::conditional_t<(A + B > 16), uint32_t, std::conditional_t<(A + B > 8), uint16_t, uint8_t>>;
};

template <int A, int B>
using combinedType_t = typename combinedType<A, B>::type;

template <typename value_T, size_t shift>
class ShiftFunctor
{
 public:
  template <typename iterA_T, typename iterB_T>
  inline value_T operator()(iterA_T iterA, iterB_T iterB) const
  {
    return *iterB + (static_cast<value_T>(*iterA) << shift);
  };

  template <typename iterA_T, typename iterB_T>
  inline void operator()(iterA_T iterA, iterB_T iterB, value_T value) const
  {
    *iterA = value >> shift;
    *iterB = value & ((0x1 << shift) - 0x1);
  };
};

template <typename iterA_T, typename iterB_T, typename F>
auto makeInputIterators(iterA_T iterA, iterB_T iterB, size_t nElements, F functor)
{
  using namespace o2::rans;

  auto advanceIter = [](auto iter, size_t nElements) {
    auto tmp = iter;
    std::advance(tmp, nElements);
    return tmp;
  };

  return std::make_tuple(CombinedInputIterator{iterA, iterB, functor},
                         CombinedInputIterator{advanceIter(iterA, nElements), advanceIter(iterB, nElements), functor});
};

template <int bits_A, int bits_B>
struct MergedColumnsDecoder {

  using combined_t = combinedType_t<bits_A, bits_B>;

  template <typename iterA_T, typename iterB_T, typename F>
  static void decode(iterA_T iterA, iterB_T iterB, CTF::Slots slot, F decodingFunctor)
  {
    using namespace o2::rans;
    ShiftFunctor<combined_t, bits_B> f{};
    auto iter = CombinedOutputIteratorFactory<combined_t>::makeIter(iterA, iterB, f);

    decodingFunctor(iter, slot);
  }
};

} // namespace detail

class CTFCoder : public o2::ctf::CTFCoderBase
{
 public:
  CTFCoder(o2::ctf::CTFCoderBase::OpType op) : o2::ctf::CTFCoderBase(op, CTF::getNBlocks(), o2::detectors::DetID::TPC) {}
  ~CTFCoder() final = default;

  /// entropy-encode compressed clusters to flat buffer
  template <typename VEC>
  o2::ctf::CTFIOSize encode(VEC& buff, const CompressedClusters& ccl, const CompressedClusters& cclFiltered, std::vector<bool>* rejectHits = nullptr, std::vector<bool>* rejectTracks = nullptr, std::vector<bool>* rejectTrackHits = nullptr, std::vector<bool>* rejectTrackHitsReduced = nullptr);

  template <typename VEC>
  o2::ctf::CTFIOSize decode(const CTF::base& ec, VEC& buff);

  void createCoders(const std::vector<char>& bufVec, o2::ctf::CTFCoderBase::OpType op) final;

  size_t estimateCompressedSize(const CompressedClusters& ccl);

  static size_t constexpr Alignment = 16;
  static size_t estimateSize(CompressedClusters& c);
  static void setCompClusAddresses(CompressedClusters& c, void*& buff);

  template <size_t ALG = Alignment, typename T>
  static size_t alignSize(T*& var, size_t n = 1)
  {
    auto sz = sizeof(T) * n;
    auto res = sz % ALG;
    return res ? sz + (ALG - res) : sz;
  }

  template <size_t ALG = Alignment, typename T>
  static void setAlignedPtr(void*& ptr, T*& var, size_t n = 1)
  {
    auto sz = sizeof(T) * n;
    auto res = sz % ALG;
    var = reinterpret_cast<T*>(ptr);
    auto& ptrR = reinterpret_cast<size_t&>(ptr);
    ptrR += res ? sz + (ALG - res) : sz;
  }

  bool getCombineColumns() const { return mCombineColumns; }
  void setCombineColumns(bool v) { mCombineColumns = v; }

 private:
  void checkDataDictionaryConsistency(const CTFHeader& h);

  template <int NU, int NL, typename CU, typename CL>
  static void splitColumns(const std::vector<detail::combinedType_t<NU, NL>>& vm, CU*& vu, CL*& vl);

  template <typename source_T>
  void buildCoder(ctf::CTFCoderBase::OpType coderType, const CTF::container_t& ctf, CTF::Slots slot);

  bool mCombineColumns = false; // combine correlated columns
};

template <typename source_T>
void CTFCoder::buildCoder(ctf::CTFCoderBase::OpType coderType, const CTF::container_t& ctf, CTF::Slots slot)
{
  this->createCoder(coderType, ctf.getFrequencyTable<source_T>(slot, mANSVersion), static_cast<int>(slot));
}

/// entropy-encode clusters to buffer with CTF
template <typename VEC>
o2::ctf::CTFIOSize CTFCoder::encode(VEC& buff, const CompressedClusters& ccl, const CompressedClusters& cclFiltered, std::vector<bool>* rejectHits, std::vector<bool>* rejectTracks, std::vector<bool>* rejectTrackHits, std::vector<bool>* rejectTrackHitsReduced)
{
  __itt_resume();

  using MD = o2::ctf::Metadata::OptStore;
  using namespace detail;
  // what to do which each field: see o2::ctf::Metadata explanation
  constexpr MD optField[CTF::getNBlocks()] = {
    MD::EENCODE_OR_PACK, // qTotA
    MD::EENCODE_OR_PACK, // qMaxA
    MD::EENCODE_OR_PACK, // flagsA
    MD::EENCODE_OR_PACK, // rowDiffA
    MD::EENCODE_OR_PACK, // sliceLegDiffA
    MD::EENCODE_OR_PACK, // padResA
    MD::EENCODE_OR_PACK, // timeResA
    MD::EENCODE_OR_PACK, // sigmaPadA
    MD::EENCODE_OR_PACK, // sigmaTimeA
    MD::EENCODE_OR_PACK, // qPtA
    MD::EENCODE_OR_PACK, // rowA
    MD::EENCODE_OR_PACK, // sliceA
    MD::EENCODE_OR_PACK, // timeA
    MD::EENCODE_OR_PACK, // padA
    MD::EENCODE_OR_PACK, // qTotU
    MD::EENCODE_OR_PACK, // qMaxU
    MD::EENCODE_OR_PACK, // flagsU
    MD::EENCODE_OR_PACK, // padDiffU
    MD::EENCODE_OR_PACK, // timeDiffU
    MD::EENCODE_OR_PACK, // sigmaPadU
    MD::EENCODE_OR_PACK, // sigmaTimeU
    MD::EENCODE_OR_PACK, // nTrackClusters
    MD::EENCODE_OR_PACK  // nSliceRowClusters
  };

  // book output size with some margin
  auto szIni = estimateCompressedSize(cclFiltered);
  buff.resize(szIni);

  auto ec = CTF::create(buff);
  uint32_t flags = 0;
  if (mCombineColumns) {
    flags |= CTFHeader::CombinedColumns;
  }
  ec->setHeader(CTFHeader{o2::detectors::DetID::TPC, 0, 1, 0, // dummy timestamp, version 1.0
                          cclFiltered, flags});
  assignDictVersion(static_cast<o2::ctf::CTFDictHeader&>(ec->getHeader()));
  ec->setANSHeader(mANSVersion);

  o2::ctf::CTFIOSize iosize;
  auto encodeTPC = [&buff, &optField, &coders = mCoders, mfc = this->getMemMarginFactor(), &iosize](auto begin, auto end, CTF::Slots slot, size_t probabilityBits, std::vector<bool>* reject = nullptr) {
    // at every encoding the buffer might be autoexpanded, so we don't work with fixed pointer ec
    const auto slotVal = static_cast<int>(slot);
    if (reject && begin != end) {
      std::vector<std::decay_t<decltype(*begin)>> tmp;
      tmp.reserve(std::distance(begin, end));
      for (auto i = begin; i != end; i++) {
        if (!(*reject)[std::distance(begin, i)]) {
          tmp.emplace_back(*i);
        }
      }
      iosize += CTF::get(buff.data())->encode(tmp.begin(), tmp.end(), slotVal, probabilityBits, optField[slotVal], &buff, coders[slotVal], mfc);
    } else {
      iosize += CTF::get(buff.data())->encode(begin, end, slotVal, probabilityBits, optField[slotVal], &buff, coders[slotVal], mfc);
    }
  };

  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskQTotA);
  if (mCombineColumns) {
    const auto [begin, end] = makeInputIterators(ccl.qTotA, ccl.qMaxA, ccl.nAttachedClusters,
                                                 ShiftFunctor<combinedType_t<CTF::NBitsQTot, CTF::NBitsQMax>, CTF::NBitsQMax>{});
    encodeTPC(begin, end, CTF::BLCqTotA, 0, rejectTrackHits);
  } else {
    encodeTPC(ccl.qTotA, ccl.qTotA + ccl.nAttachedClusters, CTF::BLCqTotA, 0, rejectTrackHits);
  }
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskQMaxA);
  encodeTPC(ccl.qMaxA, ccl.qMaxA + (mCombineColumns ? 0 : ccl.nAttachedClusters), CTF::BLCqMaxA, 0, rejectTrackHits);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskFlagsA);
  encodeTPC(ccl.flagsA, ccl.flagsA + ccl.nAttachedClusters, CTF::BLCflagsA, 0, rejectTrackHits);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskRowDiffA);
  if (mCombineColumns) {
    const auto [begin, end] = makeInputIterators(ccl.rowDiffA, ccl.sliceLegDiffA, ccl.nAttachedClustersReduced,
                                                 ShiftFunctor<combinedType_t<CTF::NBitsRowDiff, CTF::NBitsSliceLegDiff>, CTF::NBitsSliceLegDiff>{});
    encodeTPC(begin, end, CTF::BLCrowDiffA, 0, rejectTrackHitsReduced);
  } else {
    encodeTPC(ccl.rowDiffA, ccl.rowDiffA + ccl.nAttachedClustersReduced, CTF::BLCrowDiffA, 0, rejectTrackHitsReduced);
  }
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskSliceLegDiffA);
  encodeTPC(ccl.sliceLegDiffA, ccl.sliceLegDiffA + (mCombineColumns ? 0 : ccl.nAttachedClustersReduced), CTF::BLCsliceLegDiffA, 0, rejectTrackHitsReduced);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskPadResA);
  encodeTPC(ccl.padResA, ccl.padResA + ccl.nAttachedClustersReduced, CTF::BLCpadResA, 0, rejectTrackHitsReduced);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskTimeResA);
  encodeTPC(ccl.timeResA, ccl.timeResA + ccl.nAttachedClustersReduced, CTF::BLCtimeResA, 0, rejectTrackHitsReduced);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskSigmaPadA);
  if (mCombineColumns) {
    const auto [begin, end] = makeInputIterators(ccl.sigmaPadA, ccl.sigmaTimeA, ccl.nAttachedClusters,
                                                 ShiftFunctor<combinedType_t<CTF::NBitsSigmaPad, CTF::NBitsSigmaTime>, CTF::NBitsSigmaTime>{});
    encodeTPC(begin, end, CTF::BLCsigmaPadA, 0, rejectTrackHits);
  } else {
    encodeTPC(ccl.sigmaPadA, ccl.sigmaPadA + ccl.nAttachedClusters, CTF::BLCsigmaPadA, 0, rejectTrackHits);
  }
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskSigmaTimeA);
  encodeTPC(ccl.sigmaTimeA, ccl.sigmaTimeA + (mCombineColumns ? 0 : ccl.nAttachedClusters), CTF::BLCsigmaTimeA, 0, rejectTrackHits);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskQPtA);
  encodeTPC(ccl.qPtA, ccl.qPtA + ccl.nTracks, CTF::BLCqPtA, 0, rejectTracks);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskRowA);
  encodeTPC(ccl.rowA, ccl.rowA + ccl.nTracks, CTF::BLCrowA, 0, rejectTracks);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskSliceA);
  encodeTPC(ccl.sliceA, ccl.sliceA + ccl.nTracks, CTF::BLCsliceA, 0, rejectTracks);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskTimeA);
  encodeTPC(ccl.timeA, ccl.timeA + ccl.nTracks, CTF::BLCtimeA, 0, rejectTracks);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskPadA);
  encodeTPC(ccl.padA, ccl.padA + ccl.nTracks, CTF::BLCpadA, 0, rejectTracks);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskQTotU);
  if (mCombineColumns) {
    const auto [begin, end] = makeInputIterators(ccl.qTotU, ccl.qMaxU, ccl.nUnattachedClusters,
                                                 ShiftFunctor<combinedType_t<CTF::NBitsQTot, CTF::NBitsQMax>, CTF::NBitsQMax>{});
    encodeTPC(begin, end, CTF::BLCqTotU, 0, rejectHits);
  } else {
    encodeTPC(ccl.qTotU, ccl.qTotU + ccl.nUnattachedClusters, CTF::BLCqTotU, 0, rejectHits);
  }
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskQMaxU);
  encodeTPC(ccl.qMaxU, ccl.qMaxU + (mCombineColumns ? 0 : ccl.nUnattachedClusters), CTF::BLCqMaxU, 0, rejectHits);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskFlagsU);
  encodeTPC(ccl.flagsU, ccl.flagsU + ccl.nUnattachedClusters, CTF::BLCflagsU, 0, rejectHits);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskPadDiffU);
  encodeTPC(cclFiltered.padDiffU, cclFiltered.padDiffU + cclFiltered.nUnattachedClusters, CTF::BLCpadDiffU, 0);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskTimeDiffU);
  encodeTPC(cclFiltered.timeDiffU, cclFiltered.timeDiffU + cclFiltered.nUnattachedClusters, CTF::BLCtimeDiffU, 0);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskSigmaPadU);
  if (mCombineColumns) {
    const auto [begin, end] = makeInputIterators(ccl.sigmaPadU, ccl.sigmaTimeU, ccl.nUnattachedClusters,
                                                 ShiftFunctor<combinedType_t<CTF::NBitsSigmaPad, CTF::NBitsSigmaTime>, CTF::NBitsSigmaTime>{});
    encodeTPC(begin, end, CTF::BLCsigmaPadU, 0, rejectHits);
  } else {
    encodeTPC(ccl.sigmaPadU, ccl.sigmaPadU + ccl.nUnattachedClusters, CTF::BLCsigmaPadU, 0, rejectHits);
  }
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskSigmaTimeU);
  encodeTPC(ccl.sigmaTimeU, ccl.sigmaTimeU + (mCombineColumns ? 0 : ccl.nUnattachedClusters), CTF::BLCsigmaTimeU, 0, rejectHits);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskNTrackClusters);
  encodeTPC(ccl.nTrackClusters, ccl.nTrackClusters + ccl.nTracks, CTF::BLCnTrackClusters, 0, rejectTracks);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, encodeTaskNSliceRowClusters);
  encodeTPC(ccl.nSliceRowClusters, ccl.nSliceRowClusters + ccl.nSliceRows, CTF::BLCnSliceRowClusters, 0);
  __itt_task_end(tpcDomain);
  CTF::get(buff.data())->print(getPrefix(), mVerbosity);
  finaliseCTFOutput<CTF>(buff);
  iosize.rawIn = iosize.ctfIn;

  __itt_pause();

  return iosize;
}

/// decode entropy-encoded bloks to TPC CompressedClusters into the externally provided vector (e.g. PMR vector from DPL)
template <typename VEC>
o2::ctf::CTFIOSize CTFCoder::decode(const CTF::base& ec, VEC& buffVec)
{
  __itt_resume();

  using namespace detail;
  CompressedClusters cc;
  CompressedClustersCounters& ccCount = cc;
  auto& header = ec.getHeader();
  checkDictVersion(static_cast<const o2::ctf::CTFDictHeader&>(header));
  checkDataDictionaryConsistency(header);
  ccCount = static_cast<const CompressedClustersCounters&>(header);
  CompressedClustersFlat* ccFlat = nullptr;
  size_t sizeCFlatBody = alignSize(ccFlat);
  size_t sz = sizeCFlatBody + estimateSize(cc);                                             // total size of the buffVec accounting for the alignment
  size_t vsz = sizeof(typename std::remove_reference<decltype(buffVec)>::type::value_type); // size of the element of the buffer
  buffVec.resize(sz / vsz);
  ccFlat = reinterpret_cast<CompressedClustersFlat*>(buffVec.data());                           // RS? do we need to align this pointer, or PMR vector will be already aligned?
  auto buff = reinterpret_cast<void*>(reinterpret_cast<char*>(buffVec.data()) + sizeCFlatBody); // will be the start of the CompressedClustersFlat payload

  setCompClusAddresses(cc, buff);
  ccFlat->set(sz, cc); // set offsets
  ec.print(getPrefix(), mVerbosity);

  // decode encoded data directly to destination buff
  o2::ctf::CTFIOSize iosize;
  auto decodeTPC = [&ec, &coders = mCoders, &iosize](auto begin, CTF::Slots slot) {
    const auto slotVal = static_cast<int>(slot);
    iosize += ec.decode(begin, slotVal, coders[slotVal]);
  };

  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskQTotA);
  if (mCombineColumns) {
    detail::MergedColumnsDecoder<CTF::NBitsQTot, CTF::NBitsQMax>::decode(cc.qTotA, cc.qMaxA, CTF::BLCqTotA, decodeTPC);
  } else {
    decodeTPC(cc.qTotA, CTF::BLCqTotA);
    __itt_task_end(tpcDomain);
    __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskQMaxA);
    decodeTPC(cc.qMaxA, CTF::BLCqMaxA);
    __itt_task_end(tpcDomain);
  }

  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskFlagsA);
  decodeTPC(cc.flagsA, CTF::BLCflagsA);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskRowDiffA);
  if (mCombineColumns) {
    detail::MergedColumnsDecoder<CTF::NBitsRowDiff, CTF::NBitsSliceLegDiff>::decode(cc.rowDiffA, cc.sliceLegDiffA, CTF::BLCrowDiffA, decodeTPC);
  } else {
    decodeTPC(cc.rowDiffA, CTF::BLCrowDiffA);
    __itt_task_end(tpcDomain);
    __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskSliceLegDiffA);
    decodeTPC(cc.sliceLegDiffA, CTF::BLCsliceLegDiffA);
    __itt_task_end(tpcDomain);
  }

  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskPadResA);
  decodeTPC(cc.padResA, CTF::BLCpadResA);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskTimeResA);
  decodeTPC(cc.timeResA, CTF::BLCtimeResA);
  __itt_task_end(tpcDomain);

  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskSigmaPadA);
  if (mCombineColumns) {
    detail::MergedColumnsDecoder<CTF::NBitsSigmaPad, CTF::NBitsSigmaTime>::decode(cc.sigmaPadA, cc.sigmaTimeA, CTF::BLCsigmaPadA, decodeTPC);
  } else {
    decodeTPC(cc.sigmaPadA, CTF::BLCsigmaPadA);
    __itt_task_end(tpcDomain);
    __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskSigmaTimeA);
    decodeTPC(cc.sigmaTimeA, CTF::BLCsigmaTimeA);
    __itt_task_end(tpcDomain);
  }

  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskQPtA);
  decodeTPC(cc.qPtA, CTF::BLCqPtA);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskRowA);
  decodeTPC(cc.rowA, CTF::BLCrowA);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskSliceA);
  decodeTPC(cc.sliceA, CTF::BLCsliceA);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskTimeA);
  decodeTPC(cc.timeA, CTF::BLCtimeA);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskPadA);
  decodeTPC(cc.padA, CTF::BLCpadA);
  __itt_task_end(tpcDomain);

  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskQTotU);
  if (mCombineColumns) {
    detail::MergedColumnsDecoder<CTF::NBitsQTot, CTF::NBitsQMax>::decode(cc.qTotU, cc.qMaxU, CTF::BLCqTotU, decodeTPC);
  } else {
    decodeTPC(cc.qTotU, CTF::BLCqTotU);
    __itt_task_end(tpcDomain);
    __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskQMaxU);
    decodeTPC(cc.qMaxU, CTF::BLCqMaxU);
    __itt_task_end(tpcDomain);
  }

  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskFlagsU);
  decodeTPC(cc.flagsU, CTF::BLCflagsU);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskPadDiffU);
  decodeTPC(cc.padDiffU, CTF::BLCpadDiffU);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskTimeDiffU);
  decodeTPC(cc.timeDiffU, CTF::BLCtimeDiffU);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskSigmaPadU);

  if (mCombineColumns) {
    detail::MergedColumnsDecoder<CTF::NBitsSigmaPad, CTF::NBitsSigmaTime>::decode(cc.sigmaPadU, cc.sigmaTimeU, CTF::BLCsigmaPadU, decodeTPC);
  } else {
    decodeTPC(cc.sigmaPadU, CTF::BLCsigmaPadU);
    __itt_task_end(tpcDomain);
    __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskSigmaTimeU);
    decodeTPC(cc.sigmaTimeU, CTF::BLCsigmaTimeU);
    __itt_task_end(tpcDomain);
  }

  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskNTrackClusters);
  decodeTPC(cc.nTrackClusters, CTF::BLCnTrackClusters);
  __itt_task_end(tpcDomain);
  __itt_task_begin(tpcDomain, __itt_null, __itt_null, decodeTaskNSliceRowClusters);
  decodeTPC(cc.nSliceRowClusters, CTF::BLCnSliceRowClusters);
  __itt_task_end(tpcDomain);

  __itt_pause();
  iosize.rawIn = iosize.ctfIn;
  return iosize;
}

} // namespace tpc
} // namespace o2

#endif // O2_TPC_CTFCODER_H
