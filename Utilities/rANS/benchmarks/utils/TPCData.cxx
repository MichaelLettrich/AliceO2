#include "TPCData.h"

#include <cassert>
#include <tuple>
#include <filesystem>

#include <ROOT/RDataFrame.hxx>
#include <boost/mp11.hpp>

namespace bmp = boost::mp11;

namespace o2::rans::benchmarks
{

TPCData::TPCData(size_t nTracks,
                 size_t nAttachedClusters,
                 size_t nUnattachedClusters,
                 size_t nAttachedClustersReduced,
                 size_t nTriggers)
{
  resizeTracks(nTracks);
  resizeAttachedClusters(nAttachedClusters);
  resizeUnattachedClusters(nUnattachedClusters);
  resizeAttachedClustersReduced(nAttachedClustersReduced);
  resizeTriggers(nTriggers);
  nSliceRowClusters.resize(nSliceRows);
};

void TPCData::resizeTracks(size_t n)
{
  nTracks = n;

  qPtA.resize(this->nTracks);
  rowA.resize(this->nTracks);
  sliceA.resize(this->nTracks);
  timeA.resize(this->nTracks);
  padA.resize(this->nTracks);
  nTrackClusters.resize(this->nTracks);
}

void TPCData::resizeAttachedClusters(size_t n)
{
  nAttachedClusters = n;

  qTotA.resize(this->nAttachedClusters);
  qMaxA.resize(this->nAttachedClusters);
  flagsA.resize(this->nAttachedClusters);
  sigmaPadA.resize(this->nAttachedClusters);
  sigmaTimeA.resize(this->nAttachedClusters);
}

void TPCData::resizeUnattachedClusters(size_t n)
{
  nUnattachedClusters = n;

  qTotU.resize(this->nUnattachedClusters);
  qMaxU.resize(this->nUnattachedClusters);
  flagsU.resize(this->nUnattachedClusters);
  padDiffU.resize(this->nUnattachedClusters);
  timeDiffU.resize(this->nUnattachedClusters);
  sigmaPadU.resize(this->nUnattachedClusters);
  sigmaTimeU.resize(this->nUnattachedClusters);
}

void TPCData::resizeAttachedClustersReduced(size_t n)
{
  nAttachedClustersReduced = n;

  rowDiffA.resize(this->nAttachedClustersReduced);
  sliceLegDiffA.resize(this->nAttachedClustersReduced);
  padResA.resize(this->nAttachedClustersReduced);
  timeResA.resize(this->nAttachedClustersReduced);
}

void TPCData::resizeTriggers(size_t n)
{
  nTriggers = n;

  deltaOrbit.resize(this->nTriggers);
  deltaBC.resize(this->nTriggers);
  triggerType.resize(this->nTriggers);
}

TPCData loadTPCData(const std::filesystem::path& fileName)
{
  TPCData data;
  ROOT::RDataFrame df("TPC", fileName.c_str());

  // makes a tuple of RResultPtrs, with one entry per field
  auto result_ptrs = []<auto... Field>(bmp::mp_list_v<Field...>, ROOT::RDataFrame& df) {
    return std::make_tuple(
      df.Take<std::vector<typename traits::FieldTrait<Field>::value_type>>(
        std::string(traits::FieldTrait<Field>::name))...);
  }(traits::TPCFields{}, df);

  [&]<auto... Field>(auto& results, TPCData& data, bmp::mp_list_v<Field...>) {
    [&]<std::size_t... Is>(std::index_sequence<Is...>) {
      // For each field, move the vector from the result pointer (at index 0) to the member.
      ((data.*std::get<0>(traits::FieldTrait<Field>::members) = std::move(std::get<Is>(results)->at(0))), ...);
    }(std::make_index_sequence<sizeof...(Field)>{});
  }(result_ptrs, data, traits::TPCFields{});

  // Set derived sizes
  data.nTracks = data.qPtA.size();
  data.nAttachedClusters = data.qTotA.size();
  data.nUnattachedClusters = data.qTotU.size();
  data.nAttachedClustersReduced = data.rowDiffA.size();
  data.nTriggers = data.deltaOrbit.size();

  return data;
}

bool operator==(const TPCData& lhs, const TPCData& rhs)
{
  auto lhs_sizes = std::tie(lhs.nTracks,
                            lhs.nAttachedClusters,
                            lhs.nUnattachedClusters,
                            lhs.nAttachedClustersReduced,
                            lhs.nTriggers);
  auto rhs_sizes = std::tie(rhs.nTracks,
                            rhs.nAttachedClusters,
                            rhs.nUnattachedClusters,
                            rhs.nAttachedClustersReduced,
                            rhs.nTriggers);
  if (!(lhs_sizes == rhs_sizes)) {
    return false;
  }

  return []<auto... Fields>(const TPCData& lhs, const TPCData& rhs, boost::mp11::mp_list_v<Fields...>) {
    return ((
              lhs.*std::get<0>(traits::FieldTrait<Fields>::members) ==
              rhs.*std::get<0>(traits::FieldTrait<Fields>::members)) &&
            ...);
  }(lhs, rhs, traits::TPCFields{});
}

} // namespace o2::rans::benchmarks