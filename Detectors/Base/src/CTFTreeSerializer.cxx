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

#include "DetectorsBase/CTFTreeSerializer.h"

#include <TFile.h>
#include <fmt/format.h>

namespace o2::ctf
{

CTFTreeSerializer::CTFTreeSerializer(const std::string& detName, const std::string& prefix) : mDetName{detName},
                                                                                              mPrefix{prefix}
{
}

void CTFTreeSerializer::initTree()
{
  mTree = std::make_unique<TTree>(mDetName.c_str(), fmt::format("CTF for {}", mDetName).c_str());
}

void CTFTreeSerializer::writeTree()
{
  const auto& ti = mTimingInfo;

  auto filename = fmt::format("{}ctf_{:08d}_{:010d}_{:010d}_{}.root", mPrefix, ti.runNumber, ti.firstTForbit, ti.tfCounter, mDetName);
  TFile f(filename.c_str(), "RECREATE");
  mTree->Fill();
  mTree->Write();
  f.Close();

  mTree.reset();
}

TTree* CTFTreeSerializer::getTree()
{
  return mTree.get();
}

void CTFTreeSerializer::updateTimingInfo(const o2::framework::TimingInfo& timingInfo)
{
  mTimingInfo = timingInfo;
}

} // namespace o2::ctf