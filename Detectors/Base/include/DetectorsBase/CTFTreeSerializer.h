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

#ifndef _ALICEO2_CTFTreeSerializer_BASE_H_
#define _ALICEO2_CTFTreeSerializer_BASE_H_

#include <memory>
#include <TTree.h>

namespace o2::ctf
{

class CTFTreeSerializer
{
 public:
  CTFTreeSerializer(const std::string& detNam, const std::string& prefix = "");

  void initTree();
  TTree* getTree();
  void writeTree();

 private:
  std::unique_ptr<TTree> mTree{};
  std::string mDetName{};
  std::string mPrefix{};
  uint32_t mCounter{};
};

} // namespace o2::ctf

#endif /* _ALICEO2_CTFTreeSerializer_BASE_H_ */
