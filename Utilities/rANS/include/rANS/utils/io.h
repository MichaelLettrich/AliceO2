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

/// \file io.h
/// \brief
/// \author michael.lettrich@cern.ch

#include <fstream>

#include <rapidjson/rapidjson.h>
#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/prettywriter.h>

namespace o2
{
namespace rans
{

namespace utils
{

template <typename IT>
rapidjson::Value toJSON(IT begin, IT end, uint32_t min, uint32_t max, rapidjson::Document& jsonDocument)
{
  auto& allocator = jsonDocument.GetAllocator();
  rapidjson::Value dict{rapidjson::kObjectType};

  dict.AddMember("min", min, allocator);
  dict.AddMember("max", max, allocator);
  dict.AddMember(
    "values",
    [&]() {
      rapidjson::Value tmpArray{rapidjson::kArrayType};
      for (auto iter = begin; iter != end; ++iter) {
        tmpArray.PushBack(*iter, allocator);
      }
      return tmpArray;
    }(),
    allocator);

  return dict;
};

} // namespace utils
} // namespace rans
} // namespace o2