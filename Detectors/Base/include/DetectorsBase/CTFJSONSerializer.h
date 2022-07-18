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

/// \file   CTFJSONSerializer.h
/// \author ruben.shahoyan@cern.ch
/// \brief class for entropy encoding/decoding of TPC compressed clusters data

#ifndef O2_CTF_JSONSERIALIZER_H
#define O2_CTF_JSONSERIALIZER_H

#include <string>
#include <fstream>

#include <rapidjson/rapidjson.h>
#include <rapidjson/writer.h>
#include <rapidjson/ostreamwrapper.h>

#include <fairlogger/Logger.h>

namespace o2
{
namespace ctf
{

class CTFJSONSerializer
{
 public:
  inline CTFJSONSerializer(std::string filename) : mFileStream{fmt::format("{}-{}.json", filename, ++id)},
                                                   mStreamWrapper{mFileStream},
                                                   mWriter{mStreamWrapper}
  {
    mWriter.StartObject();
  };

  inline ~CTFJSONSerializer()
  {
    mWriter.EndObject();
    mWriter.Flush();
    mFileStream.close();
  };

  inline void startDetector(std::string name)
  {
    mWriter.Key(name.c_str());
    mWriter.StartObject();
  }

  inline void endDetector()
  {
    mWriter.EndObject();
  };

  template <typename IT>
  void writeDataset(std::string name, IT begin, IT end)
  {
    mWriter.Key(name.c_str());
    mWriter.StartArray();
    for (IT iter = begin; iter != end; ++iter) {
      if constexpr (std::is_signed_v<typename std::iterator_traits<IT>::value_type>) {
        mWriter.Int(*iter);
      } else {
        mWriter.Uint(*iter);
      }
    }
    mWriter.EndArray();
  };

 private:
  static size_t id;

  std::ofstream mFileStream;
  rapidjson::OStreamWrapper mStreamWrapper;
  rapidjson::Writer<rapidjson::OStreamWrapper> mWriter;
};

inline size_t CTFJSONSerializer::id = 0;

} // namespace ctf
} // namespace o2

#endif // O2_CTF_JSONSERIALIZER_H