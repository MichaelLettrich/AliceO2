// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   AliAlgDetITS.h
/// @author ruben.shahoyan@cern.ch, michael.lettrich@cern.ch
/// @since  2021-02-01
/// @brief  ITS detector wrapper

#ifndef ALIALGDETITS_H
#define ALIALGDETITS_H

#include "Align/AliAlgDet.h"
#include "Align/AliAlgAux.h"
#include "DataFormatsITS/TrackITS.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"

namespace o2
{
namespace align
{

class AliAlgDetITS : public AliAlgDet
{
 public:
  //
  enum ITSSel_t { All,
                  NITSSelTypes };

  AliAlgDetITS(const char* title = "");
  virtual ~AliAlgDetITS();
  //
  virtual void DefineVolumes();
  //
  //  Bool_t AcceptTrack(const AliESDtrack* trc, Int_t trtype) const; // TODO RS templatize this

  void SetAddErrorLr(int ilr, double sigY, double sigZ);
  void SetSkipLr(int ilr);
  //
  /// virtual void UpdatePointByTrackInfo(AliAlgPoint* pnt, const AliExternalTrackParam* t) const; // TODO RS
  void SetITSSelPattern(AliAlgAux::TrackType tp, ITSSel_t sel) { fITSPatt[tp] = sel; }
  void SetITSSelPatternColl(ITSSel_t sel = All) { SetITSSelPattern(AliAlgAux::kColl, sel); }
  void SetITSSelPatternCosm(ITSSel_t sel = All) { SetITSSelPattern(AliAlgAux::kCosm, sel); }

  Int_t GetITSSelPattern(AliAlgAux::TrackType tp) const { return fITSPatt[tp]; }
  Int_t GetITSSelPatternColl() const { return fITSPatt[AliAlgAux::kColl]; }
  Int_t GetITSSelPatternCosm() const { return fITSPatt[AliAlgAux::kCosm]; }

  auto& getChipMapping() const { return mMapping; }

  //
  virtual void Print(const Option_t* opt = "") const;
  //
  static Bool_t CheckHitPattern(const o2::its::TrackITS& trc, ITSSel_t sel);
  static const char* GetITSPattName(ITSSel_t sel) { return fgkHitsSel[sel]; }
  //
 protected:
  //
  ITSSel_t fITSPatt[AliAlgAux::kNTrackTypes]; // ITS hits selection pattern for coll/cosm tracks
  //
  static const Char_t* fgkHitsSel[NITSSelTypes]; // ITS selection names // TODO RS Change to constexpr string_view
  //
  o2::itsmft::ChipMappingITS mMapping;

  ClassDef(AliAlgDetITS, 1);
};
} // namespace align
} // namespace o2
#endif
