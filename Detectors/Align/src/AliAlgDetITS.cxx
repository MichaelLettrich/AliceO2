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

#include "ITSBase/GeometryTGeo.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"

#include "Align/AliAlgDetITS.h"
#include "Align/AliAlgVol.h"
#include "Align/AliAlgSensITS.h"
#include "Align/AliAlgSteer.h"
//#include "AliITSgeomTGeo.h"
//#include "AliGeomManager.h"
//#include "AliESDtrack.h"
//#include "AliCheb3DCalc.h"
#include <TMath.h>
#include <stdio.h>

using namespace TMath;
using namespace o2::align::AliAlgAux;

ClassImp(o2::align::AliAlgDetITS);

namespace o2
{
namespace align
{

using namespace o2::its;

const Char_t* AliAlgDetITS::fgkHitsSel[AliAlgDetITS::NITSSelTypes] =
  {"All"};

//____________________________________________
AliAlgDetITS::AliAlgDetITS(const char* title)
{
  // default c-tor
  SetNameTitle(AliAlgSteer::GetDetNameByDetID(AliAlgSteer::kITS), title);
  SetDetID(AliAlgSteer::kITS);
  SetITSSelPatternColl();
  SetITSSelPatternCosm();
}

//____________________________________________
AliAlgDetITS::~AliAlgDetITS()
{
  // d-tor
}

//____________________________________________
void AliAlgDetITS::DefineVolumes()
{
  // define ITS volumes
  //
  int labNS = GetDetLabel() + mMapping.getNChips(); // non-sensors labels will start from this
  auto volITS = AddVolume(new AliAlgVol(GeometryTGeo::composeSymNameITS(), GetDetLabel()));
  using MP = o2::itsmft::ChipMappingITS;

  auto addVolumeChip = [this](int ilr, int ist, int isst, int imd, int ich, AliAlgVol* parent) {
    auto volChip = AddVolume(new AliAlgSensITS(GeometryTGeo::composeSymNameChip(ilr, ist, isst, imd, ich), ilr, this->fSensors.GetEntriesFast()));
    volChip->SetParent(parent);
  };

  auto addVolumeModule = [&, this](int ilr, int ist, int isst, int imd, AliAlgVol* parent) {
    auto lrType = MP::layer2RUType(ilr);
    if (imd >= 0) {
      auto volMod = AddVolume(new AliAlgVol(GeometryTGeo::composeSymNameModule(ilr, ist, isst, imd), labNS++));
      volMod->SetParent(parent);
      parent = volMod;
    }
    for (int ich = 0; ich < MP::getNChipsPerModuleSB(lrType); ich++) {
      addVolumeChip(ilr, ist, isst, imd, ich, parent);
    }
  };

  auto addVolumeSStave = [&, this](int ilr, int ist, int isst, AliAlgVol* parent) {
    auto lrType = MP::layer2RUType(ilr);
    if (isst >= 0) {
      auto volSStave = AddVolume(new AliAlgVol(GeometryTGeo::composeSymNameHalfStave(ilr, ist, isst), labNS++));
      volSStave->SetParent(parent);
      parent = volSStave;
    }
    for (int imd = 0; imd < (lrType == MP::IB ? 1 : MP::getNModulesAlongStaveSB(lrType)); imd++) {
      addVolumeModule(ilr, ist, isst, imd, parent);
    }
  };

  auto addVolumeStave = [&, this](int ilr, int ist, AliAlgVol* parent) {
    auto lrType = MP::layer2RUType(ilr);
    auto volStave = AddVolume(new AliAlgVol(GeometryTGeo::composeSymNameStave(ilr, ist), labNS++));
    volStave->SetParent(parent);
    for (int isst = 0; isst < (lrType == MP::IB ? 1 : 2); isst++) { // substave volume exist only for layers of OB
      addVolumeSStave(ilr, ist, isst, volStave);
    }
  };

  auto addVolumeLayer = [&, this](int ilr, AliAlgVol* parent) {
    auto volLayer = AddVolume(new AliAlgVol(GeometryTGeo::composeSymNameLayer(ilr), labNS++));
    volLayer->SetParent(parent);
    for (int ist = 0; ist < MP::getNStavesOnLr(ilr); ist++) {
      addVolumeStave(ilr, ist, volLayer);
    }
  };

  for (int ilr = 0; ilr < MP::NLayers; ilr++) {
    addVolumeLayer(ilr, volITS);
  }
}

//____________________________________________
void AliAlgDetITS::Print(const Option_t* opt) const
{
  AliAlgDet::Print(opt);
  printf("Sel.pattern   Collisions: %7s | Cosmic: %7s\n",
         GetITSPattName(fITSPatt[kColl]), GetITSPattName(fITSPatt[kCosm]));
}

/* TODO RS templatize this
//____________________________________________
Bool_t AliAlgDetITS::AcceptTrack(const AliESDtrack* trc, Int_t trtype) const
{
  // test if detector had seed this track
  if (!CheckFlags(trc, trtype))
    return kFALSE;
  if (trc->GetNcls(0) < fNPointsSel[trtype])
    return kFALSE;
  if (!CheckHitPattern(trc, GetITSSelPattern(trtype)))
    return kFALSE;
  //
  return kTRUE;
}
*/

//____________________________________________
void AliAlgDetITS::SetAddErrorLr(int ilr, double sigY, double sigZ)
{
  // set syst. errors for specific layer
  for (int isn = GetNSensors(); isn--;) {
    AliAlgSensITS* sens = (AliAlgSensITS*)GetSensor(isn);
    if (sens->getLayer() > ilr) {
      continue;
    } else if (sens->getLayer() < ilr) {
      break;
    }
    sens->SetAddError(sigY, sigZ);
  }
}

//____________________________________________
void AliAlgDetITS::SetSkipLr(int ilr)
{
  // exclude sensor of the layer from alignment
  for (int isn = GetNSensors(); isn--;) {
    AliAlgSensITS* sens = (AliAlgSensITS*)GetSensor(isn);
    if (sens->getLayer() > ilr) {
      continue;
    } else if (sens->getLayer() < ilr) {
      break;
    }
    sens->SetSkip();
  }
}

//_________________________________________________
Bool_t AliAlgDetITS::CheckHitPattern(const o2::its::TrackITS& trc, ITSSel_t sel)
{
  // check if track hit pattern is ok
  // TODO RS
  /*
  switch (sel) {
    case kSPDBoth:
      if (!trc->HasPointOnITSLayer(0) || !trc->HasPointOnITSLayer(1))
        return kFALSE;
      break;
    case kSPDAny:
      if (!trc->HasPointOnITSLayer(0) && !trc->HasPointOnITSLayer(1))
        return kFALSE;
      break;
    case kSPD0:
      if (!trc->HasPointOnITSLayer(0))
        return kFALSE;
      break;
    case kSPD1:
      if (!trc->HasPointOnITSLayer(1))
        return kFALSE;
      break;
    default:
      break;
  }
  */
  return kTRUE;
}

/*
// TODO RS
//_________________________________________________
void AliAlgDetITS::UpdatePointByTrackInfo(AliAlgPoint* pnt, const AliExternalTrackParam* t) const
{
  // update point using specific error parameterization
  // the track must be in the detector tracking frame
  const AliAlgSens* sens = pnt->GetSensor();
  int vid = sens->GetVolID();
  int lr = AliGeomManager::VolUIDToLayer(vid) - 1;
  double angPol = ATan(t->GetTgl());
  double angAz = ASin(t->GetSnp());
  double errY, errZ;
  GetErrorParamAngle(lr, angPol, angAz, errY, errZ);
  const double* sysE = sens->GetAddError(); // additional syst error
  //
  pnt->SetYZErrTracking(errY * errY + sysE[0] * sysE[0], 0, errZ * errZ + sysE[1] * sysE[1]);
  pnt->Init();
  //
}

*/

} // namespace align
} // namespace o2
