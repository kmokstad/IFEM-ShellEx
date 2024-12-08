// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_ALL_FE_PARTS_H
#define FFL_ALL_FE_PARTS_H
#ifndef FFL_INIT_ONLY

#include "FFlNode.H"
#include "FFlRGD.H"
#include "FFlRBAR.H"
#include "FFlWAVGM.H"
#include "FFlCMASS.H"
#include "FFlSPRING.H"
#include "FFlBUSH.H"
#include "FFlBeams.H"
#include "FFlShells.H"
#include "FFlSolids.H"
#include "FFlLoads.H"
#include "FFlPMASS.H"
#include "FFlPSPRING.H"
#include "FFlPBUSHCOEFF.H"
#include "FFlPBUSHECCENT.H"
#include "FFlPORIENT.H"
#include "FFlPBEAMECCENT.H"
#include "FFlPBEAMSECTION.H"
#include "FFlPBEAMPIN.H"
#include "FFlPEFFLENGTH.H"
#include "FFlPSPRING.H"
#include "FFlPTHICK.H"
#include "FFlPCOMP.H"
#include "FFlPNSM.H"
#include "FFlPRGD.H"
#include "FFlPWAVGM.H"
#include "FFlPRBAR.H"
#include "FFlPMAT.H"
#include "FFlPCOORDSYS.H"
#ifdef FT_USE_STRAINCOAT
#include "FFlSTRCoat.H"
#include "FFlPSTRC.H"
#include "FFlPHEIGHT.H"
#include "FFlPTHICKREF.H"
#include "FFlPFATIGUE.H"
#endif
#ifdef FT_USE_VISUALS
#include "FFlVAppearance.H"
#include "FFlVDetail.H"
#endif
#include "FFlLib/FFlGroup.H"

#endif


namespace FFl
{
  //! \brief Initializes all element singleton objects.
  void initAllElements();
  //! \brief Releases all heap-allocated memory related to singleton objects.
  void releaseAllElements();
}

#endif
