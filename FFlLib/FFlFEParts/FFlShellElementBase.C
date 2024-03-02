// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlShellElementBase.H"
#include "FFlLib/FFlFEParts/FFlPTHICK.H"
#include "FFlLib/FFlFEParts/FFlPCOMP.H"
#include "FFlLib/FFlFEParts/FFlPMAT.H"
#include "FFlLib/FFlFEParts/FFlPNSM.H"


double FFlShellElementBase::offPlaneTol = 0.1;


double FFlShellElementBase::getThickness() const
{
  FFlPTHICK* pthk = dynamic_cast<FFlPTHICK*>(this->getAttribute("PTHICK"));
  if (pthk)
    return pthk->thickness.getValue();

  FFlPCOMP* pcomp = dynamic_cast<FFlPCOMP*>(this->getAttribute("PCOMP"));
  if (pcomp)
    return -2.0*pcomp->Z0.getValue(); //TODO,bh: This is quick and dirty

  return 0.0; // Should not happen
}


double FFlShellElementBase::getMassDensity() const
{
  FFlPMAT* pmat = dynamic_cast<FFlPMAT*>(this->getAttribute("PMAT"));
  double rho = pmat ? pmat->materialDensity.getValue() : 0.0;

  // Check for non-structural mass
  FFlPNSM* pnsm = dynamic_cast<FFlPNSM*>(this->getAttribute("PNSM"));
  if (!pnsm) return rho;

  // Shell thickness
  double th = this->getThickness();
  if (th < 1.0e-16) return rho; // avoid division by zero

  // Modify the mass density to account for the non-structural mass
  return rho + pnsm->NSM.getValue() / th;
}
