// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlBUSH.H"
#include "FFlLib/FFlFEParts/FFlPBUSHECCENT.H"


void FFlBUSH::init()
{
  FFlBUSHTypeInfoSpec::instance()->setTypeName("BUSH");
  FFlBUSHTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::OTHER_ELM);

  ElementFactory::instance()->registerCreator(FFlBUSHTypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlBUSH::create,int,FFlElementBase*&));

  FFlBUSHAttributeSpec::instance()->addLegalAttribute("PBUSHCOEFF", false);
  FFlBUSHAttributeSpec::instance()->addLegalAttribute("PBUSHECCENT", false);
  FFlBUSHAttributeSpec::instance()->addLegalAttribute("PBUSHORIENT", false);
  FFlBUSHAttributeSpec::instance()->addLegalAttribute("PORIENT", false);
  FFlBUSHAttributeSpec::instance()->addLegalAttribute("PCOORDSYS", false);

  FFlBUSHElementTopSpec::instance()->setNodeCount(2);
  FFlBUSHElementTopSpec::instance()->setNodeDOFs(6);
}


int FFlBUSH::getNodalCoor(double* X, double* Y, double* Z) const
{
  int ierr = this->FFlElementBase::getNodalCoor(X,Y,Z);
  if (ierr < 0) return ierr;

  X[2] = X[1];
  Y[2] = Y[1];
  Z[2] = Z[1];
  X[1] = X[0];
  Y[1] = Y[0];
  Z[1] = Z[0];

  // Get bushing eccentricity, if any
  FFlPBUSHECCENT* curEcc = dynamic_cast<FFlPBUSHECCENT*>(this->getAttribute("PBUSHECCENT"));
  if (!curEcc) return ierr;

  FaVec3 e1 = curEcc->offset.getValue();
  X[0] += e1[0];
  Y[0] += e1[1];
  Z[0] += e1[2];

  return ierr;
}
