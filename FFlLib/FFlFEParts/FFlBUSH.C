// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlBUSH.H"
#include "FFlLib/FFlFEParts/FFlPBUSHECCENT.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


void FFlBUSH::init()
{
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlBUSH>;
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlBUSH>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlBUSH>;

  TypeInfoSpec::instance()->setTypeName("BUSH");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::OTHER_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlBUSH::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PBUSHCOEFF",false);
  AttributeSpec::instance()->addLegalAttribute("PBUSHECCENT",false);
  AttributeSpec::instance()->addLegalAttribute("PBUSHORIENT",false);
  AttributeSpec::instance()->addLegalAttribute("PORIENT",false);
  AttributeSpec::instance()->addLegalAttribute("PCOORDSYS",false);

  ElementTopSpec::instance()->setNodeCount(2);
  ElementTopSpec::instance()->setNodeDOFs(6);
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

#ifdef FF_NAMESPACE
} // namespace
#endif
