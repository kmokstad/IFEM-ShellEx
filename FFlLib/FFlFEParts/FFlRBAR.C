// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlRBAR.H"


void FFlRBAR::init()
{
  FFlRBARTypeInfoSpec::instance()->setTypeName("RBAR");
  FFlRBARTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::CONSTRAINT_ELM);

  ElementFactory::instance()->registerCreator(FFlRBARTypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlRBAR::create,int,FFlElementBase*&));

  FFlRBARAttributeSpec::instance()->addLegalAttribute("PRBAR");

  FFlRBARElementTopSpec::instance()->setNodeCount(2);
  FFlRBARElementTopSpec::instance()->setNodeDOFs(6);
  FFlRBARElementTopSpec::instance()->setSlaveStatus(true); // Both nodes have slave DOFs

  FFlRBARElementTopSpec::instance()->addExplicitEdge(EdgeType(1,2));
  FFlRBARElementTopSpec::instance()->myExplEdgePattern = 0xf0f0;
}
