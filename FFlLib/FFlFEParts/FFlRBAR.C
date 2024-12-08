// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlRBAR.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


void FFlRBAR::init()
{
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlRBAR>;
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlRBAR>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlRBAR>;

  TypeInfoSpec::instance()->setTypeName("RBAR");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::CONSTRAINT_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlRBAR::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PRBAR");

  ElementTopSpec::instance()->setNodeCount(2);
  ElementTopSpec::instance()->setNodeDOFs(6);
  ElementTopSpec::instance()->setSlaveStatus(true); // Both nodes have slave DOFs

  ElementTopSpec::instance()->addExplicitEdge(EdgeType(1,2));
  ElementTopSpec::instance()->myExplEdgePattern = 0xf0f0;
}

#ifdef FF_NAMESPACE
} // namespace
#endif
