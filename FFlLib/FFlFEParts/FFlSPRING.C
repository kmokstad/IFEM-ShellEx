// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlSPRING.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


void FFlSPRING::init()
{
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlSPRING>;
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlSPRING>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlSPRING>;

  TypeInfoSpec::instance()->setTypeName("SPRING");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::OTHER_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlSPRING::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PSPRING");

  ElementTopSpec::instance()->setNodeCount(2);
  ElementTopSpec::instance()->setNodeDOFs(3);
}


void FFlRSPRING::init()
{
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlRSPRING>;
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlRSPRING>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlRSPRING>;

  TypeInfoSpec::instance()->setTypeName("RSPRING");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::OTHER_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlRSPRING::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PSPRING");

  ElementTopSpec::instance()->setNodeCount(2);
  ElementTopSpec::instance()->setNodeDOFs(6);
}

#ifdef FF_NAMESPACE
} // namespace
#endif
