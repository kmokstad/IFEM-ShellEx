// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlSPRING.H"


void FFlSPRING::init()
{
  FFlSPRINGTypeInfoSpec::instance()->setTypeName("SPRING");
  FFlSPRINGTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::OTHER_ELM);

  ElementFactory::instance()->registerCreator(FFlSPRINGTypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlSPRING::create,int,FFlElementBase*&));

  FFlSPRINGAttributeSpec::instance()->addLegalAttribute("PSPRING");
  FFlSPRINGElementTopSpec::instance()->setNodeCount(2);
  FFlSPRINGElementTopSpec::instance()->setNodeDOFs(3);
}


void FFlRSPRING::init()
{
  FFlRSPRINGTypeInfoSpec::instance()->setTypeName("RSPRING");
  FFlRSPRINGTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::OTHER_ELM);

  ElementFactory::instance()->registerCreator(FFlRSPRINGTypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlRSPRING::create,int,FFlElementBase*&));

  FFlRSPRINGAttributeSpec::instance()->addLegalAttribute("PSPRING");
  FFlRSPRINGElementTopSpec::instance()->setNodeCount(2);
  FFlRSPRINGElementTopSpec::instance()->setNodeDOFs(6);
}
