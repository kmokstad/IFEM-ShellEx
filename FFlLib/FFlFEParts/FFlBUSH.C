// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlBUSH.H"


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
