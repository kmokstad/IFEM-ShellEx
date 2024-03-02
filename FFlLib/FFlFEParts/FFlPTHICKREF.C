// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlPTHICKREF.H"


FFlPTHICKREF::FFlPTHICKREF(int id) : FFlAttributeBase(id)
{
  this->addField(factor);
  factor = 0.5;
}


FFlPTHICKREF::FFlPTHICKREF(const FFlPTHICKREF& obj) : FFlAttributeBase(obj)
{
  this->addField(factor);
  factor = obj.factor;
}


void FFlPTHICKREF::init()
{
  FFlPTHICKREFTypeInfoSpec::instance()->setTypeName("PTHICKREF");
  FFlPTHICKREFTypeInfoSpec::instance()->setDescription("Strain coat heights");
  FFlPTHICKREFTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::STRC_PROP);

  FFlPTHICKREFAttributeSpec::instance()->addLegalAttribute("PTHICK", true);

  AttributeFactory::instance()->registerCreator
    (FFlPTHICKREFTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPTHICKREF::create,int,FFlAttributeBase*&));
}
