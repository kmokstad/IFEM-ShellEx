// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlPSTRC.H"


FFlPSTRC::FFlPSTRC(int id) : FFlAttributeBase(id)
{
  this->addField(name);
}


FFlPSTRC::FFlPSTRC(const FFlPSTRC& obj) : FFlAttributeBase(obj)
{
  this->addField(name);
  name = obj.name;
}


void FFlPSTRC::init()
{
  FFlPSTRCTypeInfoSpec::instance()->setTypeName("PSTRC");
  FFlPSTRCTypeInfoSpec::instance()->setDescription("Strain coat properties");
  FFlPSTRCTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::STRC_PROP);

  FFlPSTRCAttributeSpec::instance()->addLegalAttribute("PMAT"     , false);
  FFlPSTRCAttributeSpec::instance()->addLegalAttribute("PTHICKREF", false);
  FFlPSTRCAttributeSpec::instance()->addLegalAttribute("PHEIGHT"  , false);

  AttributeFactory::instance()->registerCreator
    (FFlPSTRCTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPSTRC::create,int,FFlAttributeBase*&));
}
