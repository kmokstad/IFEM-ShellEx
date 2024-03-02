// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlPBEAMPIN.H"


FFlPBEAMPIN::FFlPBEAMPIN(int id) : FFlAttributeBase(id)
{
  this->addField(PA);
  this->addField(PB);

  PA = 0;
  PB = 0;
}


FFlPBEAMPIN::FFlPBEAMPIN(const FFlPBEAMPIN& obj) : FFlAttributeBase(obj)
{
  this->addField(PA);
  this->addField(PB);

  PA = obj.PA;
  PB = obj.PB;
}


bool FFlPBEAMPIN::isIdentic(const FFlAttributeBase* otherAttrib) const
{
  const FFlPBEAMPIN* other = dynamic_cast<const FFlPBEAMPIN*>(otherAttrib);
  if (!other) return false;

  if (this->PA != other->PA) return false;
  if (this->PB != other->PB) return false;

  return true;
}


void FFlPBEAMPIN::init()
{
  FFlPBEAMPINTypeInfoSpec::instance()->setTypeName("PBEAMPIN");
  FFlPBEAMPINTypeInfoSpec::instance()->setDescription("Beam pin flags");
  FFlPBEAMPINTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::OTHER_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPBEAMPINTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPBEAMPIN::create,int,FFlAttributeBase*&));
}


