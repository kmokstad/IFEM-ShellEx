// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlPRGD.H"


FFlPRGD::FFlPRGD(int id) : FFlAttributeBase(id)
{
  this->addField(dependentDofs);
  dependentDofs = 0;
}


FFlPRGD::FFlPRGD(const FFlPRGD& obj) : FFlAttributeBase(obj)
{
  this->addField(dependentDofs);
  dependentDofs = obj.dependentDofs;
}


bool FFlPRGD::isIdentic(const FFlAttributeBase* otherAttrib) const
{
  const FFlPRGD* other = dynamic_cast<const FFlPRGD*>(otherAttrib);
  if (!other) return false;

  if (this->dependentDofs != other->dependentDofs) return false;

  return true;
}


void FFlPRGD::init()
{
  FFlPRGDTypeInfoSpec::instance()->setTypeName("PRGD");
  FFlPRGDTypeInfoSpec::instance()->setDescription("Rigid element properties");
  FFlPRGDTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::OTHER_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPRGDTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPRGD::create,int,FFlAttributeBase*&));
}
