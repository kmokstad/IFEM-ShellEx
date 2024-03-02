// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlPBUSHCOEFF.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"


FFlPBUSHCOEFF::FFlPBUSHCOEFF(int id) : FFlAttributeBase(id)
{
  for (int i = 0; i < 6; i++)
  {
    this->addField(K[i]);
    K[i] = 0.0;
  }
}


FFlPBUSHCOEFF::FFlPBUSHCOEFF(const FFlPBUSHCOEFF& obj) : FFlAttributeBase(obj)
{
  for (int i = 0; i < 6; i++)
  {
    this->addField(K[i]);
    K[i] = obj.K[i];
  }
}


void FFlPBUSHCOEFF::convertUnits(const FFaUnitCalculator* convCal)
{
  int i; // Round to 10 significant digits
  for (i = 0; i < 3; i++) convCal->convert(K[i].data(),"FORCE/LENGTH",10);
  for (i = 3; i < 6; i++) convCal->convert(K[i].data(),"FORCE*LENGTH",10);
}


void FFlPBUSHCOEFF::init()
{
  FFlPBUSHCOEFFTypeInfoSpec::instance()->setTypeName("PBUSHCOEFF");
  FFlPBUSHCOEFFTypeInfoSpec::instance()->setDescription("Spring element coefficients");

  FFlPBUSHCOEFFTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::OTHER_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPBUSHCOEFFTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPBUSHCOEFF::create,int,FFlAttributeBase*&));
}
