// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPSPRING.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"


FFlPSPRING::FFlPSPRING(int id) : FFlAttributeBase(id)
{
  this->addField(type);
  type = TRANS_SPRING;
  for (int i = 0; i < PSPRING_size; i++)
  {
    this->addField(K[i]);
    K[i] = 0.0;
  }
}


FFlPSPRING::FFlPSPRING(const FFlPSPRING& obj) : FFlAttributeBase(obj)
{
  this->addField(type);
  type = obj.type;
  for (int i = 0; i < PSPRING_size; i++)
  {
    this->addField(K[i]);
    K[i] = obj.K[i];
  }
}


void FFlPSPRING::convertUnits(const FFaUnitCalculator* convCal)
{
  const std::string unit = (type.data() == TRANS_SPRING ? "FORCE/LENGTH" : "FORCE*LENGTH");

  // Round to 10 significant digits
  for (int i = 0; i < PSPRING_size; i++)
    convCal->convert(K[i].data(),unit,10);
}


void FFlPSPRING::init()
{
  FFlPSPRINGTypeInfoSpec::instance()->setTypeName("PSPRING");
  FFlPSPRINGTypeInfoSpec::instance()->setDescription("Spring element properties");
  FFlPSPRINGTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::OTHER_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPSPRINGTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPSPRING::create,int,FFlAttributeBase*&));
}
