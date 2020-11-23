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
  type.setValue(TRANS_SPRING);

  for (FFlField<double>& field : K)
    this->addField(field);
}


FFlPSPRING::FFlPSPRING(const FFlPSPRING& obj) : FFlAttributeBase(obj)
{
  this->addField(type);
  type.setValue(obj.type.getValue());

  for (size_t i = 0; i < K.size(); i++)
  {
    this->addField(K[i]);
    K[i].setValue(obj.K[i].getValue());
  }
}


void FFlPSPRING::convertUnits(const FFaUnitCalculator* convCal)
{
  const char* unit = type.getValue() == TRANS_SPRING ? "FORCE/LENGTH" : "FORCE*LENGTH";

  // Round to 10 significant digits
  for (FFlField<double>& field : K)
    convCal->convert(field.data(),unit,10);
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
