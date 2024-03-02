// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlPHEIGHT.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"


FFlPHEIGHT::FFlPHEIGHT(int id) : FFlAttributeBase(id)
{
  this->addField(height);
  height = 0.1;
}


FFlPHEIGHT::FFlPHEIGHT(const FFlPHEIGHT& obj)
  : FFlAttributeBase(obj)
{
  this->addField(height);
  height = obj.height;
}


void FFlPHEIGHT::convertUnits(const FFaUnitCalculator* convCal)
{
  convCal->convert(height.data(),"LENGTH");
}


void FFlPHEIGHT::init()
{
  FFlPHEIGHTTypeInfoSpec::instance()->setTypeName("PHEIGHT");
  FFlPHEIGHTTypeInfoSpec::instance()->setDescription("Strain coat heights");
  FFlPHEIGHTTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPHEIGHTTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPHEIGHT::create,int,FFlAttributeBase*&));

}
