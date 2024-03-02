// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlPTHICK.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"


FFlPTHICK::FFlPTHICK(int id) : FFlAttributeBase(id)
{
  this->addField(thickness);
  thickness = 0.1;
}


FFlPTHICK::FFlPTHICK(const FFlPTHICK& obj) : FFlAttributeBase(obj)
{
  this->addField(thickness);
  thickness = obj.thickness;
}


void FFlPTHICK::convertUnits(const FFaUnitCalculator* convCal)
{
  // Round to 10 significant digits
  convCal->convert(thickness.data(),"LENGTH",10);
}


void FFlPTHICK::init()
{
  FFlPTHICKTypeInfoSpec::instance()->setTypeName("PTHICK");
  FFlPTHICKTypeInfoSpec::instance()->setDescription("Shell thicknesses");
  FFlPTHICKTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPTHICKTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPTHICK::create,int,FFlAttributeBase*&));
}
