// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPSPRING.H"
#include "FFlLib/FFlFEParts/FFlPEFFLENGTH.H"

#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"


FFlPEFFLENGTH::FFlPEFFLENGTH(int id) : FFlAttributeBase(id)
{
  this->addField(length);
  length = 0.0;
}


FFlPEFFLENGTH::FFlPEFFLENGTH(const FFlPEFFLENGTH& obj) : FFlAttributeBase(obj)
{
  this->addField(length);
  length = obj.length;
}


void FFlPEFFLENGTH::convertUnits(const FFaUnitCalculator* convCal)
{
  convCal->convert(length.data(), "LENGTH");
}


bool FFlPEFFLENGTH::isIdentic(const FFlAttributeBase* otherAttrib) const
{
  const FFlPEFFLENGTH* other = dynamic_cast<const FFlPEFFLENGTH*>(otherAttrib);
  return other ? (this->length == other->length) : false;
}


void FFlPEFFLENGTH::init()
{
  FFlPEFFLENGTHTypeInfoSpec::instance()->setTypeName("PEFFLENGTH");
  FFlPEFFLENGTHTypeInfoSpec::instance()->setDescription("Effective beam lengths");
  FFlPEFFLENGTHTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPEFFLENGTHTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPEFFLENGTH::create,int,FFlAttributeBase*&));
}
