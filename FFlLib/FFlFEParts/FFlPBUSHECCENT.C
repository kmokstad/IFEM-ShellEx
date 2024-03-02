// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPBUSHECCENT.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"


FFlPBUSHECCENT::FFlPBUSHECCENT(int id) : FFlAttributeBase(id)
{
  this->addField(offset);
}


FFlPBUSHECCENT::FFlPBUSHECCENT(const FFlPBUSHECCENT& ob) : FFlAttributeBase(ob)
{
  this->addField(offset);

  offset = ob.offset;
}


void FFlPBUSHECCENT::convertUnits(const FFaUnitCalculator* convCal)
{
  convCal->convert(offset.data(),"LENGTH");
}


void FFlPBUSHECCENT::init()
{
  FFlPBUSHECCENTTypeInfoSpec::instance()->setTypeName("PBUSHECCENT");
  FFlPBUSHECCENTTypeInfoSpec::instance()->setDescription("Spring element eccentricities");
  FFlPBUSHECCENTTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPBUSHECCENTTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPBUSHECCENT::create,int,FFlAttributeBase*&));
}
