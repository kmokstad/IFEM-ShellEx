// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPBUSHECCENT.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


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
  using TypeInfoSpec = FFaSingelton<FFlTypeInfoSpec,FFlPBUSHECCENT>;

  TypeInfoSpec::instance()->setTypeName("PBUSHECCENT");
  TypeInfoSpec::instance()->setDescription("Spring element eccentricities");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPBUSHECCENT::create,int,FFlAttributeBase*&));
}

#ifdef FF_NAMESPACE
} // namespace
#endif
