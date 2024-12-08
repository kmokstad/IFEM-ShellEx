// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPCOORDSYS.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


FFlPCOORDSYS::FFlPCOORDSYS(int id) : FFlAttributeBase(id)
{
  this->addField(Origo);
  this->addField(Zaxis);
  this->addField(XZpnt);

  Origo = FaVec3(0,0,0);
  Zaxis = FaVec3(0,0,1);
  XZpnt = FaVec3(1,0,0);
}


FFlPCOORDSYS::FFlPCOORDSYS(const FFlPCOORDSYS& ob) : FFlAttributeBase(ob)
{
  this->addField(Origo);
  this->addField(Zaxis);
  this->addField(XZpnt);

  Origo = ob.Origo;
  Zaxis = ob.Zaxis;
  XZpnt = ob.XZpnt;
}


void FFlPCOORDSYS::convertUnits(const FFaUnitCalculator* convCal)
{
  convCal->convert(Origo.data(),"LENGTH");
  convCal->convert(Zaxis.data(),"LENGTH");
  convCal->convert(XZpnt.data(),"LENGTH");
}


void FFlPCOORDSYS::init()
{
  using TypeInfoSpec = FFaSingelton<FFlTypeInfoSpec,FFlPCOORDSYS>;

  TypeInfoSpec::instance()->setTypeName("PCOORDSYS");
  TypeInfoSpec::instance()->setDescription("Local coordinate systems");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPCOORDSYS::create,int,FFlAttributeBase*&));
}

#ifdef FF_NAMESPACE
} // namespace
#endif
