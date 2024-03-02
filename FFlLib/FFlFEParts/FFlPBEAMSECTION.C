// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPBEAMSECTION.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"


FFlPBEAMSECTION::FFlPBEAMSECTION(int id) : FFlAttributeBase(id)
{
  this->setupFields();

  crossSectionArea = 1.963E-3;
  Iy = Iz = 3.068E-7;
  It = 6.136E-7;
  Kxy = Kxz = Sy = Sz = 0.0;
}


FFlPBEAMSECTION::FFlPBEAMSECTION(const FFlPBEAMSECTION& obj)
  : FFlAttributeBase(obj)
{
  this->setupFields();

  crossSectionArea = obj.crossSectionArea;
  Iy = obj.Iy ;
  Iz = obj.Iz;
  It = obj.It;
  Kxy = obj.Kxy;
  Kxz = obj.Kxz;
  Sy = obj.Sy;
  Sz = obj.Sz;
}


void FFlPBEAMSECTION::convertUnits(const FFaUnitCalculator* convCal)
{
  // Round to 10 significant digits
  convCal->convert(crossSectionArea.data(),"AREA",10);
  convCal->convert(Iy.data(),"LENGTH^4",10);
  convCal->convert(Iz.data(),"LENGTH^4",10);
  convCal->convert(It.data(),"LENGTH^4",10);
  convCal->convert(Sy.data(),"LENGTH",10);
  convCal->convert(Sz.data(),"LENGTH",10);
}


void FFlPBEAMSECTION::setupFields()
{
  this->addField(crossSectionArea);
  this->addField(Iy);
  this->addField(Iz);
  this->addField(It);
  this->addField(Kxy);
  this->addField(Kxz);
  this->addField(Sy);
  this->addField(Sz);
}


void FFlPBEAMSECTION::init()
{
  FFlPBEAMSECTIONTypeInfoSpec::instance()->setTypeName("PBEAMSECTION");
  FFlPBEAMSECTIONTypeInfoSpec::instance()->setDescription("Beam cross sections");
  FFlPBEAMSECTIONTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPBEAMSECTIONTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPBEAMSECTION::create,int,FFlAttributeBase*&));
}
