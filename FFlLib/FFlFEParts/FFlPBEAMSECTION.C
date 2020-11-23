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

  // Default values: Circular bar with radius R=0.025
  crossSectionArea.setValue(0.001963);
  Iy.setValue(3.068e-7);
  Iz.setValue(3.068e-7);
  It.setValue(6.136e-7);
}


FFlPBEAMSECTION::FFlPBEAMSECTION(const FFlPBEAMSECTION& obj)
  : FFlAttributeBase(obj)
{
  this->setupFields();

  crossSectionArea.setValue(obj.crossSectionArea.getValue());
  Iy.setValue(obj.Iy.getValue());
  Iz.setValue(obj.Iz.getValue());
  It.setValue(obj.It.getValue());
  Kxy.setValue(obj.Kxy.getValue());
  Kxz.setValue(obj.Kxz.getValue());
  Sy.setValue(obj.Sy.getValue());
  Sz.setValue(obj.Sz.getValue());
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
