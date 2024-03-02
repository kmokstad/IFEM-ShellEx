// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPMAT.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"


FFlPMAT::FFlPMAT(int id) : FFlAttributeBase(id)
{
  this->addField(youngsModule);
  this->addField(shearModule);
  this->addField(poissonsRatio);
  this->addField(materialDensity);

  // default values - Typical for steel
  youngsModule    = 205.00E+9;
  shearModule     =  80.00E+9;
  poissonsRatio   =   0.29E+0;
  materialDensity =   7.85E+3;
}


FFlPMAT::FFlPMAT(const FFlPMAT& obj) : FFlAttributeBase(obj)
{
  this->addField(youngsModule);
  this->addField(shearModule);
  this->addField(poissonsRatio);
  this->addField(materialDensity);

  youngsModule    = obj.youngsModule;
  shearModule     = obj.shearModule;
  poissonsRatio   = obj.poissonsRatio;
  materialDensity = obj.materialDensity;
}


void FFlPMAT::convertUnits(const FFaUnitCalculator* convCal)
{
  // Round to 10 significant digits
  convCal->convert(youngsModule.data(),"FORCE/AREA",10);
  convCal->convert(shearModule.data(),"FORCE/AREA",10);
  convCal->convert(materialDensity.data(),"MASS/VOLUME",10);
}


void FFlPMAT::init()
{
  FFlPMATTypeInfoSpec::instance()->setTypeName("PMAT");
  FFlPMATTypeInfoSpec::instance()->setDescription("Material properties");
  FFlPMATTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::MATERIAL_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPMATTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPMAT::create,int,FFlAttributeBase*&));
}
