// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPMATSHELL.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"


FFlPMATSHELL::FFlPMATSHELL(int id) : FFlAttributeBase(id)
{
  this->addField(E1);
  this->addField(E2);
  this->addField(NU12);
  this->addField(G12);
  this->addField(G1Z);
  this->addField(G2Z);
  this->addField(materialDensity);

  // default values - Typical for steel
  double E            = 205.00E+9;
  double shearModule  =  80.00E+9;
  double nu           =   0.29E+0;

  E1 = E;
  E2 = E;
  NU12 = nu;
  G12 = shearModule;
  G1Z = shearModule;
  G2Z = shearModule;
  materialDensity =   7.85E+3;
}


FFlPMATSHELL::FFlPMATSHELL(const FFlPMATSHELL& obj) : FFlAttributeBase(obj)
{
  this->addField(E1);
  this->addField(E2);
  this->addField(NU12);
  this->addField(G12);
  this->addField(G1Z);
  this->addField(G2Z);
  this->addField(materialDensity);

  E1 = obj.E1;
  E2 = obj.E2;
  NU12 = obj.NU12;
  G12 = obj.G12;
  G1Z = obj.G1Z;
  G2Z = obj.G2Z;
  materialDensity = obj.materialDensity;
}


void FFlPMATSHELL::convertUnits(const FFaUnitCalculator* convCal)
{
  // Round to 10 significant digits
  convCal->convert(E1.data(),"FORCE/AREA",10);
  convCal->convert(E2.data(),"FORCE/AREA",10);
  convCal->convert(NU12.data(),"FORCE/AREA",10);
  convCal->convert(G12.data(),"FORCE/AREA",10);
  convCal->convert(G1Z.data(),"FORCE/AREA",10);
  convCal->convert(G2Z.data(),"FORCE/AREA",10);
  convCal->convert(materialDensity.data(),"MASS/VOLUME",10);
}


void FFlPMATSHELL::init()
{
  FFlPMATSHELLTypeInfoSpec::instance()->setTypeName("PMATSHELL");
  FFlPMATSHELLTypeInfoSpec::instance()->setDescription("Material properties");
  FFlPMATSHELLTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::MATERIAL_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPMATSHELLTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPMATSHELL::create,int,FFlAttributeBase*&));
}
