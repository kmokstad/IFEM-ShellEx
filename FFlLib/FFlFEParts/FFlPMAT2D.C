// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlPMAT2D.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"


//TODO,bh: use lower triangle here. Similar to other matrices
FFlPMAT2D::FFlPMAT2D(int id) : FFlAttributeBase(id)
{
  this->addField(C11);
  this->addField(C12);
  this->addField(C13);
  this->addField(C22);
  this->addField(C23);
  this->addField(C33);
  this->addField(materialDensity);

  // default values - Typical for steel
  double E  = 205.00E+9;
  double nu =   0.29E+0;
  double f1 = E / ( 1.0+nu*nu );
  double f2 = f1 * nu;
  double f3 = f1 * (1.0-nu) / 2.0;

  C11 = f1;
  C12 = f2;
  C13 = 0.0;
  C22 = f1;
  C23 = 0.0;
  C33 = f3;
  materialDensity =   7.85E+3;
}


FFlPMAT2D::FFlPMAT2D(const FFlPMAT2D& obj) : FFlAttributeBase(obj)
{
  this->addField(C11);
  this->addField(C12);
  this->addField(C13);
  this->addField(C22);
  this->addField(C23);
  this->addField(C33);
  this->addField(materialDensity);

  C11 = obj.C11;
  C12 = obj.C12;
  C13 = obj.C13;
  C22 = obj.C22;
  C23 = obj.C23;
  C33 = obj.C33;
  materialDensity = obj.materialDensity;
}


void FFlPMAT2D::convertUnits(const FFaUnitCalculator* convCal)
{
  // Round to 10 significant digits
  convCal->convert(C11.data(),"FORCE/AREA",10);
  convCal->convert(C12.data(),"FORCE/AREA",10);
  convCal->convert(C13.data(),"FORCE/AREA",10);
  convCal->convert(C22.data(),"FORCE/AREA",10);
  convCal->convert(C23.data(),"FORCE/AREA",10);
  convCal->convert(C33.data(),"FORCE/AREA",10);
  convCal->convert(materialDensity.data(),"MASS/VOLUME",10);
}


void FFlPMAT2D::init()
{
  FFlPMAT2DTypeInfoSpec::instance()->setTypeName("PMAT2D");
  FFlPMAT2DTypeInfoSpec::instance()->setDescription("Material properties");
  FFlPMAT2DTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::MATERIAL_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPMAT2DTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPMAT2D::create,int,FFlAttributeBase*&));
}
