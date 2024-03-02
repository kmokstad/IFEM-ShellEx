// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlPMAT3D.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"

//TODO,bh: use lower triangle here. Similar to other matrices

FFlPMAT3D::FFlPMAT3D(int id) : FFlAttributeBase(id)
{
  this->addField(C11);
  this->addField(C12);
  this->addField(C13);
  this->addField(C14);
  this->addField(C15);
  this->addField(C16);
  this->addField(C22);
  this->addField(C23);
  this->addField(C24);
  this->addField(C25);
  this->addField(C26);
  this->addField(C33);
  this->addField(C34);
  this->addField(C35);
  this->addField(C36);
  this->addField(C44);
  this->addField(C45);
  this->addField(C46);
  this->addField(C55);
  this->addField(C56);
  this->addField(C66);
  this->addField(materialDensity);

  // default values - Typical for steel
  double E  = 205.00E+9;
  double nu =  0.29E+0;
  double f1 = E * (1.0-nu) / ( (1.0+nu) * (1.0-2.0*nu) );
  double f2 = f1 * nu / (1.0-nu);
  double f3 = f1 * (1.0-1.0*nu) / ( 2.0 * (1.0-nu) );
  C11 = f1;
  C12 = f2;
  C13 = f2;
  C14 = 0.0;
  C15 = 0.0;
  C16 = 0.0;
  C22 = f1;
  C23 = f2;
  C24 = 0.0;
  C25 = 0.0;
  C26 = 0.0;
  C33 = f1;
  C34 = 0.0;
  C35 = 0.0;
  C36 = 0.0;
  C44 = f3;
  C45 = 0.0;
  C46 = 0.0;
  C55 = f3;
  C56 = 0.0;
  C66 = f3;
  materialDensity =   7.85E+3;
}


FFlPMAT3D::FFlPMAT3D(const FFlPMAT3D& obj) : FFlAttributeBase(obj)
{
  this->addField(C11);
  this->addField(C12);
  this->addField(C13);
  this->addField(C14);
  this->addField(C15);
  this->addField(C16);
  this->addField(C22);
  this->addField(C23);
  this->addField(C24);
  this->addField(C25);
  this->addField(C26);
  this->addField(C33);
  this->addField(C34);
  this->addField(C35);
  this->addField(C36);
  this->addField(C44);
  this->addField(C45);
  this->addField(C46);
  this->addField(C55);
  this->addField(C56);
  this->addField(C66);
  this->addField(materialDensity);

  C11 = obj.C11;
  C12 = obj.C12;
  C13 = obj.C13;
  C14 = obj.C14;
  C15 = obj.C15;
  C16 = obj.C16;
  C22 = obj.C22;
  C23 = obj.C23;
  C24 = obj.C24;
  C25 = obj.C25;
  C26 = obj.C26;
  C33 = obj.C33;
  C34 = obj.C34;
  C35 = obj.C35;
  C36 = obj.C36;
  C44 = obj.C44;
  C45 = obj.C45;
  C46 = obj.C46;
  C55 = obj.C55;
  C56 = obj.C56;
  C66 = obj.C66;
  materialDensity = obj.materialDensity;
}


void FFlPMAT3D::convertUnits(const FFaUnitCalculator* convCal)
{
  // Round to 10 significant digits
  convCal->convert(C11.data(),"FORCE/AREA",10);
  convCal->convert(C12.data(),"FORCE/AREA",10);
  convCal->convert(C13.data(),"FORCE/AREA",10);
  convCal->convert(C14.data(),"FORCE/AREA",10);
  convCal->convert(C15.data(),"FORCE/AREA",10);
  convCal->convert(C16.data(),"FORCE/AREA",10);
  convCal->convert(C22.data(),"FORCE/AREA",10);
  convCal->convert(C23.data(),"FORCE/AREA",10);
  convCal->convert(C24.data(),"FORCE/AREA",10);
  convCal->convert(C25.data(),"FORCE/AREA",10);
  convCal->convert(C26.data(),"FORCE/AREA",10);
  convCal->convert(C33.data(),"FORCE/AREA",10);
  convCal->convert(C34.data(),"FORCE/AREA",10);
  convCal->convert(C35.data(),"FORCE/AREA",10);
  convCal->convert(C36.data(),"FORCE/AREA",10);
  convCal->convert(C44.data(),"FORCE/AREA",10);
  convCal->convert(C45.data(),"FORCE/AREA",10);
  convCal->convert(C46.data(),"FORCE/AREA",10);
  convCal->convert(C55.data(),"FORCE/AREA",10);
  convCal->convert(C56.data(),"FORCE/AREA",10);
  convCal->convert(C66.data(),"FORCE/AREA",10);
  convCal->convert(materialDensity.data(),"MASS/VOLUME",10);
}


void FFlPMAT3D::init()
{
  FFlPMAT3DTypeInfoSpec::instance()->setTypeName("PMAT3D");
  FFlPMAT3DTypeInfoSpec::instance()->setDescription("Material properties");
  FFlPMAT3DTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::MATERIAL_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPMAT3DTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPMAT3D::create,int,FFlAttributeBase*&));
}
