// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPMAT3D.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"


FFlPMAT3D::FFlPMAT3D(int id) : FFlAttributeBase(id)
{
  for (FFlField<double>& field : C)
    this->addField(field);

  this->addField(materialDensity);

  // default values - Typical for steel
  double E  = 205.0e9;
  double nu = 0.29;
  double f1 = E * (1.0-nu) / ((1.0+nu)*(1.0-2.0*nu));
  double f2 = f1 * nu / (1.0-nu);
  double f3 = f1 * (0.5-nu);

  C[ 0].setValue(f1); // C(1,1)
  C[ 1].setValue(f2); // C(1,2)
  C[ 2].setValue(f2); // C(1,3)
  C[ 6].setValue(f1); // C(2,2)
  C[ 7].setValue(f2); // C(2,3)
  C[11].setValue(f1); // C(3,3)
  C[15].setValue(f3); // C(4,4)
  C[18].setValue(f3); // C(5,5)
  C[20].setValue(f3); // C(6,6)

  materialDensity.setValue(7850.0);
}


FFlPMAT3D::FFlPMAT3D(const FFlPMAT3D& obj) : FFlAttributeBase(obj)
{
  for (size_t i = 0; i < C.size(); i++)
  {
    this->addField(C[i]);
    C[i].setValue(obj.C[i].getValue());
  }

  this->addField(materialDensity);
  materialDensity.setValue(obj.materialDensity.getValue());
}


void FFlPMAT3D::convertUnits(const FFaUnitCalculator* convCal)
{
  // Round to 10 significant digits
  for (FFlField<double>& field : C)
    convCal->convert(field.data(),"FORCE/AREA",10);
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
