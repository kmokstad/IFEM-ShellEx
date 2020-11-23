// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPMAT2D.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"


FFlPMAT2D::FFlPMAT2D(int id) : FFlAttributeBase(id)
{
  for (FFlField<double>& field : C)
    this->addField(field);

  this->addField(materialDensity);

  // default values - Typical for steel
  double E  = 205.0e9;
  double nu = 0.29;
  double f0 = E / (1.0-nu*nu);

  C[0].setValue(f0);            // C(1,1)
  C[1].setValue(f0*nu);         // C(1,2)
  C[3].setValue(f0);            // C(2,2)
  C[5].setValue(f0*(1-nu)*0.5); // C(3,3)

  materialDensity.setValue(7850.0);
}


FFlPMAT2D::FFlPMAT2D(const FFlPMAT2D& obj) : FFlAttributeBase(obj)
{
  for (size_t i = 0; i < C.size(); i++)
  {
    this->addField(C[i]);
    C[i].setValue(obj.C[i].getValue());
  }

  this->addField(materialDensity);
  materialDensity.setValue(obj.materialDensity.getValue());
}


void FFlPMAT2D::convertUnits(const FFaUnitCalculator* convCal)
{
  // Round to 10 significant digits
  for (FFlField<double>& field : C)
    convCal->convert(field.data(),"FORCE/AREA",10);
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
