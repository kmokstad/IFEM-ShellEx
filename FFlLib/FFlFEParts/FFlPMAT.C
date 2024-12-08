// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPMAT.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


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
  using TypeInfoSpec = FFaSingelton<FFlTypeInfoSpec,FFlPMAT>;

  TypeInfoSpec::instance()->setTypeName("PMAT");
  TypeInfoSpec::instance()->setDescription("Material properties");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::MATERIAL_PROP);

  AttributeFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPMAT::create,int,FFlAttributeBase*&));
}


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
  using TypeInfoSpec = FFaSingelton<FFlTypeInfoSpec,FFlPMAT2D>;

  TypeInfoSpec::instance()->setTypeName("PMAT2D");
  TypeInfoSpec::instance()->setDescription("Material properties");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::MATERIAL_PROP);

  AttributeFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPMAT2D::create,int,FFlAttributeBase*&));
}


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
  using TypeInfoSpec = FFaSingelton<FFlTypeInfoSpec,FFlPMAT3D>;

  TypeInfoSpec::instance()->setTypeName("PMAT3D");
  TypeInfoSpec::instance()->setDescription("Material properties");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::MATERIAL_PROP);

  AttributeFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPMAT3D::create,int,FFlAttributeBase*&));
}


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
  materialDensity = 7.85E+3;
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
  using TypeInfoSpec = FFaSingelton<FFlTypeInfoSpec,FFlPMATSHELL>;

  TypeInfoSpec::instance()->setTypeName("PMATSHELL");
  TypeInfoSpec::instance()->setDescription("Material properties");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::MATERIAL_PROP);

  AttributeFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPMATSHELL::create,int,FFlAttributeBase*&));
}

#ifdef FF_NAMESPACE
} // namespace
#endif
