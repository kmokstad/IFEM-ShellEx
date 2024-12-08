// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PMAT_H
#define FFL_PMAT_H

#include <array>

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


class FFlPMAT : public FFlAttributeBase
{
public:
  FFlPMAT(int ID);
  FFlPMAT(const FFlPMAT& obj);
  virtual ~FFlPMAT() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPMAT);
  FFL_FE_ATTRIBUTE_SPEC(FFlPMAT);
  FFL_TYPE_INFO(FFlPMAT);

  FFlField<double> youngsModule;
  FFlField<double> shearModule;
  FFlField<double> poissonsRatio;
  FFlField<double> materialDensity;
};


class FFlPMAT2D : public FFlAttributeBase
{
public:
  FFlPMAT2D(int ID);
  FFlPMAT2D(const FFlPMAT2D& obj);
  virtual ~FFlPMAT2D() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPMAT2D);
  FFL_FE_ATTRIBUTE_SPEC(FFlPMAT2D);
  FFL_TYPE_INFO(FFlPMAT2D);

  std::array<FFlField<double>,6> C; //!< Upper triangle of 3x3 matrix

  FFlField<double> materialDensity;
};


class FFlPMAT3D : public FFlAttributeBase
{
public:
  FFlPMAT3D(int ID);
  FFlPMAT3D(const FFlPMAT3D& obj);
  virtual ~FFlPMAT3D() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPMAT3D);
  FFL_FE_ATTRIBUTE_SPEC(FFlPMAT3D);
  FFL_TYPE_INFO(FFlPMAT3D);

  std::array<FFlField<double>,21> C; //!< Upper triangle of 6x6 matrix

  FFlField<double> materialDensity;
};


class FFlPMATSHELL : public FFlAttributeBase
{
public:
  FFlPMATSHELL(int ID);
  FFlPMATSHELL(const FFlPMATSHELL& obj);
  virtual ~FFlPMATSHELL() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPMATSHELL);
  FFL_FE_ATTRIBUTE_SPEC(FFlPMATSHELL);
  FFL_TYPE_INFO(FFlPMATSHELL);

  FFlField<double> E1;
  FFlField<double> E2;
  FFlField<double> NU12;
  FFlField<double> G12;
  FFlField<double> G1Z;
  FFlField<double> G2Z;
  FFlField<double> materialDensity;
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
