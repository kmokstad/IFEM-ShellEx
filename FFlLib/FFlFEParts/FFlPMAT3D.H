// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PMAT_3D_H
#define FFL_PMAT_3D_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"

#include <array>


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


typedef FFaSingelton<FFlTypeInfoSpec,FFlPMAT3D>    FFlPMAT3DTypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPMAT3D> FFlPMAT3DAttributeSpec;

#endif
