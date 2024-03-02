// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PMAT_2D_H
#define FFL_PMAT_2D_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"


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

  FFlField<double> C11;
  FFlField<double> C12;
  FFlField<double> C13;
  FFlField<double> C22;
  FFlField<double> C23;
  FFlField<double> C33;
  FFlField<double> materialDensity;
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlPMAT2D>    FFlPMAT2DTypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPMAT2D> FFlPMAT2DAttributeSpec;

#endif
