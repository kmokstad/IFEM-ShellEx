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

  FFlField<double> C11;
  FFlField<double> C12;
  FFlField<double> C13;
  FFlField<double> C14;
  FFlField<double> C15;
  FFlField<double> C16;
  FFlField<double> C22;
  FFlField<double> C23;
  FFlField<double> C24;
  FFlField<double> C25;
  FFlField<double> C26;
  FFlField<double> C33;
  FFlField<double> C34;
  FFlField<double> C35;
  FFlField<double> C36;
  FFlField<double> C44;
  FFlField<double> C45;
  FFlField<double> C46;
  FFlField<double> C55;
  FFlField<double> C56;
  FFlField<double> C66;
  FFlField<double> materialDensity;
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlPMAT3D>    FFlPMAT3DTypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPMAT3D> FFlPMAT3DAttributeSpec;

#endif
