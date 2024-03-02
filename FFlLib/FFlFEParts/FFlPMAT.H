// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PMAT_H
#define FFL_PMAT_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"


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

typedef FFaSingelton<FFlTypeInfoSpec,FFlPMAT>    FFlPMATTypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPMAT> FFlPMATAttributeSpec;

#endif
