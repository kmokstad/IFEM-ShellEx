// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PMAT_SHELL_H
#define FFL_PMAT_SHELL_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"


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


typedef FFaSingelton<FFlTypeInfoSpec,FFlPMATSHELL>    FFlPMATSHELLTypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPMATSHELL> FFlPMATSHELLAttributeSpec;

#endif
