// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PNSM_H
#define FFL_PNSM_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlPNSM : public FFlAttributeBase
{
public:
  FFlPNSM(int ID);
  FFlPNSM(const FFlPNSM& obj);
  virtual ~FFlPNSM() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPNSM);
  FFL_FE_ATTRIBUTE_SPEC(FFlPNSM);
  FFL_TYPE_INFO(FFlPNSM);

  FFlField<double> NSM;
  FFlField<int> isShell; // used to determine the proper unit conversion
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlPNSM>    FFlPNSMTypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPNSM> FFlPNSMAttributeSpec;

#endif
