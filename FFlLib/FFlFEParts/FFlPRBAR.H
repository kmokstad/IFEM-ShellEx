// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PRBAR_H
#define FFL_PRBAR_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlPRBAR : public FFlAttributeBase
{
public:
  FFlPRBAR(int ID);
  FFlPRBAR(const FFlPRBAR& obj);
  virtual ~FFlPRBAR() {}

  virtual bool isIdentic(const FFlAttributeBase* otherAttrib) const;

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPRBAR);
  FFL_FE_ATTRIBUTE_SPEC(FFlPRBAR);
  FFL_TYPE_INFO(FFlPRBAR);

  FFlField<int> CNA;
  FFlField<int> CNB;
  FFlField<int> CMA;
  FFlField<int> CMB;
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlPRBAR>    FFlPRBARTypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPRBAR> FFlPRBARAttributeSpec;

#endif
