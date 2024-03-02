// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PFATIGUE_H
#define FFL_PFATIGUE_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlPFATIGUE : public FFlAttributeBase
{
public:
  FFlPFATIGUE(int ID);
  FFlPFATIGUE(const FFlPFATIGUE& obj);
  virtual ~FFlPFATIGUE() {}

  virtual bool isIdentic(const FFlAttributeBase* otherAttrib) const;

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPFATIGUE);
  FFL_FE_ATTRIBUTE_SPEC(FFlPFATIGUE);
  FFL_TYPE_INFO(FFlPFATIGUE);

  FFlField<int>    snCurveStd;
  FFlField<int>    snCurveIndex;
  FFlField<double> stressConcentrationFactor;
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlPFATIGUE>    FFlPFATIGUETypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPFATIGUE> FFlPFATIGUEAttributeSpec;

#endif
