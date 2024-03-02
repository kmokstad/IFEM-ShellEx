// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PRGD_H
#define FFL_PRGD_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlPRGD : public FFlAttributeBase
{
public:
  FFlPRGD(int ID);
  FFlPRGD(const FFlPRGD& obj);
  virtual ~FFlPRGD() {}

  virtual bool isIdentic(const FFlAttributeBase* otherAttrib) const;

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPRGD);
  FFL_FE_ATTRIBUTE_SPEC(FFlPRGD);
  FFL_TYPE_INFO(FFlPRGD);

  FFlField<int> dependentDofs;
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlPRGD>    FFlPRGDTypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPRGD> FFlPRGDAttributeSpec;

#endif
