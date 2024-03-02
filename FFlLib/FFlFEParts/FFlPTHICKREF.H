// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PTHICK_REF_H
#define FFL_PTHICK_REF_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlPTHICKREF : public FFlAttributeBase
{
public:
  FFlPTHICKREF(int ID);
  FFlPTHICKREF(const FFlPTHICKREF& obj);
  virtual ~FFlPTHICKREF() {}

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPTHICKREF);
  FFL_FE_ATTRIBUTE_SPEC(FFlPTHICKREF);
  FFL_TYPE_INFO(FFlPTHICKREF);

  FFlField<double> factor;
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlPTHICKREF>    FFlPTHICKREFTypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPTHICKREF> FFlPTHICKREFAttributeSpec;

#endif
