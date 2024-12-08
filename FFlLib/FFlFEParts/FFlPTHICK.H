// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PTHICK_H
#define FFL_PTHICK_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


class FFlPTHICK : public FFlAttributeBase
{		
public:
  FFlPTHICK(int ID);
  FFlPTHICK(const FFlPTHICK& obj);
  virtual ~FFlPTHICK() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPTHICK);
  FFL_FE_ATTRIBUTE_SPEC(FFlPTHICK);
  FFL_TYPE_INFO(FFlPTHICK);

  FFlField<double> thickness;
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
