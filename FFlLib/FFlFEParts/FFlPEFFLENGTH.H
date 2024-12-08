// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PEFF_LENGTH_H_
#define FFL_PEFF_LENGTH_H_

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


class FFlPEFFLENGTH : public FFlAttributeBase
{		
public:
  FFlPEFFLENGTH(int ID);
  FFlPEFFLENGTH(const FFlPEFFLENGTH& obj);
  virtual ~FFlPEFFLENGTH() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);
  virtual bool isIdentic(const FFlAttributeBase* otherAttrib) const;

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPEFFLENGTH);
  FFL_FE_ATTRIBUTE_SPEC(FFlPEFFLENGTH);
  FFL_TYPE_INFO(FFlPEFFLENGTH);

  FFlField<double> length;
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
