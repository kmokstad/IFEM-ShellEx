// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PSPRING_H
#define FFL_PSPRING_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"

#include <array>

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


class FFlPSPRING : public FFlAttributeBase
{
public:
  FFlPSPRING(int ID);
  FFlPSPRING(const FFlPSPRING& obj);
  virtual ~FFlPSPRING() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPSPRING);
  FFL_FE_ATTRIBUTE_SPEC(FFlPSPRING);
  FFL_TYPE_INFO(FFlPSPRING);

  enum SpringType { TRANS_SPRING = 1, ROT_SPRING = 2 };

  std::array<FFlField<double>,21> K; //!< Upper triangle of 6x6 matrix

  FFlField<int> type; //!< Needed to determine correct unit conversion
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
