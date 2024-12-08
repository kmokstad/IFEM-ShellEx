// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PBUSH_COEFF_H
#define FFL_PBUSH_COEFF_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"

#include <array>

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


class FFlPBUSHCOEFF : public FFlAttributeBase
{
public:
  FFlPBUSHCOEFF(int ID);
  FFlPBUSHCOEFF(const FFlPBUSHCOEFF& obj);
  virtual ~FFlPBUSHCOEFF() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPBUSHCOEFF);
  FFL_FE_ATTRIBUTE_SPEC(FFlPBUSHCOEFF);
  FFL_TYPE_INFO(FFlPBUSHCOEFF);

  std::array<FFlField<double>,6> K; //!< Diagonal 6x6 matrix
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
