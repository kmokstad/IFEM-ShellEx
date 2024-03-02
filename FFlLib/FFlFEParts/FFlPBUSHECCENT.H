// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PBUSH_ECCENT_H
#define FFL_PBUSH_ECCENT_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlPBUSHECCENT : public FFlAttributeBase
{
public:
  FFlPBUSHECCENT(int ID);
  FFlPBUSHECCENT(const FFlPBUSHECCENT& obj);
  virtual ~FFlPBUSHECCENT() {}

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPBUSHECCENT);
  FFL_FE_ATTRIBUTE_SPEC(FFlPBUSHECCENT);
  FFL_TYPE_INFO(FFlPBUSHECCENT);

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  FFlField<FaVec3> offset;
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlPBUSHECCENT>    FFlPBUSHECCENTTypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPBUSHECCENT> FFlPBUSHECCENTAttributeSpec;

#endif
