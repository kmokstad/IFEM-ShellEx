// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PBEAM_ECCENT_H
#define FFL_PBEAM_ECCENT_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlPBEAMECCENT : public FFlAttributeBase
{
public:
  FFlPBEAMECCENT(int ID);
  FFlPBEAMECCENT(const FFlPBEAMECCENT& obj);
  virtual ~FFlPBEAMECCENT() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  virtual void resize(size_t n);

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPBEAMECCENT);
  FFL_FE_ATTRIBUTE_SPEC(FFlPBEAMECCENT);
  FFL_TYPE_INFO(FFlPBEAMECCENT);

  FFlField<FaVec3> node1Offset;
  FFlField<FaVec3> node2Offset;
  FFlField<FaVec3> node3Offset;
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlPBEAMECCENT>    FFlPBEAMECCENTTypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPBEAMECCENT> FFlPBEAMECCENTAttributeSpec;

#endif
