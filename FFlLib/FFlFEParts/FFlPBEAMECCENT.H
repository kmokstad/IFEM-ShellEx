// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PBEAM_ECCENT_H
#define FFL_PBEAM_ECCENT_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


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
  FFL_TYPE_INFO(FFlPBEAMECCENT);

  FFlField<FaVec3> node1Offset;
  FFlField<FaVec3> node2Offset;
  FFlField<FaVec3> node3Offset;
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
