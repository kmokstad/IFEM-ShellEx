// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PCOORD_SYS_H
#define FFL_PCOORD_SYS_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include "FFaLib/FFaPatterns/FFaSingelton.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


class FFlPCOORDSYS : public FFlAttributeBase
{		
public:
  FFlPCOORDSYS(int ID);
  FFlPCOORDSYS(const FFlPCOORDSYS& obj);
  virtual ~FFlPCOORDSYS() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPCOORDSYS);
  FFL_FE_ATTRIBUTE_SPEC(FFlPCOORDSYS);
  FFL_TYPE_INFO(FFlPCOORDSYS);

  FFlField<FaVec3> Origo;
  FFlField<FaVec3> Zaxis;
  FFlField<FaVec3> XZpnt;
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
