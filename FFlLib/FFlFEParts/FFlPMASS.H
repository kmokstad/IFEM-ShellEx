// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PMASS_H
#define FFL_PMASS_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


class FFlPMASS : public FFlAttributeBase
{
  using DoubleVec = std::vector<double>;

public:
  FFlPMASS(int ID);
  FFlPMASS(const FFlPMASS& obj);
  virtual ~FFlPMASS() {}

  virtual void calculateChecksum(FFaCheckSum* cs, int csMask) const;
  virtual void convertUnits(const FFaUnitCalculator* convCal);
  virtual bool isIdentic(const FFlAttributeBase* otherAttrib) const;

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPMASS);
  FFL_TYPE_INFO(FFlPMASS);

  FFlField<DoubleVec> M;
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
