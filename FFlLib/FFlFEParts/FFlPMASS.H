// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PMASS_H
#define FFL_PMASS_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlPMASS : public FFlAttributeBase
{
public:
  FFlPMASS(int ID);
  FFlPMASS(const FFlPMASS& obj);
  virtual ~FFlPMASS() {}

  virtual void calculateChecksum(FFaCheckSum* cs, int csMask) const;
  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPMASS);
  FFL_FE_ATTRIBUTE_SPEC(FFlPMASS);
  FFL_TYPE_INFO(FFlPMASS);

  FFlField< std::vector<double> > M;
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlPMASS>    FFlPMASSTypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPMASS> FFlPMASSAttributeSpec;

#endif
