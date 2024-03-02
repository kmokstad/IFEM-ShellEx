// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PHEIGHT_H
#define FFL_PHEIGHT_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlPHEIGHT : public FFlAttributeBase
{		
public:
  FFlPHEIGHT(int ID);
  FFlPHEIGHT(const FFlPHEIGHT& obj);
  virtual ~FFlPHEIGHT() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPHEIGHT);
  FFL_FE_ATTRIBUTE_SPEC(FFlPHEIGHT);
  FFL_TYPE_INFO(FFlPHEIGHT);

  FFlField<double> height;
};

typedef FFaSingelton<FFlTypeInfoSpec,FFlPHEIGHT>    FFlPHEIGHTTypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPHEIGHT> FFlPHEIGHTAttributeSpec;

#endif
