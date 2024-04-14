// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PBEAM_SECTION_H
#define FFL_PBEAM_SECTION_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlPBEAMSECTION : public FFlAttributeBase
{
public:
  FFlPBEAMSECTION(int ID);
  FFlPBEAMSECTION(const FFlPBEAMSECTION& obj);
  virtual ~FFlPBEAMSECTION() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPBEAMSECTION);
  FFL_FE_ATTRIBUTE_SPEC(FFlPBEAMSECTION);
  FFL_TYPE_INFO(FFlPBEAMSECTION);

public:
  FFlField<double> crossSectionArea;
  FFlField<double> Iy;
  FFlField<double> Iz;
  FFlField<double> It;
  FFlField<double> Kxy;
  FFlField<double> Kxz;
  FFlField<double> Sy;
  FFlField<double> Sz;
  FFlField<double> phi;

private:
  void setupFields();
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlPBEAMSECTION>    FFlPBEAMSECTIONTypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPBEAMSECTION> FFlPBEAMSECTIONAttributeSpec;

#endif
