// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PCOMP_H
#define FFL_PCOMP_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlPCOMP : public FFlAttributeBase
{
public:
  class Ply
  {
  public:
    int MID;      // material id
    double T;     // thickness
    double thetaInDeg; // Orientation angle ply longitudinal direction
    Ply() : MID(1), T(20.0), thetaInDeg(30.0) {}
    friend bool operator==(const Ply& a, const Ply& b) { return &a == &b; }
  };
  typedef std::vector<Ply> PlyVec;

  friend std::ostream& operator<< (std::ostream& os, const PlyVec& val);

  FFlPCOMP(int ID);
  FFlPCOMP(const FFlPCOMP& obj);
  virtual ~FFlPCOMP() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPCOMP);
  FFL_FE_ATTRIBUTE_SPEC(FFlPCOMP);
  FFL_TYPE_INFO(FFlPCOMP);

  virtual void calculateChecksum(FFaCheckSum* cs, int csMask) const;

  FFlField<double> Z0;
  FFlField<PlyVec> plySet;
};


template<> void FFaCheckSum::add(const FFlPCOMP::PlyVec& v);


typedef FFaSingelton<FFlTypeInfoSpec,FFlPCOMP>    FFlPCOMPTypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlPCOMP> FFlPCOMPAttributeSpec;

#endif
