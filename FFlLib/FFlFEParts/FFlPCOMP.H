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


//! \brief A struct representing one ply in a composite material layer.
struct FFlPly
{
  int    MID   = 1;    //!< Material Id
  double T     = 20.0; //!< Ply thickness
  double theta = 30.0; //!< Ply orientation angle [deg] longitudinal direction
};

using FFlPlyVec = std::vector<FFlPly>; //!< Convenience definition

std::ostream& operator<<(std::ostream& os, const FFlPlyVec& val);

template<> void FFaCheckSum::add(const FFlPlyVec& val);


#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif

class FFlPCOMP : public FFlAttributeBase
{
public:
  FFlPCOMP(int ID);
  FFlPCOMP(const FFlPCOMP& obj);
  virtual ~FFlPCOMP() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPCOMP);
  FFL_FE_ATTRIBUTE_SPEC(FFlPCOMP);
  FFL_TYPE_INFO(FFlPCOMP);

  virtual void calculateChecksum(FFaCheckSum* cs, int csMask) const;

  FFlField<double>    Z0;
  FFlField<FFlPlyVec> plySet;
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
