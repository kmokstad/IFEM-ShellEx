// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_CMASS_H
#define FFL_CMASS_H

#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


class FFlCMASS : public FFlElementBase
{
public:
  FFlCMASS(int ID) : FFlElementBase(ID) {}
  virtual ~FFlCMASS() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlCMASS);
  FFL_FE_ELEMENT_TOP_SPEC(FFlCMASS);
  FFL_FE_ATTRIBUTE_SPEC(FFlCMASS);
  FFL_TYPE_INFO(FFlCMASS);

  virtual double getMassDensity() const { return 1.0; }
  virtual bool getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const;

private:
  static bool RtMR(const std::vector<double>& M, double R[6][6], double II[6][6]);

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
