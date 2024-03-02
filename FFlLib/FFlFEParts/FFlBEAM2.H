// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_BEAM2_H
#define FFL_BEAM2_H

#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif
#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlBEAM2 : public FFlElementBase
{
public:
  FFlBEAM2(int ID) : FFlElementBase(ID) {}
  virtual ~FFlBEAM2() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlBEAM2);
  FFL_FE_ELEMENT_TOP_SPEC(FFlBEAM2);
  FFL_FE_ATTRIBUTE_SPEC(FFlBEAM2);
  FFL_TYPE_INFO(FFlBEAM2);

  virtual double getMassDensity() const;
  virtual bool getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const;

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


typedef FFaSingelton<FFlFEAttributeSpec,FFlBEAM2>  FFlBEAM2AttributeSpec;
typedef FFaSingelton<FFlFEElementTopSpec,FFlBEAM2> FFlBEAM2ElementTopSpec;
typedef FFaSingelton<FFlTypeInfoSpec,FFlBEAM2>     FFlBEAM2TypeInfoSpec;

#endif
