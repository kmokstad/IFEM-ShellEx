// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_BEAM3_H
#define FFL_BEAM3_H

#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif
#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlBEAM3 : public FFlElementBase
{
public:
  FFlBEAM3(int ID) : FFlElementBase(ID) {}
  virtual ~FFlBEAM3() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlBEAM3);
  FFL_FE_ELEMENT_TOP_SPEC(FFlBEAM3);
  FFL_FE_ATTRIBUTE_SPEC(FFlBEAM3);
  FFL_TYPE_INFO(FFlBEAM3);

  virtual double getMassDensity() const;
  virtual bool getVolumeAndInertia(double& volume, FaVec3& cog,
                                   FFaTensor3& inertia) const;

  //! \brief Splits the 3-noded beam element into two 2-noded elements.
  virtual bool split(Elements& newElm, FFlLinkHandler* owner, int);

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


typedef FFaSingelton<FFlFEAttributeSpec,FFlBEAM3>  FFlBEAM3AttributeSpec;
typedef FFaSingelton<FFlFEElementTopSpec,FFlBEAM3> FFlBEAM3ElementTopSpec;
typedef FFaSingelton<FFlTypeInfoSpec,FFlBEAM3>     FFlBEAM3TypeInfoSpec;

#endif
