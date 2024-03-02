// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_WEDG6_H
#define FFL_WEDG6_H

#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif
#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlWEDG6 : public FFlElementBase
{
public:
  FFlWEDG6(int ID) : FFlElementBase(ID) {}
  virtual ~FFlWEDG6() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlWEDG6);
  FFL_FE_ELEMENT_TOP_SPEC(FFlWEDG6);
  FFL_FE_ATTRIBUTE_SPEC(FFlWEDG6);
  FFL_TYPE_INFO(FFlWEDG6);

  virtual bool getFaceNormals(std::vector<FaVec3>& normals, short int face = 1,
			      bool switchNormal = false) const;

  virtual bool getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const;

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlWEDG6>     FFlWEDG6TypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlWEDG6>  FFlWEDG6AttributeSpec;
typedef FFaSingelton<FFlFEElementTopSpec,FFlWEDG6> FFlWEDG6ElementTopSpec;

#endif
