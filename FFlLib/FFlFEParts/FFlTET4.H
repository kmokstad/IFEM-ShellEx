// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_TET4_H
#define FFL_TET4_H

#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif
#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlTET4 : public FFlElementBase
{
public:
  FFlTET4(int ID) : FFlElementBase(ID) {}
  virtual ~FFlTET4() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlTET4);
  FFL_FE_ELEMENT_TOP_SPEC(FFlTET4);
  FFL_FE_ATTRIBUTE_SPEC(FFlTET4);
  FFL_TYPE_INFO(FFlTET4);

  virtual bool getFaceNormals(std::vector<FaVec3>& normals, short int face = 1,
			      bool switchNormal = false) const;

  virtual bool getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const;

  virtual int checkOrientation(bool fixIt = false);

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlTET4>     FFlTET4TypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlTET4>  FFlTET4AttributeSpec;
typedef FFaSingelton<FFlFEElementTopSpec,FFlTET4> FFlTET4ElementTopSpec;

#endif
