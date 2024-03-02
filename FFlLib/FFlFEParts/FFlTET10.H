// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_TET10_H
#define FFL_TET10_H

#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif
#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlTET10 : public FFlElementBase
{
public:
  FFlTET10(int ID) : FFlElementBase(ID) {}
  virtual ~FFlTET10() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlTET10);
  FFL_FE_ELEMENT_TOP_SPEC(FFlTET10);
  FFL_FE_ATTRIBUTE_SPEC(FFlTET10);
  FFL_TYPE_INFO(FFlTET10);

  virtual bool getFaceNormals(std::vector<FaVec3>& normals, short int face = 1,
			      bool switchNormal = false) const;

  virtual bool getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const;

  virtual int checkOrientation(bool fixIt = false);

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlTET10>     FFlTET10TypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlTET10>  FFlTET10AttributeSpec;
typedef FFaSingelton<FFlFEElementTopSpec,FFlTET10> FFlTET10ElementTopSpec;

#endif
