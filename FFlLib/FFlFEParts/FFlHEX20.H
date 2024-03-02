// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_HEX20_H
#define FFL_HEX20_H

#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif
#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlHEX20 : public FFlElementBase
{
public:
  FFlHEX20(int ID) : FFlElementBase(ID) {}
  virtual ~FFlHEX20() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlHEX20);
  FFL_FE_ELEMENT_TOP_SPEC(FFlHEX20);
  FFL_FE_ATTRIBUTE_SPEC(FFlHEX20);
  FFL_TYPE_INFO(FFlHEX20);

  virtual bool getFaceNormals(std::vector<FaVec3>& normals, short int face = 1,
			      bool switchNormal = false) const;

  virtual bool getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const;

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlHEX20>     FFlHEX20TypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlHEX20>  FFlHEX20AttributeSpec;
typedef FFaSingelton<FFlFEElementTopSpec,FFlHEX20> FFlHEX20ElementTopSpec;

#endif
