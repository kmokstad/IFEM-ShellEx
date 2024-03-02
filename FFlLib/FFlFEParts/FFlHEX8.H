// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_HEX8_H
#define FFL_HEX8_H

#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif
#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlHEX8 : public FFlElementBase
{
public:
  FFlHEX8(int ID) : FFlElementBase(ID) {}
  virtual ~FFlHEX8() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlHEX8);
  FFL_FE_ELEMENT_TOP_SPEC(FFlHEX8);
  FFL_FE_ATTRIBUTE_SPEC(FFlHEX8);
  FFL_TYPE_INFO(FFlHEX8);

  virtual bool getFaceNormals(std::vector<FaVec3>& normals, short int face = 1,
			      bool switchNormal = false) const;

  virtual bool getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const;

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlHEX8>     FFlHEX8TypeInfoSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlHEX8>  FFlHEX8AttributeSpec;
typedef FFaSingelton<FFlFEElementTopSpec,FFlHEX8> FFlHEX8ElementTopSpec;

#endif
