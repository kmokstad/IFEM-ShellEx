// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_QUAD8_H
#define FFL_QUAD8_H

#include "FFlLib/FFlFEParts/FFlShellElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif
#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlQUAD8 : public FFlShellElementBase
{
public:
  FFlQUAD8(int ID) : FFlShellElementBase(ID) {}
  virtual ~FFlQUAD8() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlQUAD8);
  FFL_FE_ELEMENT_TOP_SPEC(FFlQUAD8);
  FFL_FE_ATTRIBUTE_SPEC(FFlQUAD8);
  FFL_TYPE_INFO(FFlQUAD8);

  //! \brief Splits the 8-noded shell element into four 4-noded elements.
  virtual bool split(Elements& newElm, FFlLinkHandler* owner, int centerNode);

  virtual FaMat33 getGlobalizedElmCS() const;
  virtual bool getFaceNormals(std::vector<FaVec3>& normals, short int face = 1,
			      bool switchNormal = false) const;

  virtual bool getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const;

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


typedef FFaSingelton<FFlFEAttributeSpec,FFlQUAD8>  FFlQUAD8AttributeSpec;
typedef FFaSingelton<FFlFEElementTopSpec,FFlQUAD8> FFlQUAD8ElementTopSpec;
typedef FFaSingelton<FFlTypeInfoSpec,FFlQUAD8>     FFlQUAD8TypeInfoSpec;

#endif
