// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_TRI6_H
#define FFL_TRI6_H

#include "FFlLib/FFlFEParts/FFlShellElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif
#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlTRI6 : public FFlShellElementBase
{
public:
  FFlTRI6(int ID) : FFlShellElementBase(ID) {}
  virtual ~FFlTRI6() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlTRI6);
  FFL_FE_ELEMENT_TOP_SPEC(FFlTRI6);
  FFL_FE_ATTRIBUTE_SPEC(FFlTRI6);
  FFL_TYPE_INFO(FFlTRI6);

  //! \brief Splits the 6-noded shell element into four 3-noded elements.
  virtual bool split(Elements& newElm, FFlLinkHandler* owner, int);

  virtual FaMat33 getGlobalizedElmCS() const;
  virtual bool getFaceNormals(std::vector<FaVec3>& normals, short int face = 1,
			      bool switchNormal = false) const;

  virtual bool getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const;
  virtual int getNodalCoor(double* X, double* Y, double* Z) const;

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


typedef FFaSingelton<FFlFEAttributeSpec,FFlTRI6>  FFlTRI6AttributeSpec;
typedef FFaSingelton<FFlFEElementTopSpec,FFlTRI6> FFlTRI6ElementTopSpec;
typedef FFaSingelton<FFlTypeInfoSpec,FFlTRI6>     FFlTRI6TypeInfoSpec;

#endif
