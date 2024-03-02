// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_TRI3_H
#define FFL_TRI3_H

#include "FFlShellElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif
#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlTRI3 : public FFlShellElementBase
{
public:
  FFlTRI3(int ID) : FFlShellElementBase(ID) {}
  virtual ~FFlTRI3() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlTRI3);
  FFL_FE_ELEMENT_TOP_SPEC(FFlTRI3);
  FFL_FE_ATTRIBUTE_SPEC(FFlTRI3);
  FFL_TYPE_INFO(FFlTRI3);

  virtual FaMat33 getGlobalizedElmCS() const;
  virtual bool getFaceNormals(std::vector<FaVec3>& normals, short int face = 1,
			      bool switchNormal = false) const;

  virtual bool getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const;

  virtual FaVec3 mapping(double xi, double eta, double) const;
  virtual bool invertMapping(const FaVec3& X, double* Xi) const;
  virtual FaVec3 interpolate(const double* Xi, const std::vector<FaVec3>& v) const;

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


typedef FFaSingelton<FFlFEAttributeSpec,FFlTRI3>  FFlTRI3AttributeSpec;
typedef FFaSingelton<FFlFEElementTopSpec,FFlTRI3> FFlTRI3ElementTopSpec;
typedef FFaSingelton<FFlTypeInfoSpec,FFlTRI3>     FFlTRI3TypeInfoSpec;

#endif
