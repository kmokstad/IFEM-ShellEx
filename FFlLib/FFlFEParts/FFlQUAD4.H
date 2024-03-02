// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_QUAD4_H
#define FFL_QUAD4_H

#include "FFlLib/FFlFEParts/FFlShellElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif
#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlQUAD4 : public FFlShellElementBase
{
public:
  FFlQUAD4(int ID) : FFlShellElementBase(ID) {}
  virtual ~FFlQUAD4() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlQUAD4);
  FFL_FE_ELEMENT_TOP_SPEC(FFlQUAD4);
  FFL_FE_ATTRIBUTE_SPEC(FFlQUAD4);
  FFL_TYPE_INFO(FFlQUAD4);

  virtual FaMat33 getGlobalizedElmCS() const;
  virtual bool getFaceNormals(std::vector<FaVec3>& normals, short int face = 1,
			      bool switchNormal = false) const;

  virtual bool getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const;

  virtual FaVec3 mapping(double xi, double eta, double = 0.0) const;
  virtual bool invertMapping(const FaVec3& X, double* Xi) const;
  virtual FaVec3 interpolate(const double* Xi, const std::vector<FaVec3>& v) const;

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


typedef FFaSingelton<FFlFEAttributeSpec,FFlQUAD4>  FFlQUAD4AttributeSpec;
typedef FFaSingelton<FFlFEElementTopSpec,FFlQUAD4> FFlQUAD4ElementTopSpec;
typedef FFaSingelton<FFlTypeInfoSpec,FFlQUAD4>     FFlQUAD4TypeInfoSpec;

#endif
