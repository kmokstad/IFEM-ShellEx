// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_SHELLS_H
#define FFL_SHELLS_H

#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


class FFlShellElementBase : public FFlElementBase
{
public:
  FFlShellElementBase(int ID) : FFlElementBase(ID) {}
  virtual ~FFlShellElementBase() {}

  double getThickness() const;

  virtual double getMassDensity() const;

  static double offPlaneTol; //!< Relative out-of-plane point search tolerance
};


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

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
