// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_SOLIDS_H
#define FFL_SOLIDS_H

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


class FFlWEDG15 : public FFlElementBase
{
public:
  FFlWEDG15(int ID) : FFlElementBase(ID) {}
  virtual ~FFlWEDG15() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlWEDG15);
  FFL_FE_ELEMENT_TOP_SPEC(FFlWEDG15);
  FFL_FE_ATTRIBUTE_SPEC(FFlWEDG15);
  FFL_TYPE_INFO(FFlWEDG15);

  virtual bool getFaceNormals(std::vector<FaVec3>& normals, short int face = 1,
                              bool switchNormal = false) const;

  virtual bool getVolumeAndInertia(double& volume, FaVec3& cog,
                                   FFaTensor3& inertia) const;

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


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

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
