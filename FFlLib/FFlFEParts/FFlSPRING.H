// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_SPRING_H
#define FFL_SPRING_H

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


/* TODO,kmo: This is not very elegant. The only difference between the two
   classes FFlSPRING and FFlRSPRING is that FFlSPRING has three DOFs per node
   whereas FFlRSPRING has six. Find a better way to handle this. */

// For translational DOFs only
class FFlSPRING : public FFlElementBase
{
public:
  FFlSPRING(int ID) : FFlElementBase(ID) {}
  virtual ~FFlSPRING() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlSPRING);
  FFL_FE_ELEMENT_TOP_SPEC(FFlSPRING);
  FFL_FE_ATTRIBUTE_SPEC(FFlSPRING);
  FFL_TYPE_INFO(FFlSPRING);

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


// For rotational DOFs
class FFlRSPRING : public FFlElementBase
{
public:
  FFlRSPRING(int ID) : FFlElementBase(ID) {}
  virtual ~FFlRSPRING() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlRSPRING);
  FFL_FE_ELEMENT_TOP_SPEC(FFlRSPRING);
  FFL_FE_ATTRIBUTE_SPEC(FFlRSPRING);
  FFL_TYPE_INFO(FFlRSPRING);

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
