// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_RBAR_H
#define FFL_RBAR_H

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


class FFlRBAR : public FFlElementBase
{
public:
  FFlRBAR(int ID) : FFlElementBase(ID) {}
  virtual ~FFlRBAR() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlRBAR);
  FFL_FE_ELEMENT_TOP_SPEC(FFlRBAR);
  FFL_FE_ATTRIBUTE_SPEC(FFlRBAR);
  FFL_TYPE_INFO(FFlRBAR);

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
