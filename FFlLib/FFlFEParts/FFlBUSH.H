// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_BUSH_H
#define FFL_BUSH_H

#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif
#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlBUSH : public FFlElementBase
{
public:
  FFlBUSH(int ID) : FFlElementBase(ID) {}
  virtual ~FFlBUSH() {}

  static void init();

  FFL_FE_ELEMENT_FACTORY_INIT(FFlBUSH);
  FFL_FE_ELEMENT_TOP_SPEC(FFlBUSH);
  FFL_FE_ATTRIBUTE_SPEC(FFlBUSH);
  FFL_TYPE_INFO(FFlBUSH);

  virtual int getNodalCoor(double* X, double* Y, double* Z) const;

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


typedef FFaSingelton<FFlFEAttributeSpec,FFlBUSH>  FFlBUSHAttributeSpec;
typedef FFaSingelton<FFlFEElementTopSpec,FFlBUSH> FFlBUSHElementTopSpec;
typedef FFaSingelton<FFlTypeInfoSpec,FFlBUSH>     FFlBUSHTypeInfoSpec;

#endif
