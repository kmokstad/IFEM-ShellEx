// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PORIENT_H
#define FFL_PORIENT_H

#include <array>

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


class FFlPORIENT : public FFlAttributeBase
{		
public:
  FFlPORIENT(int ID);
  FFlPORIENT(const FFlPORIENT& obj);
  virtual ~FFlPORIENT() {}

  virtual bool isIdentic(const FFlAttributeBase* otherAttrib) const;

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPORIENT);
  FFL_TYPE_INFO(FFlPORIENT);

  FFlField<FaVec3> directionVector;
};


/*!
  \brief Class for nodal orientation vectors for 3-noded beam elements.
*/

class FFlPORIENT3 : public FFlAttributeBase
{
public:
  FFlPORIENT3(int ID);
  FFlPORIENT3(const FFlPORIENT3& obj);
  virtual ~FFlPORIENT3() {}

  virtual bool isIdentic(const FFlAttributeBase* otherAttrib) const;

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPORIENT3);
  FFL_TYPE_INFO(FFlPORIENT3);

  std::array<FFlField<FaVec3>,3> directionVector;
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
