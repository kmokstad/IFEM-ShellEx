// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PBEAM_PIN_H
#define FFL_PBEAM_PIN_H

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


class FFlPBEAMPIN : public FFlAttributeBase
{
public:
  FFlPBEAMPIN(int ID);
  FFlPBEAMPIN(const FFlPBEAMPIN& obj);
  virtual ~FFlPBEAMPIN() {}

  virtual bool isIdentic(const FFlAttributeBase* otherAttrib) const;

  static void init();

  FFL_FE_ATTRIBUTE_FACTORY_INIT(FFlPBEAMPIN);
  FFL_TYPE_INFO(FFlPBEAMPIN);

  FFlField<int> PA;
  FFlField<int> PB;
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
