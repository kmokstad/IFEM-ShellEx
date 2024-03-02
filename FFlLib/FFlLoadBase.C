// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlLoadBase.H"
#include "FFlLib/FFlFieldBase.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFaLib/FFaAlgebra/FFaCheckSum.H"


const std::string& FFlLoadBase::getTypeName() const
{
  return this->getTypeInfoSpec()->getTypeName();
}


void FFlLoadBase::calculateChecksum(FFaCheckSum* cs, int csMask) const
{
  this->FFlPartBase::checksum(cs,csMask);
  this->FFlFEAttributeRefs::checksum(cs);
  for (FFlFieldBase* field : myFields)
    field->calculateChecksum(cs);
}
