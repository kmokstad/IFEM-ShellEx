// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFieldBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFaLib/FFaAlgebra/FFaCheckSum.H"

#include <cstring>
#include <iostream>
#if FFL_DEBUG > 2
#include <iomanip>
#endif

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


const std::string& FFlAttributeBase::getTypeName() const
{
  return this->getTypeInfoSpec()->getTypeName();
}


const std::string& FFlAttributeBase::getDescription() const
{
  return this->getTypeInfoSpec()->getDescription();
}


void FFlAttributeBase::calculateChecksum(FFaCheckSum* cs, int csMask) const
{
  this->FFlNamedPartBase::checksum(cs,csMask);
  this->FFlFEAttributeRefs::checksum(cs);

  for (FFlFieldBase* field : myFields)
  {
    field->calculateChecksum(cs);
#if FFL_DEBUG > 2
    std::cout <<"Checksum after "<< this->getTypeName()
              <<" field "<< std::setprecision(16) << *field
              <<" : "<< cs->getCurrent() << std::endl;
#endif
  }
}


bool FFlAttributeBase::isIdentic(const FFlAttributeBase* other) const
{
#ifdef FFL_DEBUG
  std::cout <<"  ** FFlAttributeBase::isIdentic() is not implemented for "
            << this->getTypeName() <<", comparing pointers only."<< std::endl;
#endif
  return this == other;
}


void FFlAttributeBase::print(const char* prefix) const
{
  std::cout << prefix << this->getTypeName() <<", ID = "<< this->getID() <<",";
  if (!this->getName().empty())
  {
    std::cout <<" Name: "<< this->getName();
    size_t nindent = strlen(prefix) + this->getTypeName().size() + 8;
    for (int id = this->getID(); id > 0; id /= 10) nindent++;
    std::cout <<"\n" << std::string(nindent,' ');
  }
  std::cout <<" Fields:";
  for (FFlFieldBase* field : myFields)
    std::cout <<" "<< *field;
  std::cout << std::endl;
}

#ifdef FF_NAMESPACE
} // namespace
#endif
