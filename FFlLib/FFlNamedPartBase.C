// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlNamedPartBase.H"
#include "FFlLib/FFlLinkCSMask.H"
#include "FFaLib/FFaAlgebra/FFaCheckSum.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


void FFlNamedPartBase::setName(const std::string& groupName)
{
  userName = groupName;

  // Erase all instances of the "-character in the userName (if any).
  // The FTL-file reader does not cope with their presense (TT #2910).
  size_t i, j = 0;
  while (j < userName.size() && (i = userName.find('"',j)) < userName.size())
  {
    userName.erase(i,1);
    j = i;
  }
}


void FFlNamedPartBase::checksum(FFaCheckSum* cs, int csMask) const
{
  FFlPartBase::checksum(cs,csMask);

  if (!userName.empty() && (csMask & FFl::CS_GROUPMASK) != FFl::CS_NOGROUPINFO)
    cs->add(userName);
}

#ifdef FF_NAMESPACE
} // namespace
#endif
