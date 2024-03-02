// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlPartBase.H"
#include "FFlLib/FFlLinkCSMask.H"
#include "FFaLib/FFaAlgebra/FFaCheckSum.H"


void FFlPartBase::checksum(FFaCheckSum* cs, int csMask) const
{
  if ((csMask & FFl::CS_IDMASK) != FFl::CS_NOIDINFO)
    cs->add(myID);
}
