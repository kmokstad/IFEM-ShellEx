// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_LINK_CS_MASK_H
#define FFL_LINK_CS_MASK_H


namespace FFl
{
  enum CSMask { CS_NOEXTINFO    = 0x000001,
                CS_EXTMASK      = 0x00000f,
                CS_NOGROUPINFO  = 0x000010,
                CS_GROUPMASK    = 0x0000f0,
                CS_NOSTRCINFO   = 0x000100,
                CS_STRCMASK     = 0x000f00,
                CS_NOVISUALINFO = 0x001000,
                CS_VISUALMASK   = 0x00f000,
                CS_NOIDINFO     = 0x010000,
                CS_IDMASK       = 0x0f0000 };
}

#endif
