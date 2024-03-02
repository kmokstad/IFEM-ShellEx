// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_MEM_POOL_H
#define FFL_MEM_POOL_H

class FFlLinkHandler;


namespace FFlMemPool
{
  void deleteAllLinkMemPools();
  void deleteVisualsMemPools();
  void setAsMemPoolPart(FFlLinkHandler*);
  void freeMemPoolPart(FFlLinkHandler*);
  void resetMemPoolPart();
}

#endif
