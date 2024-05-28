// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlMemPool.H"
#ifdef FT_USE_MEMPOOL
#ifdef FT_USE_VERTEX
#include "FFlLib/FFlVertex.H"
#endif
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlFEParts/FFlCMASS.H"
#include "FFlLib/FFlFEParts/FFlSPRING.H"
#include "FFlLib/FFlFEParts/FFlBUSH.H"
#include "FFlLib/FFlFEParts/FFlBEAM2.H"
#include "FFlLib/FFlFEParts/FFlBEAM3.H"
#include "FFlLib/FFlFEParts/FFlTRI3.H"
#include "FFlLib/FFlFEParts/FFlTRI6.H"
#include "FFlLib/FFlFEParts/FFlQUAD4.H"
#include "FFlLib/FFlFEParts/FFlQUAD8.H"
#include "FFlLib/FFlFEParts/FFlTET4.H"
#include "FFlLib/FFlFEParts/FFlTET10.H"
#include "FFlLib/FFlFEParts/FFlWEDG6.H"
#include "FFlLib/FFlFEParts/FFlWEDG15.H"
#include "FFlLib/FFlFEParts/FFlHEX8.H"
#include "FFlLib/FFlFEParts/FFlHEX20.H"
#include "FFlLib/FFlFEParts/FFlRBAR.H"
#include "FFlLib/FFlFEParts/FFlRGD.H"
#include "FFlLib/FFlFEParts/FFlWAVGM.H"
#ifdef FT_USE_STRAINCOAT
#include "FFlLib/FFlFEParts/FFlSTRCoat.H"
#endif

typedef FFaSingelton<FFaMemPoolMgr,FFlElementBase> FFlMemPoolMgr;

FFaMemPool FFlNode   ::ourMemPool(sizeof(FFlNode));
#ifdef FT_USE_VERTEX
FFaMemPool FFlVertex ::ourMemPool(sizeof(FFlVertex));
#endif
FFaMemPool FFlCMASS  ::ourMemPool(sizeof(FFlCMASS)  ,FFlMemPoolMgr::instance());
FFaMemPool FFlSPRING ::ourMemPool(sizeof(FFlSPRING) ,FFlMemPoolMgr::instance());
FFaMemPool FFlRSPRING::ourMemPool(sizeof(FFlRSPRING),FFlMemPoolMgr::instance());
FFaMemPool FFlBUSH   ::ourMemPool(sizeof(FFlBUSH)   ,FFlMemPoolMgr::instance());
FFaMemPool FFlBEAM2  ::ourMemPool(sizeof(FFlBEAM2)  ,FFlMemPoolMgr::instance());
FFaMemPool FFlBEAM3  ::ourMemPool(sizeof(FFlBEAM3)  ,FFlMemPoolMgr::instance());
FFaMemPool FFlTRI3   ::ourMemPool(sizeof(FFlTRI3)   ,FFlMemPoolMgr::instance());
FFaMemPool FFlTRI6   ::ourMemPool(sizeof(FFlTRI6)   ,FFlMemPoolMgr::instance());
FFaMemPool FFlQUAD4  ::ourMemPool(sizeof(FFlQUAD4)  ,FFlMemPoolMgr::instance());
FFaMemPool FFlQUAD8  ::ourMemPool(sizeof(FFlQUAD8)  ,FFlMemPoolMgr::instance());
FFaMemPool FFlTET4   ::ourMemPool(sizeof(FFlTET4)   ,FFlMemPoolMgr::instance());
FFaMemPool FFlTET10  ::ourMemPool(sizeof(FFlTET10)  ,FFlMemPoolMgr::instance());
FFaMemPool FFlWEDG6  ::ourMemPool(sizeof(FFlWEDG6)  ,FFlMemPoolMgr::instance());
FFaMemPool FFlWEDG15 ::ourMemPool(sizeof(FFlWEDG15) ,FFlMemPoolMgr::instance());
FFaMemPool FFlHEX8   ::ourMemPool(sizeof(FFlHEX8)   ,FFlMemPoolMgr::instance());
FFaMemPool FFlHEX20  ::ourMemPool(sizeof(FFlHEX20)  ,FFlMemPoolMgr::instance());
FFaMemPool FFlRBAR   ::ourMemPool(sizeof(FFlRBAR)   ,FFlMemPoolMgr::instance());
FFaMemPool FFlRGD    ::ourMemPool(sizeof(FFlRGD)    ,FFlMemPoolMgr::instance());
FFaMemPool FFlWAVGM  ::ourMemPool(sizeof(FFlWAVGM)  ,FFlMemPoolMgr::instance());
#ifdef FT_USE_STRAINCOAT
FFaMemPool FFlSTRCT3 ::ourMemPool(sizeof(FFlSTRCT3) ,FFlMemPoolMgr::instance());
FFaMemPool FFlSTRCQ4 ::ourMemPool(sizeof(FFlSTRCQ4) ,FFlMemPoolMgr::instance());
FFaMemPool FFlSTRCT6 ::ourMemPool(sizeof(FFlSTRCT6) ,FFlMemPoolMgr::instance());
FFaMemPool FFlSTRCQ8 ::ourMemPool(sizeof(FFlSTRCQ8) ,FFlMemPoolMgr::instance());
#endif
#endif


void FFlMemPool::deleteAllLinkMemPools()
{
#ifdef FT_USE_MEMPOOL
  FFlNode::freePool();
#ifdef FT_USE_VERTEX
  FFlVertex::freePool();
#endif
  FFlMemPoolMgr::removeInstance();
#endif
}


void FFlMemPool::setAsMemPoolPart(FFlLinkHandler*
#ifdef FT_USE_MEMPOOL
                                  link
#endif
                                  ){
#ifdef FT_USE_MEMPOOL
  FFlNode   ::usePartOfPool(link);
#ifdef FT_USE_VERTEX
  FFlVertex ::usePartOfPool(link);
#endif
  FFlCMASS  ::usePartOfPool(link);
  FFlSPRING ::usePartOfPool(link);
  FFlRSPRING::usePartOfPool(link);
  FFlBUSH   ::usePartOfPool(link);
  FFlBEAM2  ::usePartOfPool(link);
  FFlBEAM3  ::usePartOfPool(link);
  FFlTRI3   ::usePartOfPool(link);
  FFlTRI6   ::usePartOfPool(link);
  FFlQUAD4  ::usePartOfPool(link);
  FFlQUAD8  ::usePartOfPool(link);
  FFlTET4   ::usePartOfPool(link);
  FFlTET10  ::usePartOfPool(link);
  FFlWEDG6  ::usePartOfPool(link);
  FFlWEDG15 ::usePartOfPool(link);
  FFlHEX8   ::usePartOfPool(link);
  FFlHEX20  ::usePartOfPool(link);
  FFlRBAR   ::usePartOfPool(link);
  FFlRGD    ::usePartOfPool(link);
  FFlWAVGM  ::usePartOfPool(link);
#ifdef FT_USE_STRAINCOAT
  FFlSTRCT3 ::usePartOfPool(link);
  FFlSTRCQ4 ::usePartOfPool(link);
  FFlSTRCT6 ::usePartOfPool(link);
  FFlSTRCQ8 ::usePartOfPool(link);
#endif
#endif
}


void FFlMemPool::freeMemPoolPart(FFlLinkHandler*
#ifdef FT_USE_MEMPOOL
                                 link
#endif
                                 ){
#ifdef FT_USE_MEMPOOL
  FFlNode   ::freePartOfPool(link);
#ifdef FT_USE_VERTEX
  FFlVertex ::freePartOfPool(link);
#endif
  FFlCMASS  ::freePartOfPool(link);
  FFlSPRING ::freePartOfPool(link);
  FFlRSPRING::freePartOfPool(link);
  FFlBUSH   ::freePartOfPool(link);
  FFlBEAM2  ::freePartOfPool(link);
  FFlBEAM3  ::freePartOfPool(link);
  FFlTRI3   ::freePartOfPool(link);
  FFlTRI6   ::freePartOfPool(link);
  FFlQUAD4  ::freePartOfPool(link);
  FFlQUAD8  ::freePartOfPool(link);
  FFlTET4   ::freePartOfPool(link);
  FFlTET10  ::freePartOfPool(link);
  FFlWEDG6  ::freePartOfPool(link);
  FFlWEDG15 ::freePartOfPool(link);
  FFlHEX8   ::freePartOfPool(link);
  FFlHEX20  ::freePartOfPool(link);
  FFlRBAR   ::freePartOfPool(link);
  FFlRGD    ::freePartOfPool(link);
  FFlWAVGM  ::freePartOfPool(link);
#ifdef FT_USE_STRAINCOAT
  FFlSTRCT3 ::freePartOfPool(link);
  FFlSTRCQ4 ::freePartOfPool(link);
  FFlSTRCT6 ::freePartOfPool(link);
  FFlSTRCQ8 ::freePartOfPool(link);
#endif
#endif
}


void FFlMemPool::resetMemPoolPart()
{
#ifdef FT_USE_MEMPOOL
  FFlNode   ::useDefaultPartOfPool();
#ifdef FT_USE_VERTEX
  FFlVertex ::useDefaultPartOfPool();
#endif
  FFlCMASS  ::useDefaultPartOfPool();
  FFlSPRING ::useDefaultPartOfPool();
  FFlRSPRING::useDefaultPartOfPool();
  FFlBUSH   ::useDefaultPartOfPool();
  FFlBEAM2  ::useDefaultPartOfPool();
  FFlBEAM3  ::useDefaultPartOfPool();
  FFlTRI3   ::useDefaultPartOfPool();
  FFlTRI6   ::useDefaultPartOfPool();
  FFlQUAD4  ::useDefaultPartOfPool();
  FFlQUAD8  ::useDefaultPartOfPool();
  FFlTET4   ::useDefaultPartOfPool();
  FFlTET10  ::useDefaultPartOfPool();
  FFlWEDG6  ::useDefaultPartOfPool();
  FFlWEDG15 ::useDefaultPartOfPool();
  FFlHEX8   ::useDefaultPartOfPool();
  FFlHEX20  ::useDefaultPartOfPool();
  FFlRBAR   ::useDefaultPartOfPool();
  FFlRGD    ::useDefaultPartOfPool();
  FFlWAVGM  ::useDefaultPartOfPool();
#ifdef FT_USE_STRAINCOAT
  FFlSTRCT3 ::useDefaultPartOfPool();
  FFlSTRCQ4 ::useDefaultPartOfPool();
  FFlSTRCT6 ::useDefaultPartOfPool();
  FFlSTRCQ8 ::useDefaultPartOfPool();
#endif
#endif
}
