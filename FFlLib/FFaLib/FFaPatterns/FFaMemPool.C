// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFaLib/FFaPatterns/FFaMemPool.H"


FFaMemPool::FFaMemPool(size_t objSize, FFaMemPoolMgr* mgr, int blockSize)
{
  myObjSize   = objSize;
  myBlockSize = blockSize;
  myManager   = mgr;
  myCurrentPP = &myDefaultPP;

  if (myManager)
    myManager->insert(this);
}

FFaMemPool::FFaMemPool(const FFaMemPool& mempool)
{
  myObjSize   = mempool.myObjSize;
  myBlockSize = mempool.myBlockSize;
  myManager   = mempool.myManager;
  myCurrentPP = &myDefaultPP;

  if (myManager)
    myManager->insert(this);
}

FFaMemPool::~FFaMemPool()
{
  if (myManager)
    myManager->remove(this);
  this->freePool();
}


void* FFaMemPool::doNew(size_t size)
{
  if (size != myObjSize)
    return ::operator new(size);

  void* p = myCurrentPP->myHeadOfFreeList;

  if (p)
    {
      void** pContents = (void**)p;
      myCurrentPP->myHeadOfFreeList = pContents[0];
    }
  else
    {
      void* newBlock_void = ::operator new(myBlockSize*myObjSize);
      if (!newBlock_void)
        return newBlock_void;

      myCurrentPP->myMemBlocks.push_back(newBlock_void);

      char* newBlock_char = (char*)newBlock_void;
      void** item = NULL;
      for (size_t i = myObjSize; i < (myBlockSize-1)*myObjSize; i += myObjSize)
      {
        item = (void**)&newBlock_char[i];
        item[0] = (void*)&newBlock_char[i+myObjSize];
      }

      item = (void**)&newBlock_char[(myBlockSize-1)*myObjSize];
      item[0] = NULL;

      p = newBlock_void;
      myCurrentPP->myHeadOfFreeList = (void*)(&newBlock_char[myObjSize]);
    }

  return p;
}


void FFaMemPool::doDelete(void* deadObject, size_t size)
{
  if (deadObject == NULL) return;

  if (size != myObjSize) {
    ::operator delete(deadObject);
    return;
  }

  void** deadObjContent = (void**)deadObject;
  deadObjContent[0] = myCurrentPP->myHeadOfFreeList;
  myCurrentPP->myHeadOfFreeList = deadObject;
}


void FFaMemPool::freePool(bool release)
{
  myDefaultPP.freePoolPart();
  std::map<void*,PoolPart*>::iterator it;
  for (it = myPoolParts.begin(); it != myPoolParts.end(); ++it)
    delete it->second;

  myPoolParts.clear();
  if (release) myManager = NULL;
}


void FFaMemPool::freePartOfPool(void* objPtrAsId)
{
  std::map<void*,PoolPart*>::iterator it = myPoolParts.find(objPtrAsId);
  if (it != myPoolParts.end())
    delete it->second;

  myPoolParts.erase(it);
}


void FFaMemPool::usePartOfPool(void* objPtrAsId)
{
  std::map<void*,PoolPart*>::iterator it = myPoolParts.find(objPtrAsId);
  if (it == myPoolParts.end())
    myPoolParts[objPtrAsId] = myCurrentPP = new PoolPart();
  else
    myCurrentPP = it->second;
}


void FFaMemPool::PoolPart::freePoolPart()
{
  for (void* block : myMemBlocks)
    if (block)
      ::operator delete(block);

  myHeadOfFreeList = NULL;
  std::vector<void*> empty;
  myMemBlocks.swap(empty);
}


void FFaMemPoolMgr::freeMemPools(bool release)
{
  for (FFaMemPool* pool : myMemPools)
    pool->freePool(release);
}
