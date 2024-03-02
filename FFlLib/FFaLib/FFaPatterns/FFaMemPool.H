// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFA_MEM_POOL_H
#define FFA_MEM_POOL_H

#include <stddef.h> // Needed for size_t definition on some platforms

#include <vector>
#include <list>
#include <map>

class FFaMemPoolMgr;


/*!
  \brief Class to streamline efficient memory allocation

  \details
  This class is used to help the programmer to reimplement the operators
  new and delete for small classes that will be instantiated a lot of times.
  To use, simply make a static member of FFaMemPool in your class, reimplement
  the operators new and delete in that class by forwarding the calls to
  FFaMemPool::doNew() and FFaMemPool::doDelete(), respectively:

  \code
  #include "FFaLib/FFaPatterns/FFaMemPool.H"

  class foo
  {
  public:
    static void* operator new(size_t n) { return ourMemPool.doNew(n); }
    static void operator delete(void* p, size_t n) { ourMemPool.doDelete(p,n); }
    static void freePool() { ourMemPool.freePool(); }

  private:
    static FFaMemPool ourMemPool;
  };
  \endcode

  In the implementation file :

  \code
  FFaMemPool foo::ourMemPool(sizeof(foo));
  \endcode

  \author Jacob Storen

*///////////////////////////////////////////////////////////////////////////////

class FFaMemPool
{
public:
  FFaMemPool(size_t objSize, FFaMemPoolMgr* mgr = NULL, int blockSize = 512);
  FFaMemPool(const FFaMemPool& mempool);
  ~FFaMemPool();

  void* doNew(size_t objSize);
  void  doDelete(void* p, size_t objSize);

  void  freePool(bool release = false);
  void  freePartOfPool(void* objPtrAsId);

  void  usePartOfPool(void* objPtrAsId);
  void  useDefaultPartOfPool() { myCurrentPP = &myDefaultPP; }

private:
  struct PoolPart
  {
    PoolPart() { myHeadOfFreeList = NULL; }
    ~PoolPart() { this->freePoolPart(); }

    void freePoolPart();

    void* myHeadOfFreeList;
    std::vector<void*> myMemBlocks;

  private:
    PoolPart(const PoolPart&) = delete;
    PoolPart& operator=(const PoolPart&) = delete;
  };

  PoolPart* myCurrentPP;
  PoolPart  myDefaultPP;

  std::map<void*,PoolPart*> myPoolParts;

  int    myBlockSize;
  size_t myObjSize;

  FFaMemPoolMgr* myManager;

  FFaMemPool& operator=(const FFaMemPool&) { return *this; }
};


class FFaMemPoolMgr
{
public:
  ~FFaMemPoolMgr() { this->freeMemPools(true); }
  void freeMemPools(bool release = false);
  void remove(FFaMemPool* p) { myMemPools.remove(p); }
  void insert(FFaMemPool* p) { if (p) myMemPools.push_front(p); }

private:
  std::list<FFaMemPool*> myMemPools;
};


#define FFA_MAKE_MEMPOOL public: \
  static void* operator new(size_t s) { return ourMemPool.doNew(s); } \
  static void operator delete(void* p, size_t s) { ourMemPool.doDelete(p,s); } \
  static void freePool()              { ourMemPool.freePool(); } \
  static void freePartOfPool(void* p) { ourMemPool.freePartOfPool(p); } \
  static void usePartOfPool(void* p)  { ourMemPool.usePartOfPool(p); } \
  static void useDefaultPartOfPool()  { ourMemPool.useDefaultPartOfPool();} \
private: \
  static FFaMemPool ourMemPool

#endif
