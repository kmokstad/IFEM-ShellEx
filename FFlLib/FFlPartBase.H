// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PART_BASE_H
#define FFL_PART_BASE_H

#include <vector>
#include <algorithm>

class FFlTypeInfoSpec;
class FFaCheckSum;

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


/*!
  \brief Base class for all object types of a FE model.
*/

class FFlPartBase
{
protected:
#ifdef FFL_REFCOUNT
  FFlPartBase(int ID) : myID(ID), refCount(0) {}
  FFlPartBase(const FFlPartBase& other) : myID(other.myID), refCount(0) {}
#else
  FFlPartBase(int ID) : myID(ID) {}
  FFlPartBase(const FFlPartBase&) = default;
#endif
  FFlPartBase& operator=(const FFlPartBase&) = delete;

public:
  virtual ~FFlPartBase() = default;

  virtual FFlTypeInfoSpec* getTypeInfoSpec() const = 0;

  int getID() const { return myID; }
  void setID(int newID) { myID = newID; }

#ifdef FFL_REFCOUNT
  unsigned char getRefCount() const { return refCount; }
  unsigned char ref() { return ++refCount; }
  bool unref() { if (refCount) --refCount; return refCount > 0; }
#else
  unsigned char getRefCount() const { return 1; }
  unsigned char ref() { return 1; }
  bool unref() { return true; }
#endif

protected:
  void checksum(FFaCheckSum* cs, int csMask = 0) const;

private:
  int myID;
#ifdef FFL_REFCOUNT
  unsigned char refCount;
#endif
};


struct FFlFEPartBaseLess
{
  bool operator()(const FFlPartBase* lhs, const FFlPartBase* rhs) const
  { return lhs->getID() < rhs->getID(); }

  bool operator()(const FFlPartBase* lhs, int rhs) const
  { return lhs->getID() < rhs; }

  bool operator()(int lhs, const FFlPartBase* rhs) const
  { return lhs < rhs->getID(); }
};


template <class T> struct FFlPartFinder
{
  using RefIt = typename std::vector<T*>::const_iterator;
  static RefIt findObject(const std::vector<T*>& objs, int ID)
  {
    std::pair<RefIt,RefIt> p = std::equal_range(objs.begin(), objs.end(),
                                                ID, FFlFEPartBaseLess());
    return p.first == p.second ? objs.end() : p.first;
  }
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
