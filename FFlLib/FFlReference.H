// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef _FFLREFERENCE_H_
#define _FFLREFERENCE_H_


#include <map>
#include <vector>
#include <algorithm>
#include <utility>

#include "FFlLib/FFlPartBase.H"


/*!
 * \NOTE Forward added
 * \TODO Check if sematic is as wanted with friend defs
 */
template <class T> class FFlReference;
template <class T> bool operator<(const FFlReference<T>& lhs, const FFlReference<T>& rhs);


template <class T>
class FFlReference
{
public:
  FFlReference() : iAmResolved(true)
  {
    myResolvedRef = 0;
  }

  FFlReference(T* objPtr) : iAmResolved(true)
  {
    myResolvedRef = objPtr;
  }

  FFlReference(int objID) : iAmResolved(false)
  {
    ID = objID;
  }

  FFlReference(const FFlReference& ref) : iAmResolved(ref.iAmResolved)
  {
    if (iAmResolved)
      myResolvedRef = ref.myResolvedRef;
    else
      ID = ref.ID;
  }

  virtual ~FFlReference() {}


  FFlReference& operator=(const FFlReference& ref)
  {
    if (this != &ref)
    {
      if (iAmResolved)
        if (ref.iAmResolved)
	  myResolvedRef = ref.myResolvedRef;
        else
        {
	  ID = ref.ID;
	  iAmResolved = false;
        }
      else
	ID = ref.getID();
    }

    return *this;
  }


  FFlReference& operator=(const T* objPtr) // fast resolving manually
  {
    myResolvedRef = (T*)objPtr;
    iAmResolved = true;

    return *this;
  }


  FFlReference& operator=(const int& objID) // fast assignment manualy
  {
    ID = objID;
    iAmResolved = false;

    return *this;
  }


  bool isResolved() const { return iAmResolved && myResolvedRef != 0; }


  typedef std::map<int,T*>                                resolveRefMapType;
  typedef typename resolveRefMapType::const_iterator resolveRefMapTypeCIt;

  bool resolve(const resolveRefMapType& possibleReferences)
  {
    if (iAmResolved) return true;

    resolveRefMapTypeCIt it = possibleReferences.find(ID);
    if (it != possibleReferences.end())
    {
      myResolvedRef = it->second;
      iAmResolved = true;
    }

    return iAmResolved;
  }


  typedef std::vector<T*>                                 resolveRefVecType;
  typedef typename resolveRefVecType::const_iterator resolveRefVecTypeCIt;

  bool resolve(const resolveRefVecType& possibleReferences)
  {
    if (iAmResolved) return true;

    std::pair<resolveRefVecTypeCIt,resolveRefVecTypeCIt> ep;
    ep = equal_range(possibleReferences.begin(),
		     possibleReferences.end(),
		     ID, FFlFEPartBaseLess());

    if (ep.first != ep.second)
    {
      myResolvedRef = *ep.first;
      iAmResolved = true;
    }

    return iAmResolved;
  }


  /*!
    Unresolves the reference - used before e.g. copying of objects
  */

  void unresolve()
  {
    if (!iAmResolved) return;

    iAmResolved = false;
    if (myResolvedRef == 0) return;

    int tmpID = myResolvedRef->getID();
    ID = tmpID;
  }


  int getID() const
  {
    return iAmResolved && myResolvedRef != 0 ? myResolvedRef->getID() : ID;
  }


  T* getReference() const
  {
    if (iAmResolved) return myResolvedRef;

    char* p = 0; *p = 'd'; // force a core dump (internal logic error)
    return 0;
  }

  T* operator->() const
  {
    if (iAmResolved && myResolvedRef != 0) return myResolvedRef;

    char* p = 0; *p = 'd'; // force a core dump (internal logic error)
    return 0;
  }

private:
  bool iAmResolved;

  union {
    T* myResolvedRef;
    int ID;
  };

#if defined(_MSC_VER) && _MSC_VER < 1310
  friend bool operator< (const FFlReference<T>& lhs, const FFlReference<T>& rhs);
#else
  friend bool operator< <T> (const FFlReference<T>& lhs, const FFlReference<T>& rhs);
#endif
};


template <class T>
bool operator<(const FFlReference<T>& lhs, const FFlReference<T>& rhs)
{
  if (lhs.isResolved())
    if (rhs.isResolved())
      return lhs.myResolvedRef->getID() < rhs.myResolvedRef->getID();
    else
      return lhs.myResolvedRef->getID() < rhs.ID;
  else
    if (rhs.isResolved())
      return lhs.ID < rhs.myResolvedRef->getID();
    else
      return lhs.ID < rhs.ID;
}

#endif
