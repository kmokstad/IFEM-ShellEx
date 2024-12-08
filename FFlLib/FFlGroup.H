// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_GROUP_H
#define FFL_GROUP_H

#include "FFlLib/FFlNamedPartBase.H"
#include "FFlLib/FFlReference.H"
#include "FFlLib/FFlTypeInfoSpec.H"

class FFaCheckSum;

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif

class FFlElementBase;

using GroupElemRef = FFlReference<FFlElementBase>;
using GroupElemVec = std::vector<GroupElemRef>;


class FFlGroup : public FFlNamedPartBase
{
public:
  FFlGroup(int id, const std::string& groupName = "New group");
  virtual ~FFlGroup() {}

  static void init();

  bool resolveElemRefs(std::vector<FFlElementBase*>& possibleReferences,
		       bool suppressErrmsg = false);

  void addElement(FFlElementBase* elementRef, bool sortOnInsert = false);
  void addElement(int elementID, bool sortOnInsert = false);
  void swapElement(int oldElmID, const std::vector<int>& newElmID);
  bool removeElement(int elementID) { return this->remove(elementID); }
  bool removeElement(FFlElementBase* element) { return this->remove(element); }

  bool hasElement(int elementID) const;

  void sortElements(bool removeDuplicates = false);
  size_t size() const { return myElements.size(); }

  GroupElemVec::const_iterator begin() const { return myElements.begin(); }
  GroupElemVec::const_iterator end()   const { return myElements.end(); }

  FFlGroup* clone() const { return new FFlGroup(*this); }

  void calculateChecksum(FFaCheckSum* cs, int csMask) const;

  FFL_TYPE_INFO(FFlGroup);

protected:
  bool remove(const GroupElemRef& elmRef);

private:
  GroupElemVec myElements;
  bool         iAmSorted;
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
