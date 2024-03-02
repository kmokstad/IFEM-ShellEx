// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlGroup.H"
#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlLinkCSMask.H"
#include "FFaLib/FFaAlgebra/FFaCheckSum.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"
#include <algorithm>

#if FFL_DEBUG > 2
#include <iostream>
#endif


/*!
  \class FFlGroup FFlGroup.H
  \brief Class for grouping of elements
*/


/*!
  Constructor. Creates an empty group width id \e id and name \e groupName
*/

FFlGroup::FFlGroup(int id, const std::string& groupName) : FFlNamedPartBase(id)
{
  this->setName(groupName);
  iAmSorted = true;
}


void FFlGroup::init()
{
  FFlGroupTypeInfoSpec::instance()->setTypeName("Group");
  FFlGroupTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::USER_DEF_GROUP);
}


/*!
  Adds an element to the group.
*/

void FFlGroup::addElement(FFlElementBase* anElement, bool sortOnInsert)
{
  // Inserts if unique
  if (sortOnInsert && this->hasElement(anElement->getID()))
    return;

  myElements.push_back(anElement);
#if FFL_DEBUG > 2
  std::cout <<"Element "<< anElement->getID()
            <<" added to Group "<< this->getID() << std::endl;
#endif
  iAmSorted = false;
  if (sortOnInsert) this->sortElements();
}


/*!
  Adds an element to the group.
*/

void FFlGroup::addElement(int anElementID, bool sortOnInsert)
{
  // Inserts if unique
  if (sortOnInsert && this->hasElement(anElementID))
    return;

  myElements.push_back(anElementID);
#if FFL_DEBUG > 2
  std::cout <<"Element "<< anElementID
            <<" added to Group "<< this->getID() << std::endl;
#endif
  iAmSorted = false;
  if (sortOnInsert) this->sortElements();
}


/*!
  Replaces one element by a list of new elements.
*/

void FFlGroup::swapElement(int oldElmID, const std::vector<int>& newElmID)
{
  if (this->removeElement(oldElmID))
    for (int elmID : newElmID)
      this->addElement(elmID);
}


/*!
  Removes an element from the group.
*/

bool FFlGroup::remove(const GroupElemRef& elmRef)
{
  this->sortElements();

  std::pair<GroupElemVec::iterator,GroupElemVec::iterator> ep;
  ep = std::equal_range(myElements.begin(), myElements.end(), elmRef);

  if (ep.first == ep.second)
    return false;

  myElements.erase(ep.first, ep.second);
  return true;
}


/*!
  Resolves the element references.
  Uses the \a possibleReferences range for resolving.
*/

bool FFlGroup::resolveElemRefs(std::vector<FFlElementBase*>& possibleReferences,
			       bool suppressErrmsg)
{
  bool retVar = true;
  for (GroupElemRef& gelem : myElements)
    if (!gelem.resolve(possibleReferences))
    {
      if (!suppressErrmsg)
	ListUI <<"\n *** Error: Invalid element Id "<< gelem.getID() <<"\n";
      retVar = false;
    }

  return retVar;
}


bool FFlGroup::hasElement(int elementID) const
{
  const_cast<FFlGroup*>(this)->sortElements();

  return std::binary_search(myElements.begin(), myElements.end(), GroupElemRef(elementID));
}


void FFlGroup::sortElements(bool removeDuplicates)
{
  if (!iAmSorted)
  {
    std::sort(myElements.begin(), myElements.end());
    iAmSorted = true;
  }

  if (!removeDuplicates || myElements.size() < 2)
    return;

  // Remove duplicated elements, if any
  size_t i, n = 0;
  for (i = 1; i < myElements.size(); i++)
    if (myElements[i].getID() > myElements[n].getID())
      if (++n < i) myElements[n] = myElements[i];

  myElements.resize(1+n);
}


void FFlGroup::calculateChecksum(FFaCheckSum* cs, int csMask) const
{
  if ((csMask & FFl::CS_GROUPMASK) != FFl::CS_NOGROUPINFO)
  {
    for (const GroupElemRef& gelem : myElements)
      cs->add(gelem.getID());
    FFlNamedPartBase::checksum(cs, csMask);
  }
}
