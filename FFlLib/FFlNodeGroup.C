// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlNodeGroup.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlLinkCSMask.H"
#include "FFaLib/FFaAlgebra/FFaCheckSum.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"
#include <algorithm>

#if FFL_DEBUG > 2
#include <iostream>
#endif


/*!
  \class FFlNodeGroup FFlNodeGroup.H
  \brief Class for grouping of nodes
*/


/*!
  Constructor. Creates an empty group width id \e id and name \e gName
*/

FFlNodeGroup::FFlNodeGroup(int id, const std::string& gName) : FFlGroupBase(id,gName)
{
  iAmSorted = true;
}


void FFlNodeGroup::init()
{
  FFlNodeGroupSpec::instance()->setTypeName("NodeGroup");
  FFlNodeGroupSpec::instance()->setCathegory(FFlTypeInfoSpec::USER_DEF_GROUP);
}


/*!
  Adds a node to the group.
*/

void FFlNodeGroup::addNode(FFlNode* aNode, bool sortOnInsert)
{
  // Inserts if unique
  if (sortOnInsert && this->hasNode(aNode->getID()))
    return;

  myNodes.push_back(aNode);
#if FFL_DEBUG > 2
  std::cout <<"Node "<< aNode->getID()
            <<" added to Group "<< this->getID() << std::endl;
#endif
  iAmSorted = false;
  if (sortOnInsert) this->sortElements();
}


/*!
  Adds a node to the group.
*/

void FFlNodeGroup::addElement(int aNodeID, bool sortOnInsert)
{
  // Inserts if unique
  if (sortOnInsert && this->hasNode(aNodeID))
    return;

  myNodes.push_back(aNodeID);
#if FFL_DEBUG > 2
  std::cout <<"Node "<< aNodeID
            <<" added to Group "<< this->getID() << std::endl;
#endif
  iAmSorted = false;
  if (sortOnInsert) this->sortElements();
}


/*!
  Replaces one node by a list of new nodes.
*/

void FFlNodeGroup::swapNode(int oldElmID, const std::vector<int>& newElmID)
{
  if (this->removeNode(oldElmID))
    for (int elmID : newElmID)
      this->addNode(elmID);
}


/*!
  Removes a node from the group.
*/

bool FFlNodeGroup::remove(const GroupNodeRef& elmRef)
{
  this->sortElements();

  std::pair<GroupNodeVec::iterator,GroupNodeVec::iterator> ep;
  ep = std::equal_range(myNodes.begin(), myNodes.end(), elmRef);

  if (ep.first == ep.second)
    return false;

  myNodes.erase(ep.first, ep.second);
  return true;
}


/*!
  Resolves the node references.
  Uses the \a possibleRefs range for resolving.
*/

bool FFlNodeGroup::resolveNodeRefs(std::vector<FFlNode*>& possibleRefs, bool silence)
{
  bool retVar = true;
  for (GroupNodeRef& gnode : myNodes)
    if (!gnode.resolve(possibleRefs))
    {
      if (!silence)
        ListUI <<"\n *** Error: Invalid node Id "<< gnode.getID() <<"\n";
      retVar = false;
    }

  return retVar;
}


bool FFlNodeGroup::hasNode(int nodeID) const
{
  const_cast<FFlNodeGroup*>(this)->sortElements();

  return std::binary_search(myNodes.begin(), myNodes.end(), GroupNodeRef(nodeID));
}


void FFlNodeGroup::sortElements(bool removeDuplicates)
{
  if (!iAmSorted)
  {
    std::sort(myNodes.begin(), myNodes.end());
    iAmSorted = true;
  }

  if (!removeDuplicates || myNodes.size() < 2)
    return;

  // Remove duplicated nodes, if any
  size_t i, n = 0;
  for (i = 1; i < myNodes.size(); i++)
    if (myNodes[i].getID() > myNodes[n].getID())
      if (++n < i) myNodes[n] = myNodes[i];

  myNodes.resize(1+n);
}


void FFlNodeGroup::calculateChecksum(FFaCheckSum* cs, int csMask) const
{
  if ((csMask & FFl::CS_GROUPMASK) != FFl::CS_NOGROUPINFO)
  {
    for (const GroupNodeRef& gnode : myNodes)
      cs->add(gnode.getID());
    FFlNamedPartBase::checksum(cs, csMask);
  }
}
