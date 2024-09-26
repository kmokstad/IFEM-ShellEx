// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_NODE_GROUP_H
#define FFL_NODE_GROUP_H

#include "FFlLib/FFlGroupBase.H"
#include "FFlLib/FFlReference.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFaLib/FFaPatterns/FFaSingelton.H"

class FFlNode;
class FFaCheckSum;

typedef FFlReference<FFlNode>     GroupNodeRef;
typedef std::vector<GroupNodeRef> GroupNodeVec;


class FFlNodeGroup : public FFlGroupBase
{
public:
  FFlNodeGroup(int id, const std::string& groupName = "New group");
  virtual ~FFlNodeGroup() {}

  static void init();

  bool resolveNodeRefs(std::vector<FFlNode*>& possibleRefs, bool silence = false);

  virtual void addElement(int ID, bool sortOnInsert = false);
  void addNode(FFlNode* nodeRef, bool sortOnInsert = false);
  void addNode(int ID) { this->addElement(ID); }
  void swapNode(int oldID, const std::vector<int>& newID);
  bool removeNode(int nodeID) { return this->remove(nodeID); }
  bool removeNode(FFlNode* node) { return this->remove(node); }

  bool hasNode(int nodeID) const;

  virtual void sortElements(bool removeDuplicates = false);
  size_t size() const { return myNodes.size(); }

  GroupNodeVec::const_iterator begin() const { return myNodes.begin(); }
  GroupNodeVec::const_iterator end()   const { return myNodes.end(); }

  virtual FFlGroupBase* clone() const { return new FFlNodeGroup(*this); }

  virtual void calculateChecksum(FFaCheckSum* cs, int csMask) const;

  FFL_TYPE_INFO(FFlNodeGroup);

protected:
  bool remove(const GroupNodeRef& nodeRef);

private:
  GroupNodeVec myNodes;
  bool         iAmSorted;
};

typedef FFaSingelton<FFlTypeInfoSpec,FFlNodeGroup> FFlNodeGroupSpec;

#endif
