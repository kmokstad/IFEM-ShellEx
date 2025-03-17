// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include <algorithm>

#include "FFlLib/FFlFENodeRefs.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include "FFaLib/FFaAlgebra/FFaCheckSum.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


FFlFENodeRefs::FFlFENodeRefs(const FFlFENodeRefs& obj)
{
  mySize = 0.0;
  myNodes.reserve(obj.myNodes.size());
  for (const NodeRef& node : obj.myNodes)
    myNodes.push_back(NodeRef(node.getID()));
}


FaVec3 FFlFENodeRefs::getNodeCenter() const
{
  FaVec3 center;
  for (const NodeRef& node : myNodes)
    if (node.isResolved())
      center += node->getPos();

  center /= myNodes.size();
  return center;
}


/*!
  This method calculates the size of an element,
  as the diameter of the smallest subscribing sphere.
*/

double FFlFENodeRefs::getSize() const
{
  if (mySize > 0.0)
    return mySize; // using cached value

  int numNod = 0;
  FaVec3 Xmax, Xmin;
  for (const NodeRef& node : myNodes)
    if (!node.isResolved())
      return 0.0;
    else if (++numNod == 1)
      Xmin = Xmax = node->getPos();
    else for (int i = 0; i < 3; i++)
    {
      double x = node->getPos()[i];
      if (x < Xmin[i])
	Xmin[i] = x;
      else if (x > Xmax[i])
	Xmax[i] = x;
    }

  return mySize = (Xmax-Xmin).length();
}


void FFlFENodeRefs::checksum(FFaCheckSum* cs) const
{
  for (const NodeRef& node : myNodes)
    cs->add(node.getID());
}


void FFlFENodeRefs::initNodeVector() const
{
  if (myNodes.empty())
    myNodes.resize(this->getFEElementTopSpec()->getNodeCount());
}


FFlNode* FFlFENodeRefs::getNode(const int topPos) const
{
  this->initNodeVector();

  if (topPos < 1 || topPos > (int)myNodes.size())
  {
    ListUI <<"\n *** Error: Node "<< topPos <<" is out of range.\n";
    return NULL;
  }

  return myNodes[topPos-1].getReference();
}


int FFlFENodeRefs::getNodeID(const int topPos) const
{
  this->initNodeVector();

  if (topPos > (int)myNodes.size())
  {
    ListUI <<"\n *** Error: Node "<< topPos <<" is out of range.\n";
    return 0;
  }

  return myNodes[topPos-1].getID();
}


int FFlFENodeRefs::getTopPos(const int nodeID) const
{
  for (size_t topPos = 0; topPos < myNodes.size(); topPos++)
    if (myNodes[topPos].getID() == nodeID)
      return 1+topPos;

  return 0;
}


bool FFlFENodeRefs::setNode(const int topPos, FFlNode* aNode)
{
  this->initNodeVector();

  if (topPos < 1 || topPos > (int)myNodes.size())
  {
    ListUI <<"\n *** Error: Node "<< topPos <<" is out of range.\n";
    return false;
  }

  myNodes[topPos-1] = aNode;
  aNode->pushDOFs(this->getFEElementTopSpec()->getNodeDOFs(topPos));

  return true;
}


bool FFlFENodeRefs::setNode(const int topPos, int nodeID)
{
  this->initNodeVector();

  if (topPos < 1 || topPos > (int)myNodes.size())
  {
    ListUI <<"\n *** Error: Node "<< topPos <<" is out of range.\n";
    return false;
  }

  myNodes[topPos-1] = nodeID;

  return true;
}


bool FFlFENodeRefs::setNodes(const std::vector<int>& nodeRefs,
                             size_t offset, bool)
{
  this->initNodeVector();

  size_t lastPos = offset + nodeRefs.size();
  if (lastPos > myNodes.size())
  {
    ListUI <<"\n *** Error: Node "<< lastPos-1 <<" is out of range.\n";
    return false;
  }

  for (size_t i = 0; i < nodeRefs.size(); i++)
    myNodes[offset+i] = nodeRefs[i];

  return true;
}


bool FFlFENodeRefs::setNodes(const std::vector<FFlNode*>& nodeRefs,
                             size_t offset, bool)
{
  this->initNodeVector();

  size_t lastPos = offset + nodeRefs.size();
  if (lastPos > myNodes.size())
  {
    ListUI <<"\n *** Error: Node "<< lastPos-1 <<" is out of range.\n";
    return false;
  }

  FFlFEElementTopSpec* topSpec = this->getFEElementTopSpec();
  for (size_t i = 0; i < nodeRefs.size(); i++)
  {
    myNodes[offset+i] = nodeRefs[i];
    nodeRefs[i]->pushDOFs(topSpec->getNodeDOFs(offset+i+1));
  }

  return true;
}


/*!
  Resolves references to nodes.
  Uses the \a possibleReferences range for resolving.
*/

bool FFlFENodeRefs::resolveNodeRefs(const std::vector<FFlNode*>& possibleReferences,
				    bool suppressErrmsg)
{
  if (possibleReferences.empty())
  {
    ListUI <<"\n *** Error: No nodes!\n";
    return false;
  }

  bool retVar = true;
  FFlFEElementTopSpec* topSpec = this->getFEElementTopSpec();

  int localNode = 0;
  for (NodeRef& node : myNodes)
    if (node.resolve(possibleReferences))
    {
      FFlNode* n = node.getReference(); // might still be null
      if (!n) break;

      n->pushDOFs(topSpec->getNodeDOFs(++localNode));
      if (topSpec->isReferenceNode(localNode))
        n->setStatus(FFlNode::REFNODE); // This node cannot be external
      else if (topSpec->isSlaveNode(localNode))
        n->setStatus(FFlNode::SLAVENODE); // This node cannot be external
    }
    else
    {
      retVar = false;
      if (!suppressErrmsg)
        ListUI <<"\n *** Error: Failed to resolve node "<< node.getID()
               <<" (local node "<< ++localNode <<").\n";
    }

  if (retVar && (size_t)localNode < myNodes.size())
  {
    if (!suppressErrmsg)
      ListUI <<"\n  ** Warning: Only "<< localNode
             <<" element nodes defined (expected "<< myNodes.size() <<").\n";
    myNodes.resize(localNode);
  }

  return retVar;
}


/*!
  Returns the local face number defined by one or two global node numbers.
  Used by the Nastran Bulk Data File parser.
*/

short int FFlFENodeRefs::getFaceNum(int node1, int node2) const
{
  if (node1 > 0)
    if (!(node1 = this->getTopPos(node1)))
      return 0;

  if (node2 > 0)
    if (!(node2 = this->getTopPos(node2)))
      return 0;

  return this->getFEElementTopSpec()->getFaceNum(node1,node2);
}


/*!
  Returns the number of FE nodes on the given local face.
*/

int FFlFENodeRefs::getFaceSize(short int face) const
{
  return this->getFEElementTopSpec()->getFaceVertexCount(face-1,false);
}


/*!
  Returns the FE nodes on the given local face.
*/

bool FFlFENodeRefs::getFaceNodes(std::vector<FFlNode*>& nodes, short int face,
				 bool switchNormal) const
{
  std::vector<int> topology;
  if (!this->getFEElementTopSpec()->getFaceTopology(face-1,false,switchNormal,0,topology))
    return false;

  nodes.clear();
  nodes.reserve(topology.size());
  for (int nodeIdx : topology)
    nodes.push_back(this->getNode(nodeIdx));

  return std::find(nodes.begin(),nodes.end(),(FFlNode*)0) == nodes.end();
}


int FFlFENodeRefs::getNodalCoor(double* X, double* Y, double* Z) const
{
  size_t inod = 0;
  for (const NodeRef& node : myNodes)
    if (node.isResolved())
    {
      const FaVec3& pos = node->getPos();
      X[inod] = pos.x();
      Y[inod] = pos.y();
      Z[inod] = pos.z();
      inod++;
    }
    else
      ListUI <<" *** Error: Element node "<< node.getID() <<" not resolved\n";

  return inod == myNodes.size() ? 0 : -3;
}

#ifdef FF_NAMESPACE
} // namespace
#endif
