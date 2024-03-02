// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlRGD.H"
#include "FFlLib/FFlFEParts/FFlNode.H"


bool FFlRGDTopSpec::allowSlvAttach = true;


FFlRGD::FFlRGD(int id) : FFlElementBase(id)
{
  myRGDElemTopSpec = new FFlRGDTopSpec;
  myRGDElemTopSpec->setNodeCount(0);
  myRGDElemTopSpec->setNodeDOFs(6);
  myRGDElemTopSpec->myExplEdgePattern = 0xfCfC; // 1111 1100 1111 1100
}


FFlRGD::FFlRGD(const FFlRGD& obj) : FFlElementBase(obj)
{
  myRGDElemTopSpec = new FFlRGDTopSpec(*obj.getFEElementTopSpec());
}


FFlRGD::~FFlRGD()
{
  if (myNodes.size() > 1)
    for (NodeCIter it = this->nodesBegin()+1; it != nodesEnd(); ++it)
      if (it->isResolved())
        (*it)->setStatus(FFlNode::INTERNAL);

  delete myRGDElemTopSpec;
}


FFlFEElementTopSpec* FFlRGD::getFEElementTopSpec() const
{
  if ((int)myNodes.size() != myRGDElemTopSpec->getNodeCount())
  {
    myRGDElemTopSpec->setNodeCount(myNodes.size());
    myRGDElemTopSpec->myExplicitEdges.clear();
    for (size_t i = 2; i <= myNodes.size(); i++)
      myRGDElemTopSpec->addExplicitEdge(EdgeType(1,(int)i));
  }

  return myRGDElemTopSpec;
}


bool FFlRGD::setNode(const int topPos, FFlNode* aNode)
{
  if (topPos < 1)
    return false;
  else if ((size_t)topPos > myNodes.size())
    myNodes.resize(topPos);

  FFlFEElementTopSpec* topSpec = this->getFEElementTopSpec();
  myNodes[topPos-1] = aNode;
  aNode->pushDOFs(topSpec->getNodeDOFs(topPos));
  if (!FFlRGDTopSpec::allowSlvAttach && topPos > 1)
    aNode->setStatus(FFlNode::SLAVENODE);

  return true;
}


bool FFlRGD::setNode(const int topPos, int nodeID)
{
  if (topPos < 1)
    return false;
  else if ((size_t)topPos > myNodes.size())
    myNodes.resize(topPos);

  myNodes[topPos-1] = nodeID;

  return true;
}


bool FFlRGD::setNodes(const std::vector<int>& nodeRefs, int offset)
{
  if (offset < 0) return false;

  if (offset + nodeRefs.size() > myNodes.size())
    myNodes.resize(offset + nodeRefs.size());

  for (size_t i = 0; i < nodeRefs.size(); i++)
    myNodes[offset+i] = nodeRefs[i];

  return true;
}


bool FFlRGD::setNodes(const std::vector<FFlNode*>& nodeRefs, int offset)
{
  if (offset < 0) return false;

  if (offset + nodeRefs.size() > myNodes.size())
    myNodes.resize(offset + nodeRefs.size());

  FFlFEElementTopSpec* topSpec = this->getFEElementTopSpec();
  for (size_t i = 0; i < nodeRefs.size(); i++)
  {
    myNodes[offset+i] = nodeRefs[i];
    nodeRefs[i]->pushDOFs(topSpec->getNodeDOFs(offset+i+1));
    if (!FFlRGDTopSpec::allowSlvAttach && offset+i > 0)
      nodeRefs[i]->setStatus(FFlNode::SLAVENODE);
  }

  return true;
}


void FFlRGD::setMasterNode(FFlNode* nodeRef)
{
  this->setNode(1, nodeRef);
}


void FFlRGD::setMasterNode(int nodeRef)
{
  this->setNode(1, nodeRef);
}


void FFlRGD::addSlaveNode(FFlNode* nodeRef)
{
  myNodes.push_back(nodeRef);
}


void FFlRGD::addSlaveNode(int nodeRef)
{
  myNodes.push_back(nodeRef);
}


void FFlRGD::addSlaveNodes(const std::vector<int>& nodeRefs)
{
  this->setNodes(nodeRefs, 1);
}


void FFlRGD::addSlaveNodes(const std::vector<FFlNode*>& nodeRefs)
{
  this->setNodes(nodeRefs, 1);
}


FFlNode* FFlRGD::getMasterNode() const
{
  return myNodes.empty() ? 0 : myNodes[0].getReference();
}


void FFlRGD::getSlaveNodes(std::vector<FFlNode*>& nodeRefs) const
{
  if (myNodes.size() > 1)
    for (NodeCIter it = this->nodesBegin()+1; it != nodesEnd(); ++it)
      nodeRefs.push_back(it->getReference());
}


void FFlRGD::init()
{
  FFlRGDTypeInfoSpec::instance()->setTypeName("RGD");
  FFlRGDTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::CONSTRAINT_ELM);

  ElementFactory::instance()->registerCreator(FFlRGDTypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlRGD::create,int,FFlElementBase*&));

  FFlRGDAttributeSpec::instance()->addLegalAttribute("PRGD", false);
}
