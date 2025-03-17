// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlRGD.H"
#include "FFlLib/FFlFEParts/FFlNode.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif

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
  if (static_cast<int>(myNodes.size()) != myRGDElemTopSpec->getNodeCount())
  {
    myRGDElemTopSpec->setNodeCount(myNodes.size());
    myRGDElemTopSpec->myExplicitEdges.clear();
    for (size_t i = 2; i <= myNodes.size(); i++)
      myRGDElemTopSpec->addExplicitEdge(EdgeType(1,i));
  }

  return myRGDElemTopSpec;
}


bool FFlRGD::setNode(const int topPos, FFlNode* aNode)
{
  if (topPos < 1)
    return false;
  else if (topPos > static_cast<int>(myNodes.size()))
    myNodes.resize(topPos);

  myNodes[topPos-1] = aNode;

  aNode->pushDOFs(this->getFEElementTopSpec()->getNodeDOFs(topPos));
  if (!FFlRGDTopSpec::allowSlvAttach && topPos > 1)
    aNode->setStatus(FFlNode::SLAVENODE);

  return true;
}


bool FFlRGD::setNode(const int topPos, int nodeID)
{
  if (topPos < 1)
    return false;
  else if (topPos > static_cast<int>(myNodes.size()))
    myNodes.resize(topPos);

  myNodes[topPos-1] = nodeID;

  return true;
}


bool FFlRGD::setNodes(const std::vector<int>& nodeIDs,
                      size_t offset, bool shrink)
{
  if (shrink || offset + nodeIDs.size() > myNodes.size())
    myNodes.resize(offset + nodeIDs.size());

  for (size_t i = 0; i < nodeIDs.size(); i++)
    myNodes[offset+i] = nodeIDs[i];

  return true;
}


bool FFlRGD::setNodes(const std::vector<FFlNode*>& nodes,
                      size_t offset, bool shrink)
{
  if (shrink || offset + nodes.size() > myNodes.size())
    myNodes.resize(offset + nodes.size());

  FFlFEElementTopSpec* topSpec = this->getFEElementTopSpec();
  for (size_t i = 0; i < nodes.size(); i++)
  {
    myNodes[offset+i] = nodes[i];
    nodes[i]->pushDOFs(topSpec->getNodeDOFs(offset+i+1));
    if (!FFlRGDTopSpec::allowSlvAttach && offset+i > 0)
      nodes[i]->setStatus(FFlNode::SLAVENODE);
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
  using TypeInfoSpec  = FFaSingelton<FFlTypeInfoSpec,FFlRGD>;
  using AttributeSpec = FFaSingelton<FFlFEAttributeSpec,FFlRGD>;

  TypeInfoSpec::instance()->setTypeName("RGD");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::CONSTRAINT_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlRGD::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PRGD",false);
}

#ifdef FF_NAMESPACE
} // namespace
#endif
