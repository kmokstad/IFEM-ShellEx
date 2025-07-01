// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlWAVGM.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


FFlWAVGM::FFlWAVGM(int id) : FFlElementBase(id)
{
  myWAVGMElemTopSpec = new FFlWAVGMTopSpec;
  myWAVGMElemTopSpec->setNodeCount(0);
  myWAVGMElemTopSpec->setNodeDOFs(0);
  myWAVGMElemTopSpec->myExplEdgePattern = 0xeeee; // 1110 1110 1110 1110
}


FFlWAVGM::FFlWAVGM(const FFlWAVGM& obj) : FFlElementBase(obj)
{
  myWAVGMElemTopSpec = new FFlWAVGMTopSpec(*obj.getFEElementTopSpec());
}


FFlWAVGM::~FFlWAVGM()
{
  delete myWAVGMElemTopSpec;
}


FFlFEElementTopSpec* FFlWAVGM::getFEElementTopSpec() const
{
  if (static_cast<int>(myNodes.size()) != myWAVGMElemTopSpec->getNodeCount())
  {
    myWAVGMElemTopSpec->setNodeCount(myNodes.size());
    myWAVGMElemTopSpec->myExplicitEdges.clear();
    for (size_t i = 2; i <= myNodes.size(); i++)
      myWAVGMElemTopSpec->addExplicitEdge(1,i);
  }

  return myWAVGMElemTopSpec;
}


bool FFlWAVGM::setNode(const int topPos, FFlNode* aNode)
{
  if (topPos < 1)
    return false;
  else if (topPos > static_cast<int>(myNodes.size()))
    myNodes.resize(topPos);

  myNodes[topPos-1] = aNode;

  aNode->pushDOFs(this->getFEElementTopSpec()->getNodeDOFs(topPos));
  if (topPos == 1) aNode->setStatus(FFlNode::REFNODE);

  return true;
}


bool FFlWAVGM::setNode(const int topPos, int nodeID)
{
  if (topPos < 1)
    return false;
  else if (topPos > static_cast<int>(myNodes.size()))
    myNodes.resize(topPos);

  myNodes[topPos-1] = nodeID;

  return true;
}


bool FFlWAVGM::setNodes(const std::vector<int>& nodeRefs,
                        size_t offset, bool)
{
  if (offset + nodeRefs.size() > myNodes.size())
    myNodes.resize(offset + nodeRefs.size());

  for (size_t i = 0; i < nodeRefs.size(); i++)
    myNodes[offset+i] = nodeRefs[i];

  return true;
}


bool FFlWAVGM::setNodes(const std::vector<FFlNode*>& nodeRefs,
                        size_t offset, bool)
{
  if (offset + nodeRefs.size() > myNodes.size())
    myNodes.resize(offset + nodeRefs.size());

  FFlFEElementTopSpec* topSpec = this->getFEElementTopSpec();
  for (size_t i = 0; i < nodeRefs.size(); i++)
  {
    myNodes[offset+i] = nodeRefs[i];
    nodeRefs[i]->pushDOFs(topSpec->getNodeDOFs(offset+i+1));
    if (offset+i == 0) nodeRefs[i]->setStatus(FFlNode::REFNODE);
  }

  return true;
}


void FFlWAVGM::setSlaveNode(FFlNode* nodeRef)
{
  this->setNode(1, nodeRef);
}


void FFlWAVGM::setSlaveNode(int nodeRef)
{
  this->setNode(1, nodeRef);
}


void FFlWAVGM::addMasterNode(FFlNode* nodeRef)
{
  myNodes.push_back(nodeRef);
}


void FFlWAVGM::addMasterNode(int nodeRef)
{
  myNodes.push_back(nodeRef);
}


void FFlWAVGM::addMasterNodes(const std::vector<int>& nodeRefs)
{
  this->setNodes(nodeRefs, 1);
}


void FFlWAVGM::addMasterNodes(const std::vector<FFlNode*>& nodeRefs)
{
  this->setNodes(nodeRefs, 1);
}


FFlNode* FFlWAVGM::getSlaveNode() const
{
  return myNodes.size() > 0 ? myNodes.front().getReference() : 0;
}


void FFlWAVGM::getMasterNodes(std::vector<FFlNode*>& nodeRefs) const
{
  if (myNodes.size() > 1)
    for (NodeCIter it = this->nodesBegin()+1; it != nodesEnd(); ++it)
      nodeRefs.push_back(it->getReference());
}


void FFlWAVGM::init()
{
  using TypeInfoSpec  = FFaSingelton<FFlTypeInfoSpec,FFlWAVGM>;
  using AttributeSpec = FFaSingelton<FFlFEAttributeSpec,FFlWAVGM>;

  TypeInfoSpec::instance()->setTypeName("WAVGM");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::CONSTRAINT_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlWAVGM::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PWAVGM",false);
}


bool FFlWAVGM::removeMasterNodes(std::vector<int>& nodeRefs)
{
  bool ok = true;
  std::vector<int> nodeIdx(nodeRefs.size(),0);
  for (size_t i = 0; i < nodeRefs.size(); i++)
    if ((nodeIdx[i] = this->getTopPos(nodeRefs[i])) == 1)
    {
      ListUI <<"\n *** Error: Can not remove the reference node "<< nodeRefs[i]
             <<" from WAVGM "<< this->getID();
      ok = false;
    }
    else if (nodeIdx[i] < 1)
    {
      ListUI <<"\n *** Error: Node "<< nodeRefs[i]
             <<" is not connected to WAVGM "<< this->getID();
      ok = false;
    }
    else
      --nodeIdx[i];

  if (!ok)
  {
    ListUI <<"\n";
    return false;
  }

  std::sort(nodeIdx.begin(),nodeIdx.end());
  for (int j = nodeIdx.size()-1; j >= 0; j--)
    myNodes.erase(myNodes.begin()+nodeIdx[j]);

  nodeRefs.swap(nodeIdx);
  return true;
}

#ifdef FF_NAMESPACE
} // namespace
#endif
