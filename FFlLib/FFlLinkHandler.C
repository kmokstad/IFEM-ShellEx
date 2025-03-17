// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <cfloat>

#include "FFlLib/FFlLinkHandler.H"
#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlPartBase.H"
#ifdef FT_USE_VERTEX
#include "FFlLib/FFlVertex.H"
#endif
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlGroup.H"
#include "FFlLib/FFlMemPool.H"
#include "FFlLib/FFlLoadBase.H"
#include "FFlLib/FFlAttributeBase.H"
#ifdef FT_USE_VISUALS
#include "FFlLib/FFlFEParts/FFlVDetail.H"
#include "FFlLib/FFlFEParts/FFlVAppearance.H"
#endif
#include "FFlLib/FFlFEParts/FFlPCOORDSYS.H"
#include "FFlLib/FFlFEParts/FFlPWAVGM.H"
#include "FFlLib/FFlFEParts/FFlWAVGM.H"
#include "FFlLib/FFlFEParts/FFlRGD.H"
#include "FFlLib/FFlFEResultBase.H"
#ifdef FT_USE_CONNECTORS
#include "FFlLib/FFlConnectorItems.H"
#include "FFaLib/FFaGeometry/FFaCompoundGeometry.H"
#endif
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"
#include "FFaLib/FFaAlgebra/FFaCheckSum.H"
#include "FFaLib/FFaAlgebra/FFaTensor3.H"
#include "FFaLib/FFaAlgebra/FFaMat34.H"
#include "FFaLib/FFaAlgebra/FFaMath.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"
#ifdef FFL_TIMER
#include "FFaLib/FFaProfiler/FFaProfiler.H"
#endif

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


FFlLinkHandler::FFlLinkHandler(size_t maxNodes, size_t maxElms)
{
#ifdef FFL_TIMER
  myProfiler = new FFaProfiler("LinkProfiler");
#else
  myProfiler = NULL;
#endif
  myResults = NULL;
  myChkSum  = NULL;
  nGenDofs  = 0;
  nodeLimit = maxNodes;
  elmLimit  = maxElms;

  areElementsSorted = areNodesSorted = areLoadsSorted = true;
#ifdef FT_USE_VISUALS
  areVisualsSorted = true;
#endif
  tooLarge = hasLooseNodes = isResolved = false;
}


FFlLinkHandler::FFlLinkHandler(const FFlLinkHandler& otherLink)
{
#ifdef FFL_TIMER
  myProfiler = new FFaProfiler("LinkProfiler");
#else
  myProfiler = NULL;
#endif
  myResults = NULL;
  myChkSum  = NULL;
  nGenDofs  = otherLink.nGenDofs;
  tooLarge  = otherLink.tooLarge;
  nodeLimit = otherLink.nodeLimit;
  elmLimit  = otherLink.elmLimit;

  // copy constructor:
  // 1 - assign new elements using the clone method (prototype factory pattern)
  // 2 - resolve with the new containers

  areElementsSorted = otherLink.areElementsSorted;
  myElements.reserve(otherLink.myElements.size());
  for (FFlElementBase* elm : otherLink.myElements)
    myElements.push_back(elm->clone());

  areNodesSorted = otherLink.areNodesSorted;
  myNodes.reserve(otherLink.myNodes.size());
  for (FFlNode* node : otherLink.myNodes)
  {
    myNodes.push_back(node->clone());
#ifdef FT_USE_VERTEX
    this->addVertex(myNodes.back()->getVertex());
#endif
  }

  for (const GroupMap::value_type& g : otherLink.myGroupMap)
    myGroupMap[g.second->getID()] = g.second->clone();

  areLoadsSorted = otherLink.areLoadsSorted;
  myLoads.reserve(otherLink.myLoads.size());
  for (FFlLoadBase* load : otherLink.myLoads)
    myLoads.push_back(load->clone());

  for (const AttributeTypeMap::value_type& am : otherLink.myAttributes)
    for (const AttributeMap::value_type& attr : am.second)
      this->addAttribute(attr.second->clone());

#ifdef FT_USE_VISUALS
  areVisualsSorted = otherLink.areVisualsSorted;
  myVisuals.reserve(otherLink.myVisuals.size());
  for (FFlVisualBase* vis : otherLink.myVisuals)
    myVisuals.push_back(vis->clone());
#endif

  hasLooseNodes = isResolved = false;
  this->resolve();
}


FFlLinkHandler::FFlLinkHandler(const FFlGroup& fromGroup)
{
#ifdef FFL_TIMER
  myProfiler = new FFaProfiler("LinkProfiler");
#else
  myProfiler = NULL;
#endif
  myResults = NULL;
  myChkSum  = NULL;
  nodeLimit = elmLimit = nGenDofs = 0;
  tooLarge  = hasLooseNodes = isResolved = false;
  areElementsSorted = areNodesSorted = areLoadsSorted = true;
#ifdef FT_USE_VISUALS
  areVisualsSorted = true;
#endif

  std::set<FFlNode*,FFlFEPartBaseLess> tmpNodes;
  AttributeTypeMap                     tmpMap;

  myElements.reserve(fromGroup.size());
  for (const GroupElemRef& e : fromGroup)
  {
    myElements.push_back(e->clone());
    for (NodeCIter nit = e->nodesBegin(); nit != e->nodesEnd(); ++nit)
      tmpNodes.insert(nit->getReference());
    for (AttribsCIter a = e->attributesBegin(); a != e->attributesEnd(); ++a)
    {
      FFlAttributeBase* attr = a->second.getReference();
      tmpMap[attr->getTypeName()][attr->getID()] = attr;
    }
  }

  myNodes.reserve(tmpNodes.size());
  for (FFlNode* node : tmpNodes)
  {
    myNodes.push_back(node->clone());
#ifdef FT_USE_VERTEX
    this->addVertex(myNodes.back()->getVertex());
#endif
  }

  for (const AttributeTypeMap::value_type& am : tmpMap)
    for (const AttributeMap::value_type& attr : am.second)
      this->addAttribute(attr.second->clone());

  this->resolve();
}


FFlLinkHandler::~FFlLinkHandler()
{
  this->deleteResults();
  this->deleteGeometry();

#ifdef FFL_TIMER
  myProfiler->report();
  delete myProfiler;
#endif
  delete myChkSum;
}


void FFlLinkHandler::deleteGeometry()
{
  FFlMemPool::setAsMemPoolPart(this);

#ifdef FT_USE_VERTEX
  for (FaVec3* vtx : myVertices)
    static_cast<FFlVertex*>(vtx)->unRef();
#endif

  for (const GroupMap::value_type& g : myGroupMap) delete g.second;
  for (FFlLoadBase*    l : myLoads   ) delete l;
  for (FFlElementBase* e : myElements) delete e;
  for (FFlNode*        n : myNodes   ) delete n;
#ifdef FT_USE_VISUALS
  for (FFlVisualBase*  d : myVisuals ) delete d;
#endif
  for (const AttributeTypeMap::value_type& am : myAttributes)
    for (const AttributeMap::value_type& attr : am.second)
      delete attr.second;

  myNumElements.clear();
  myElements.clear();
  myFElements.clear();
  myBushElements.clear();
  myNodes.clear();
  myFEnodes.clear();
#ifdef FT_USE_VERTEX
  myVertices.clear();
#endif
  myGroupMap.clear();
  myLoads.clear();
  myAttributes.clear();
  uniqueAtts.clear();
#ifdef FT_USE_VISUALS
  myVisuals.clear();
#endif
  ext2intNode.clear();

  FFlMemPool::resetMemPoolPart();
  FFlMemPool::freeMemPoolPart(this);

  tooLarge = hasLooseNodes = false;

  areElementsSorted = areNodesSorted = areLoadsSorted = true;
#ifdef FT_USE_VISUALS
  areVisualsSorted = true;
#endif

  isResolved = true;
}


/*!
  Checksum calculation on "all" entries.
*/

unsigned int FFlLinkHandler::calculateChecksum(int csType, bool rndOff) const
{
  if (!myChkSum)
    myChkSum = new FFaCheckSum();

  myChkSum->reset();
  this->calculateChecksum(myChkSum,csType,rndOff);
  return myChkSum->getCurrent();
}


/*!
  Checksum calculation on "important" entries. This include node, element,
  load and attribute data. Group and visualization data are not included here.
*/

void FFlLinkHandler::calculateChecksum(FFaCheckSum* cs, bool rndOff) const
{
  int csType = FFl::CS_NOGROUPINFO | FFl::CS_NOSTRCINFO | FFl::CS_NOVISUALINFO;
  this->calculateChecksum(cs,csType,rndOff);
}


void FFlLinkHandler::calculateChecksum(FFaCheckSum* cs,
                                       int csType, bool rndOff) const
{
  if (!cs) return;

#ifdef FFL_DEBUG
  std::cout <<"\nFFlLinkHandler::calculateChecksum(csType="<< csType <<")\n";
#endif
  bool checkStrainCoat = (csType & FFl::CS_STRCMASK) != FFl::CS_NOSTRCINFO;

  cs->reset();

  for (FFlElementBase* elm : myElements)
    if (checkStrainCoat || !isStrainCoat(elm))
    {
      elm->calculateChecksum(cs,csType);
#if FFL_DEBUG > 2
      std::cout <<"Link checksum after element "<< elm->getID()
                <<" : "<< cs->getCurrent() << std::endl;
#endif
    }

#ifdef FFL_DEBUG
  std::cout <<"Link checksum after elements: "<< cs->getCurrent() << std::endl;
#endif

  bool checkExtNodeInfo = (csType & FFl::CS_EXTMASK) != FFl::CS_NOEXTINFO;
#ifdef FFL_DEBUG
  std::cout <<"Checking nodes (checkExtNodeInfo="<< std::boolalpha
            << checkExtNodeInfo <<", rndOff="<< std::boolalpha
            << rndOff <<"):"<< std::endl;
#endif

  for (FFlNode* node : myNodes)
  {
    node->calculateChecksum(cs, rndOff ? 10 : 0, checkExtNodeInfo);
#if FFL_DEBUG > 2
    std::cout <<"Link checksum after node "<< node->getID();
#if FFL_DEBUG > 3
    std::cout <<", coord = "<< std::setprecision(16) << node->getPos();
#endif
    std::cout <<" : "<< cs->getCurrent() << std::endl;
#endif
  }

#ifdef FFL_DEBUG
  std::cout <<"Link checksum after nodes: "<< cs->getCurrent() << std::endl;
#endif

  for (FFlLoadBase* load : myLoads)
  {
    load->calculateChecksum(cs,csType);
#if FFL_DEBUG > 2
    std::cout <<"Link checksum after "<< load->getTypeName()
              <<" "<< load->getID() <<" : "<< cs->getCurrent() << std::endl;
#endif
  }

#ifdef FFL_DEBUG
  std::cout <<"Link checksum after loads: "<< cs->getCurrent() << std::endl;
  std::cout <<"Checking attributes (csType="<< csType <<"):"<< std::endl;
#endif

  // Lambda function checking for strain coat properties
  auto&& isStrainCoatProp = [](const FFlAttributeBase* att)
  {
    return att->getTypeInfoSpec()->getCathegory() == FFlTypeInfoSpec::STRC_PROP;
  };

  for (const AttributeTypeMap::value_type& am : myAttributes)
    for (const AttributeMap::value_type& attr : am.second)
      if (checkStrainCoat || !isStrainCoatProp(attr.second))
      {
        attr.second->calculateChecksum(cs,csType);
#if FFL_DEBUG > 1
        std::cout <<"Link checksum after attribute "<< attr.first
                  <<" "<< attr.second->getID() <<" : "
                  << cs->getCurrent() << std::endl;
#endif
      }

#ifdef FFL_DEBUG
  std::cout <<"Link checksum after attributes: "
            << cs->getCurrent() << std::endl;
#endif

  for (const GroupMap::value_type& g : myGroupMap)
  {
    g.second->calculateChecksum(cs,csType);
#if FFL_DEBUG > 2
    std::cout <<"Link checksum after group "<< g.second->getID() <<" : "
              << cs->getCurrent() << std::endl;
#endif
  }

#ifdef FFL_DEBUG
  std::cout <<"Link checksum after groups: "<< cs->getCurrent() << std::endl;
#endif

#ifdef FT_USE_VISUALS
  if ((csType & FFl::CS_VISUALMASK) != FFl::CS_NOVISUALINFO)
    for (FFlVisualBase* vis : myVisuals)
    {
      vis->calculateChecksum(cs);
#if FFL_DEBUG > 2
      std::cout <<"Link checksum after detail "<< vis->getID()
                <<" : "<< cs->getCurrent() << std::endl;
#endif
    }
#endif

#ifdef FFL_DEBUG
  std::cout <<"Link checksum : "<< cs->getCurrent() << std::endl;
#endif
}


bool FFlLinkHandler::isStrainCoat(FFlElementBase* elm)
{
  return elm ? elm->getCathegory() == FFlTypeInfoSpec::STRC_ELM : false;
}


void FFlLinkHandler::convertUnits(const FFaUnitCalculator* convCal)
{
  if (!convCal) return;

  for (FFlNode* node :  myNodes)
    node->convertUnits(convCal);

  for (FFlLoadBase* load : myLoads)
    load->convertUnits(convCal);

  for (const AttributeTypeMap::value_type& am : myAttributes)
    for (const AttributeMap::value_type& attr : am.second)
      attr.second->convertUnits(convCal);
}


void FFlLinkHandler::initiateCalculationFlag(bool status)
{
  for (FFlElementBase* elm : myElements)
    elm->setUpForCalculations(status);
}


bool FFlLinkHandler::updateCalculationFlag(int groupID, bool status)
{
  FFlGroup* group = this->getGroup(groupID);
  if (!group)
  {
    ListUI <<" *** Error: Non-existing element group "<< groupID <<" ignored\n";
    return false;
  }

  for (const GroupElemRef& elm : *group)
    elm->setUpForCalculations(status);

  return true;
}


bool FFlLinkHandler::updateCalculationFlag(FFlPartBase* part, bool status)
{
  FFlGroup* tmpGroup = dynamic_cast<FFlGroup*>(part);
  if (tmpGroup) {
    for (const GroupElemRef& elm : *tmpGroup)
      elm->setUpForCalculations(status);
    return true;
  }

  FFlAttributeBase* attr = dynamic_cast<FFlAttributeBase*>(part);
  if (attr) {
    for (FFlElementBase* elm : myElements)
      if (elm->hasAttribute(attr))
        elm->setUpForCalculations(status);
    return true;
  }

  FFlElementBase* elm = dynamic_cast<FFlElementBase*>(part);
  if (elm) {
    elm->setUpForCalculations(status);
    return true;
  }

  return false;
}


bool FFlLinkHandler::updateCalculationFlag(const std::string& type,
                                           int id, bool status)
{
  FFlAttributeBase* attrib = this->getAttribute(type,id);
  if (attrib)
    return this->updateCalculationFlag(attrib,status);
  else
    return false;
}


#ifdef FT_USE_VISUALS
void FFlLinkHandler::updateGroupVisibilityStatus()
{
  // first reset all:
  for (const GroupMap::value_type& g : myGroupMap)
    g.second->resetVisibilityStatus();

  // run over all attributes:
  for (const AttributeTypeMap::value_type& am : myAttributes)
    for (const AttributeMap::value_type& attr : am.second)
      attr.second->resetVisibilityStatus();

  // run over all elements:
  for (FFlElementBase* e : myElements)
  {
    int visStat = e->isVisible() ?
      FFlNamedPartBase::FFL_HAS_VIS_ELM :
      FFlNamedPartBase::FFL_HAS_HIDDEN_ELM;

    for (AttribsCIter a = e->attributesBegin(); a != e->attributesEnd(); ++a)
    {
      a->second->addVisibilityStatus(FFlNamedPartBase::FFL_USED);
      a->second->addVisibilityStatus(visStat);
    }
  }

  // run over all groups:
  for (const GroupMap::value_type& g : myGroupMap)
  {
    g.second->addVisibilityStatus(FFlNamedPartBase::FFL_USED);
    for (const GroupElemRef& elm : *g.second)
      g.second->addVisibilityStatus(elm->isVisible() ?
				    FFlNamedPartBase::FFL_HAS_VIS_ELM :
				    FFlNamedPartBase::FFL_HAS_HIDDEN_ELM);
  }
}


bool FFlLinkHandler::setVisDetail(const FFlVDetail* detail)
{
  for (FFlElementBase* elm : myElements)
    elm->setDetail(detail);

  return true;
}


bool FFlLinkHandler::setVisDetail(FFlPartBase* part, const FFlVDetail* detail)
{
  FFlGroup* tmpGroup = dynamic_cast<FFlGroup*>(part);
  if (tmpGroup) {
    for (const GroupElemRef& elm : *tmpGroup)
      elm->setDetail(detail);
    return true;
  }

  FFlAttributeBase* attr = dynamic_cast<FFlAttributeBase*>(part);
  if (attr) {
    for (FFlElementBase* elm : myElements)
      if (elm->hasAttribute(attr))
        elm->setDetail(detail);
    return true;
  }

  FFlElementBase* elm = dynamic_cast<FFlElementBase*>(part);
  if (elm) {
    elm->setDetail(detail);
    return true;
  }

  return false;
}


bool FFlLinkHandler::setVisDetail(const std::vector<FFlPartBase*>& parts,
                                  const FFlVDetail* detail)
{
  std::vector<FFlAttributeBase*> attribParts;

  for (FFlPartBase* part : parts)
  {
    FFlGroup* tmpGroup = dynamic_cast<FFlGroup*>(part);
    if (tmpGroup)
      for (const GroupElemRef& elm : *tmpGroup)
        elm->setDetail(detail);
    else
    {
      FFlAttributeBase* attr = dynamic_cast<FFlAttributeBase*>(part);
      if (attr)
        attribParts.push_back(attr);
      else
      {
        FFlElementBase* elm = dynamic_cast<FFlElementBase*>(part);
        if (elm)
          elm->setDetail(detail);
      }
    }
  }

  if (!attribParts.empty())
    for (FFlElementBase* elm : myElements)
      if (elm->getDetail() != detail)
        if (elm->hasAttribute(attribParts))
          elm->setDetail(detail);

  return true;
}


bool FFlLinkHandler::setVisAppearance(FFlPartBase* part,
                                      const FFlVAppearance* app)
{
  FFlGroup* tmpGroup = dynamic_cast<FFlGroup*>(part);
  if (tmpGroup) {
    for (const GroupElemRef& elm : *tmpGroup)
      elm->setAppearance(app);
    return true;
  }

  FFlAttributeBase* attr = dynamic_cast<FFlAttributeBase*>(part);
  if (attr) {
    for (FFlElementBase* elm : myElements)
      if (elm->hasAttribute(attr))
        elm->setAppearance(app);
    return true;
  }

  FFlElementBase* elm = dynamic_cast<FFlElementBase*>(part);
  if (elm) {
    elm->setAppearance(app);
    return true;
  }

  return false;
}
#endif


bool FFlLinkHandler::addElement(FFlElementBase* anElement, bool sortOnInsert)
{
  if (!anElement) return false;

  if (elmLimit > 0 && !sortOnInsert)
    if (myElements.size() >= elmLimit)
    {
      tooLarge = true;
      ListUI <<"\n *** Error: This FE model is too large!"
             <<" It has more than the allowable "<< nodeLimit <<" elements.\n";
      return false;
    }

  if (areElementsSorted && !myElements.empty())
    if (anElement->getID() <= myElements.back()->getID())
      areElementsSorted = false;

  myNumElements.clear();
  myElements.push_back(anElement);
  if (sortOnInsert && !areElementsSorted)
    this->sortElements();

  isResolved = false;
  return true;
}


bool FFlLinkHandler::addNode(FFlNode* aNode, bool sortOnInsert)
{
  if (!aNode) return false;

  if (nodeLimit > 0 && !sortOnInsert)
    if (myNodes.size() >= nodeLimit)
    {
      tooLarge = true;
      ListUI <<"\n *** Error: This FE model is too large!"
             <<" It has more than the allowable "<< nodeLimit <<" nodes.\n";
      return false;
    }

  if (areNodesSorted && !myNodes.empty())
    if (aNode->getID() <= myNodes.back()->getID())
      areNodesSorted = false;

  myNodes.push_back(aNode);
#ifdef FT_USE_VERTEX
  this->addVertex(aNode->getVertex());
#endif
  if (sortOnInsert && !areNodesSorted)
    this->sortNodes();

  isResolved = false;
  return true;
}


void FFlLinkHandler::addLoad(FFlLoadBase* load, bool sortOnInsert)
{
  if (!load) return;

  if (areLoadsSorted && !myLoads.empty())
    if (load->getID() < myLoads.back()->getID())
      areLoadsSorted = false;

  myLoads.push_back(load);
  if (sortOnInsert && !areLoadsSorted)
    this->sortLoads();

  isResolved = false;
}


#ifdef FT_USE_VISUALS
void FFlLinkHandler::addVisual(FFlVisualBase* visual, bool sortOnInsert)
{
  if (!visual) return;

  if (areVisualsSorted && !myVisuals.empty())
    if (visual->getID() < myVisuals.back()->getID())
      areVisualsSorted = false;

  myVisuals.push_back(visual);
  if (sortOnInsert && !areVisualsSorted)
    this->sortVisuals();

  isResolved = false;
}


void FFlLinkHandler::setRunningIdxOnAppearances()
{
  FFlVAppearance* vapp = NULL; int idx = 0;
  for (FFlVisualBase* vis : myVisuals)
    if ((vapp = dynamic_cast<FFlVAppearance*>(vis)))
      vapp->runningIdx = idx++;
}
#endif


/*!
  Returns the node with the given \a ID.
  Returns NULL if no such node.
*/

FFlLinkHandler::NodesIter FFlLinkHandler::getNodeIter(int ID) const
{
  if (!areNodesSorted) this->sortNodes();

  return FFlPartFinder<FFlNode>::findObject(myNodes,ID);
}

FFlNode* FFlLinkHandler::getNode(int ID) const
{
  NodesIter nit = this->getNodeIter(ID);
  return nit == myNodes.end() ? NULL : *nit;
}


FFlNode* FFlLinkHandler::getFENode(int inod) const
{
  if (inod <= 0) return NULL;

  if (hasLooseNodes)
  {
    if (myFEnodes.empty())
    {
      if (!areNodesSorted) this->sortNodes();

      // Build a vector of FE nodes only (those that are not DOF-less).
      // Note that this method is assumed used by the solvers only.
      // Therefore, we don't include any loose external nodes here.
      for (FFlNode* node : myNodes)
        if (node->hasDOFs())
          myFEnodes.push_back(node);

      if (myFEnodes.empty())
        return NULL;
    }

    if ((size_t)inod <= myFEnodes.size())
      return myFEnodes[inod-1];
  }
  else if ((size_t)inod <= myNodes.size())
    return myNodes[inod-1]; // Assuming no loose nodes exist!

  return NULL;
}


#ifdef FT_USE_VISUALS
FFlVAppearance* FFlLinkHandler::getAppearance(int ID) const
{
  if (!areVisualsSorted) this->sortVisuals();

  std::pair<VisualsCIter,VisualsCIter> ep = equal_range(myVisuals.begin(),
                                                        myVisuals.end(), ID,
                                                        FFlFEPartBaseLess());

  FFlVAppearance* vapp = NULL;
  while (ep.first != ep.second)
    if ((vapp = dynamic_cast<FFlVAppearance*>(*ep.first)))
      break;
    else
      ep.first++;

  return vapp;
}


FFlVDetail* FFlLinkHandler::getDetail(int ID) const
{
  if (!areVisualsSorted) this->sortVisuals();

  std::pair<VisualsCIter,VisualsCIter> ep = equal_range(myVisuals.begin(),
                                                        myVisuals.end(), ID,
                                                        FFlFEPartBaseLess());

  FFlVDetail* vdet = NULL;
  while (ep.first != ep.second)
    if ((vdet = dynamic_cast<FFlVDetail*>(*ep.first)))
      break;
    else
      ep.first++;

  return vdet;
}
#endif


/*!
  Returns the element with the given \a ID.
  Returns NULL if no such element.
*/

FFlLinkHandler::ElementsIter FFlLinkHandler::getElementIter(int ID) const
{
  if (!areElementsSorted) this->sortElements();

  return FFlPartFinder<FFlElementBase>::findObject(myElements,ID);
}

FFlElementBase* FFlLinkHandler::getElement(int ID, bool internalID) const
{
  if (internalID)
    return ID >= 0 && (size_t)ID < myElements.size() ? myElements[ID] : NULL;

  ElementsIter eit = this->getElementIter(ID);
  return eit == myElements.end() ? NULL : *eit;
}


/*!
  Filters out strain coat elements and optionally the stiffness/mass-giving
  elements which do not have recovery results (like constraint elements).
*/

static bool filterFiniteElements(FFlElementBase* elm, bool resultLessElements,
                                 bool keepStrainCoat = false)
{
  switch (elm->getCathegory()) {
  case FFlTypeInfoSpec::SOLID_ELM:
  case FFlTypeInfoSpec::SHELL_ELM:
  case FFlTypeInfoSpec::BEAM_ELM:
    return true; // These element types have results
  case FFlTypeInfoSpec::STRC_ELM:
    return keepStrainCoat; // Ignore strain coat elements?
  default:
    return resultLessElements;
  }
}


/*!
  Builds the vector of finite elements that contribute to the stiffness matrix.
  Note that this function must not be invoked before resolve() is called,
  or else we'll get a crash (core dump).
*/

int FFlLinkHandler::buildFiniteElementVec(bool allFElements) const
{
  int status = 0;
  if (!areElementsSorted) this->sortElements();

  myFElements.clear();
  for (FFlElementBase* elm : myElements)
    if (filterFiniteElements(elm,allFElements)) // Skip result-less elements
    {
      std::string curTyp = elm->getTypeName();

      // Check how many element nodes do have degrees of freedom
      int nelnod = 0, lerr = 0;
      for (NodeCIter n = elm->nodesBegin(); n != elm->nodesEnd(); n++)
        if ((*n)->hasDOFs())
          nelnod++;
        else if (nelnod == 0 && (curTyp == "WAVGM" || curTyp == "CMASS"))
          break; // Reference node is not connected to other elements
        else if (curTyp == "WAVGM" && !(*n)->isExternal())
        {
          status--;
          if (++lerr == 1)
            ListUI <<" *** Error : WAVGM element "<< elm->getID()
                   <<" is invalid and must be manually corrected:\n";
          ListUI <<"             Element node "<< (*n)->getID()
                 <<" is not connected to other finite elements.\n";
        }
        else
          lerr++;

      if (curTyp == "RGD" && nelnod < 2)
      {
        if (lerr > 0)
          ListUI <<"  ** Warning : RGD element "<< elm->getID() <<" has no"
                 <<" dependent nodes connected to other elements (ignored).\n";
        else
          ListUI <<"  ** Warning : One-noded RGD element "<< elm->getID()
                 <<" (ignored).\n";
      }
      else if ((curTyp == "WAVGM" && nelnod < 2 && lerr == 0) ||
               (curTyp == "CMASS" && nelnod < 1))
        ListUI <<"  ** Warning : "<< curTyp <<" element "<< elm->getID()
               <<" has no other elements connected to its reference node "
               << (*elm->nodesBegin())->getID() <<" (ignored).\n";

      else
        myFElements.push_back(elm);
    }

  return status < 0 ? status : myFElements.size();
}


/*!
  Returns the element with the given internal element number (finite elements
  only, strain coat elements and other DOF-less elements are excluded).
  Returns NULL if no such element.
*/

FFlElementBase* FFlLinkHandler::getFiniteElement(int iel) const
{
  if (iel <= 0) return NULL;

  if (myFElements.empty())
    if (this->buildFiniteElementVec() < 1)
      return NULL;

  return (size_t)iel <= myFElements.size() ? myFElements[iel-1] : NULL;
}


/*!
  Returns the element group with the given \a ID.
  Returns NULL if no such group.
*/

FFlGroup* FFlLinkHandler::getGroup(int ID) const
{
  GroupCIter git = myGroupMap.find(ID);
  return git == myGroupMap.end() ? NULL : git->second;
}


/*!
  Returns the attribute of the given \a type and with the given \a ID.
  Returns NULL if no such attribute.
*/

FFlAttributeBase* FFlLinkHandler::getAttribute(const std::string& type,
                                               int ID) const
{
  AttributeTypeCIter atit = myAttributes.find(type);
  if (atit == myAttributes.end()) return NULL;

  AttributeMap::const_iterator ait = atit->second.find(ID);
  if (ait == atit->second.end()) return NULL;

  return ait->second;
}


/*!
  Returns the set of all attributes of the given \a type.
*/

const AttributeMap& FFlLinkHandler::getAttributes(const std::string& type) const
{
  AttributeTypeCIter atit = myAttributes.find(type);
  if (atit != myAttributes.end())
    return atit->second;

  static AttributeMap empty;
  return empty;
}


/*!
  Returns all spider reference nodes in the FE model,
  which are not connected to other elements than the spider itself.
*/

bool FFlLinkHandler::getRefNodes(std::vector<FFlNode*>& refNodes) const
{
  refNodes.clear();
  if (myNodes.empty())
    return false;

  // Find the largest external node number in the model
  int Nmax = (*std::max_element(myNodes.begin(), myNodes.end(),
                                [](FFlNode* a, FFlNode* b)
                                {return a->getID() < b->getID(); }))->getID();
#ifdef FFL_DEBUG
  std::cout <<"Largest external node number: "<< Nmax << std::endl;
#endif

  // Determine the number of elements connected to each node
  std::vector<short int> elmConnections(Nmax,0);
  for (FFlElementBase* elm : myElements)
    for (NodeCIter nit = elm->nodesBegin(); nit != elm->nodesEnd(); ++nit)
      ++elmConnections[(*nit)->getID()];

  for (FFlElementBase* elm : myElements)
    if (elm->getCathegory() == FFlTypeInfoSpec::CONSTRAINT_ELM &&
        elm->getNodeCount() > 2 && elmConnections[elm->getNodeID(1)] == 1)
      refNodes.push_back(elm->getNode(1));

  return refNodes.empty();
}


/*!
  Returns all loads with the given \a ID.
  Returns \e false if no such loads.
*/

bool FFlLinkHandler::getLoads(int ID, LoadsVec& loads) const
{
  if (!areLoadsSorted) this->sortLoads();

  loads.clear();
  std::pair<LoadsCIter,LoadsCIter> ep = equal_range(myLoads.begin(),
                                                    myLoads.end(), ID,
                                                    FFlFEPartBaseLess());
  if (ep.first == ep.second) return false;

  loads.insert(loads.end(),ep.first,ep.second);

  return true;
}


/*!
  Returns a set of the external load \a IDs.
*/

void FFlLinkHandler::getLoadCases(std::set<int>& IDs) const
{
  IDs.clear();
  for (FFlLoadBase* load : myLoads)
    IDs.insert(load->getID());
}


/*!
  Returns the number of attributes of the given \a type.
*/

int FFlLinkHandler::getAttributeCount(const std::string& type) const
{
  AttributeTypeCIter atit = myAttributes.find(type);
  return atit == myAttributes.end() ? 0 : atit->second.size();
}


/*!
  Returns the number of elements of each type.
*/

const ElmTypeCount& FFlLinkHandler::getElmTypeCount() const
{
  if (myNumElements.empty())
    this->countElements();

  return myNumElements;
}


/*!
  Returns the number of elements of the given \a type.
*/

int FFlLinkHandler::getElementCount(const std::string& type) const
{
  if (myNumElements.empty())
    this->countElements();

  ElmTypeCount::const_iterator eit = myNumElements.find(type);
  return eit == myNumElements.end() ? 0 : eit->second;
}


/*!
  Returns the number of elements of the given \a types.
  If \a checkCF is true, only the elements for which
  the calculation flag is set are counted.
*/

int FFlLinkHandler::getElementCount(int types, bool checkCF) const
{
  if (!checkCF)
  {
    if (types == FFL_ALL)
      return myElements.size();
    else if (types == FFL_FEM && !myFElements.empty())
      return myFElements.size();
  }

  int nfem  = 0; // Number of finite elements (or all if type == FFL_ALL)
  int nstrc = 0; // Number of strain coat elements
  for (FFlElementBase* elm : myElements)
    if (!checkCF || elm->doCalculations())
    {
      if (types == FFL_ALL)
        nfem++;
      else if (isStrainCoat(elm))
        nstrc++;
      else if (types == FFL_FEM || types == elm->getCathegory())
        nfem++;
    }

  return (types == FFL_STRC ? nstrc : nfem);
}


/*!
  Returns the total number of DOFs in the FE model.
  Unless \a includeExternalDofs is true, only the internal dofs are counted.
*/

int FFlLinkHandler::getDofCount(bool includeExternalDofs) const
{
  int nDOF = 0;
  for (FFlNode* node : myNodes)
    if (includeExternalDofs || !node->isExternal())
      nDOF += node->getMaxDOFs();

  return nDOF;
}


/*!
  Returns the number of nodes of given \a types.
  Loose nodes, if any, are not counted.
*/

int FFlLinkHandler::getNodeCount(int types) const
{
  if (types == FFL_ALL) return myNodes.size();

  int nnod = 0; // Count only the nodes that are connected to finite elements
  for (FFlNode* node : myNodes)
    if (node->hasDOFs()) nnod++;
#if FFL_DEBUG > 1
    else std::cout <<"Ignoring loose node "<< node->getID() << std::endl;
#endif

  return nnod;
}


/*!
  This method is assumed used by the solvers only (not not the GUI).
  Therefore, we don't include any loose external nodes here.
*/

int FFlLinkHandler::getIntNodeID(int ID) const
{
  if (ext2intNode.empty())
  {
    if (!areNodesSorted) this->sortNodes();

    int nnod = 0;
    for (FFlNode* node : myNodes)
      if (node->hasDOFs()) // ommit all loose nodes
        ext2intNode[node->getID()] = ++nnod;
  }

  std::map<int,int>::const_iterator it = ext2intNode.find(ID);
  return it == ext2intNode.end() ? -1 : it->second;
}


int FFlLinkHandler::getIntElementID(int ID) const
{
  if (myFElements.empty())
    if (this->buildFiniteElementVec() < 1)
      return 0;

  ElementsIter it = FFlPartFinder<FFlElementBase>::findObject(myFElements,ID);
  return it == myFElements.end() ? 0 : (it - myFElements.begin()) + 1;
}


int FFlLinkHandler::getNewElmID() const
{
  if (myElements.empty()) return 1;
  if (!areElementsSorted) this->sortElements();
  return (*(myElements.rbegin()))->getID() + 1;
}


int FFlLinkHandler::getNewNodeID() const
{
  if (myNodes.empty()) return 1;
  if (!areNodesSorted) this->sortNodes();
  return (*(myNodes.rbegin()))->getID() + 1;
}


int FFlLinkHandler::getNewGroupID() const
{
  if (myGroupMap.empty()) return 1;
  return (--(myGroupMap.end()))->first + 1;
}


int FFlLinkHandler::getNewAttribID(const std::string& type) const
{
  AttributeTypeCIter ait = myAttributes.find(type);
  if (ait == myAttributes.end()) return 1;
  return (--(ait->second.end()))->first + 1;
}


#ifdef FT_USE_VISUALS
int FFlLinkHandler::getNewVisualID() const
{
  if (myVisuals.empty()) return 1;
  if (!areVisualsSorted) this->sortVisuals();
  return (*(myVisuals.rbegin()))->getID() + 1;
}
#endif


bool FFlLinkHandler::addGroup(FFlGroup* group, bool silence)
{
  if (!group) return false;

  if (myGroupMap.insert({group->getID(),group}).second)
  {
    isResolved = false;
    return true;
  }

  if (!silence)
    ListUI <<"\n  ** Warning: Multiple groups with ID="<< group->getID()
	   <<" detected. Only the first one is used.\n";

  delete group;
  return false;
}


bool FFlLinkHandler::addAttribute(FFlAttributeBase* attr, bool silence,
				  const std::string& name)
{
  if (!attr) return false;

  if (myAttributes[name].insert({attr->getID(),attr}).second)
  {
    isResolved = false;
    return true;
  }

  if (!silence)
    ListUI <<"\n  ** Warning: Multiple attributes with identical ID detected ("
	   << name <<" "<< attr->getID() <<"). Only the first one is used.\n";

  delete attr;
  return false;
}


bool FFlLinkHandler::addAttribute(FFlAttributeBase* attr, bool silence)
{
  return attr ? this->addAttribute(attr,silence,attr->getTypeName()) : false;
}


int FFlLinkHandler::addUniqueAttribute(FFlAttributeBase* attr, bool silence)
{
  if (!attr) return 0;

  // Check if an identical attribute already exists
  for (const AttributeMap::value_type& attp : myAttributes[attr->getTypeName()])
    if (attp.second->isIdentic(attr))
    {
      delete attr;
      return attp.first;
    }

#ifdef FFL_DEBUG
  attr->print("Unique attribute ");
#endif

  int attId = attr->getID();
  this->addAttribute(attr,silence);
  return attId;
}


int FFlLinkHandler::addUniqueAttributeCS(FFlAttributeBase*& attr)
{
  if (!myChkSum)
    myChkSum = new FFaCheckSum();

  // Calculate the checksum of the new attribute to determine its uniqueness
  myChkSum->reset();
  attr->calculateChecksum(myChkSum,FFl::CS_NOIDINFO);
  unsigned int cs = myChkSum->getCurrent();

  // Check if an identical attribute already exists
  std::map<unsigned int,FFlAttributeBase*>::iterator cit = uniqueAtts.find(cs);
  if (cit != uniqueAtts.end())
  {
    delete attr;
    attr = cit->second;
    return attr->getID();
  }

#ifdef FFL_DEBUG
  attr->print("Unique attribute ");
#endif

  uniqueAtts[cs] = attr;
  return this->addAttribute(attr) ? attr->getID() : 0;
}


bool FFlLinkHandler::removeAttribute(const std::string& typeName,
                                     int ID, bool silence)
{
  AttributeTypeMap::iterator atit = myAttributes.find(typeName);
  if (atit == myAttributes.end()) return false;

  AttributeMap::iterator atrit = atit->second.find(ID);
  if (atrit == atit->second.end()) return false;

  if (!silence)
    ListUI <<"\n   * Note: Erasing attribute "<< typeName <<" "<< ID;
  delete atrit->second;
  atit->second.erase(atrit);
  return true;
}


#ifdef FT_USE_VISUALS

/*!
  Returns a visual detail object with detail type ON or OFF. If no such
  detail exists, a new one will be created and added to the model.
*/

FFlVDetail* FFlLinkHandler::getPredefDetail(int detailType)
{
  FFlVDetail* vdet = NULL;
  for (FFlVisualBase* vis : myVisuals)
    if ((vdet = dynamic_cast<FFlVDetail*>(vis)))
      if (vdet->detail.getValue() == detailType)
	return vdet;

  vdet = new FFlVDetail(this->getNewVisualID());
  vdet->detail.setValue(detailType);
  this->addVisual(vdet,true);

  return vdet;
}


FFlVDetail* FFlLinkHandler::getOnDetail()
{
  return this->getPredefDetail(FFlVDetail::ON);
}


FFlVDetail* FFlLinkHandler::getOffDetail()
{
  return this->getPredefDetail(FFlVDetail::OFF);
}
#endif


/*!
  Returns the transformation matrix of all internal coordinate systems.
*/

void FFlLinkHandler::getAllInternalCoordSys(std::vector<FaMat34>& mxes) const
{
  mxes.clear();
  AttributeTypeCIter mapit = myAttributes.find("PCOORDSYS");
  if (mapit == myAttributes.end()) return;

  FFlPCOORDSYS* lcs;
  mxes.reserve(mapit->second.size());
  for (const AttributeMap::value_type& attr : mapit->second)
    if ((lcs = dynamic_cast<FFlPCOORDSYS*>(attr.second)))
    {
      mxes.push_back(FaMat34());
      mxes.back().makeCS_Z_XZ(lcs->Origo.getValue(),
			      lcs->Zaxis.getValue(),
			      lcs->XZpnt.getValue());
    }
}


/*!
  Returns a node with (at least) the given number of DOFs (\a dofFilter), that
  matches the given \a point. If more than one node match the point, return
  the node among these that is either a free node or a WAVGM reference node.
  If more than one such match exist, return the closest node among them.

  TODO: This method does not work well with large tolerances.
  The tolerance is used as a "all within is coincident".
  In that regime it actually is ID-dependent which node is returned.
  It is supposed to prioritize the nodes like this:
  1. Replacement node for a WAVGM reference node
  2. WAVGM reference node
  3. Closest free node
  But if a free node is found that happens to be closer than a found
  reference node, the free node will be returned instead.
  There should be two tolerances in this method:
  Coincident tolerance and search tolerance.
  JJS 2009.03.27
*/

FFlNode* FFlLinkHandler::findFreeNodeAtPoint(const FaVec3& point,
                                             double tol, int dofFilter) const
{
  FFlNode* closest = NULL;
  double closestdist = DBL_MAX;
  double sqrTol = tol > sqrt(DBL_MAX) ? DBL_MAX : tol*tol;
  double sqrdist, xd, yd, zd;

  for (FFlNode* node : myNodes)
    if (node->hasDOFs(dofFilter) || node->isExternal() || node->isRefNode())
      if ((xd = fabs(node->getPos().x() - point.x())) < tol)
        if ((yd = fabs(node->getPos().y() - point.y())) < tol)
          if ((zd = fabs(node->getPos().z() - point.z())) < tol)
            if ((sqrdist = xd*xd+yd*yd+zd*zd) < sqrTol)
            {
              if (closestdist >= sqrTol)
              {
                // This is the first matching node
                closest     = node;
                closestdist = sqrdist;
              }
              else if (node->isAttachable() && !closest->isAttachable())
              {
                // The previous node found was dependent, replace it with this
                closest     = node;
                closestdist = sqrdist;
              }
              else if (node->isRefNode() && !closest->isSlaveNode())
              {
                // The previous node found was a free node (not connected to a
                // constraint element) whereas this node is a reference node.
                // Check if they are connected via an auto-added BUSH element.
                if (!this->areBUSHconnected(closest,node))
                {
                  // Not connected, use the found reference node instead
                  closest     = node;
                  closestdist = sqrdist;
                }
              }
	      else if (!node->isSlaveNode())
	      {
		if (closest->isRefNode())
		{
		  // The previous node found was a reference node whereas this
		  // node is a free node. Check if they are BUSH-connected.
		  if (this->areBUSHconnected(node,closest))
		  {
		    // They are connected, use the found free node instead
		    closest     = node;
		    closestdist = sqrdist;
		  }
		}
		else if (node->getMaxDOFs() > closest->getMaxDOFs())
		{
		  // This node has 6 DOFs whereas the previous one has only 3.
		  // Choose the node having the most DOFs.
		  closest     = node;
		  closestdist = sqrdist;
		}
		else if (sqrdist < closestdist)
		  if (closest->getStatus(1) == node->getStatus(1))
		  {
		    // This node matches better than the previous node found
		    closest     = node;
		    closestdist = sqrdist;
		  }
              }
            }

#ifdef FFL_DEBUG
  if (!closest)
    std::cout <<" *** No FE node found at point ="<< point
              <<", tol = "<< tol <<" dofFilter ="<< dofFilter << std::endl;
#endif

  return closest;
}


/*!
  Checks if the two nodes \a n1 and \a n2 are connected via a BUSH element.
*/

bool FFlLinkHandler::areBUSHconnected(FFlNode* n1, FFlNode* n2) const
{
  if (myBushElements.empty())
    if (this->buildBUSHelementSet() < 1)
      return false;

  for (FFlElementBase* elm : myBushElements)
    if (elm->getNodeID(1) == n1->getID() && elm->getNodeID(2) == n2->getID())
      return true;

  return false;
}


int FFlLinkHandler::buildBUSHelementSet() const
{
  for (FFlElementBase* elm : myElements)
    if (elm->getTypeName() == "BUSH")
      myBushElements.insert(elm);

  return myBushElements.size();
}


/*!
  Creates a new node at \a nodePos, a BUSH element between \a fromNode and the
  new node, and a CMASS element at \a fromNode. The new node is returned.
  The elements are created without properties - only topology.
*/

FFlNode* FFlLinkHandler::createAttachableNode(FFlNode* fromNode,
					      const FaVec3& nodePos,
					      FFlConnectorItems* cItems)
{
  if (!fromNode) return NULL;

  FFlElementBase* newElm = NULL;
  FFlNode* newNode = new FFlNode(this->getNewNodeID(),nodePos);
  ListUI <<"  -> Creating FE node "<< newNode->getID() <<"\n";
  this->addNode(newNode,true);
  if (cItems)
#ifdef FT_USE_CONNECTORS
    cItems->addNode(newNode->getID());
#else
    std::cerr <<"*** FFlLinkHandler::createAttachableNode: Logic error\n";
#endif

  newElm = ElementFactory::instance()->create("BUSH",this->getNewElmID());
  newElm->setNode(1,newNode);
  newElm->setNode(2,fromNode);
  ListUI <<"  -> Creating BUSH element "<< newElm->getID()
         <<" between nodes "<< newNode->getID()
         <<" and "<< fromNode->getID() <<"\n";
  this->addElement(newElm,true);
#ifdef FT_USE_CONNECTORS
  if (cItems)
    cItems->addElement(newElm->getID());
#endif
  if (myBushElements.empty())
    this->buildBUSHelementSet();
  myBushElements.insert(newElm);

  newElm = ElementFactory::instance()->create("CMASS",this->getNewElmID());
  newElm->setNode(1,fromNode);
  ListUI <<"  -> Creating CMASS element "<< newElm->getID()
         <<" on node "<< fromNode->getID() <<"\n";
  this->addElement(newElm,true);
#ifdef FT_USE_CONNECTORS
  if (cItems)
    cItems->addElement(newElm->getID());
#endif

  return newNode;
}


/*!
  Returns the Id of the node that is closest to \a nodePos.
  The actual coordinates for the node is then returned by updating \a nodePos.
*/

int FFlLinkHandler::findNode(FaVec3& nodePos, const FaMat34* Tlg) const
{
  if (Tlg) // Transform to local coordinate system
    nodePos = Tlg->inverse()*nodePos;

  FFlNode* node = this->findClosestNode(nodePos);
  if (!node) return 0;

  if (Tlg) // Transform to global coordinate system
    nodePos = (*Tlg) * node->getPos();
  else
    nodePos = node->getPos();

  return node->getID();
}


/*!
  Returns the node that is closest to the given \a point.
*/

FFlNode* FFlLinkHandler::findClosestNode(const FaVec3& point) const
{
  FFlNode* closest = NULL;
  double closestdist = DBL_MAX;
  double sqrdist;

  for (FFlNode* node : myNodes)
    if ((sqrdist = (point - node->getPos()).sqrLength()) <= closestdist)
    {
      closest     = node;
      closestdist = sqrdist;
    }

  return closest;
}


/*!
  Returns the element within the cathegories \a wantedTypes
  that has its node center closest to the given \a point.
  If \a wantedTypes is empty, any element is accepted.
*/

FFlElementBase* FFlLinkHandler::findClosestElement(const FaVec3& point,
                                                   const CathegoryVec& wantedTypes) const
{
  FFlElementBase* closest = NULL;
  double closestdist = DBL_MAX;
  double sqrdist;

  // OTHER_ELM = the total number of element cathegories
  std::vector<bool> typeOk(FFlTypeInfoSpec::OTHER_ELM+1,wantedTypes.empty());
  for (FFlTypeInfoSpec::Cathegory type : wantedTypes)
    typeOk[type] = true;

  for (FFlElementBase* elm : myElements)
    if (typeOk[elm->getCathegory()])
      if ((sqrdist = (point - elm->getNodeCenter()).sqrLength()) <= closestdist)
      {
        closest     = elm;
        closestdist = sqrdist;
      }

  return closest;
}


/*!
  Returns the element within a specified \a group
  that has its node center closest to the given \a point.
*/

FFlElementBase* FFlLinkHandler::findClosestElement(const FaVec3& point,
                                                   const FFlGroup& group) const
{
  FFlElementBase* closest = NULL;
  double closestdist = DBL_MAX;
  double sqrdist;

  for (const GroupElemRef& elm : group)
    if ((sqrdist = (point - elm->getNodeCenter()).sqrLength()) <= closestdist)
    {
      closest     = elm.getReference();
      closestdist = sqrdist;
    }

  return closest;
}


/*!
  Convenience method to search for the closest element among elements
  in all cathegories, or within an element group.
*/

FFlElementBase* FFlLinkHandler::findClosestElement(const FaVec3& point,
                                                   FFlGroup* group) const
{
  if (group)
    return this->findClosestElement(point,*group);
  else
    return this->findClosestElement(point,{});
}


FFlElementBase* FFlLinkHandler::findPoint(const FFlGroup& group,
                                          const FaVec3& point, double* xi) const
{
  FFlElementBase* candidate = this->findClosestElement(point,group);
  if (candidate && candidate->invertMapping(point,xi))
    return candidate;

  // Check all elements with its center closer to the target point
  // than the size of the smallest subscriping ball of the element
  for (const GroupElemRef& elm : group)
    if ((point - elm->getNodeCenter()).length() < elm->getSize())
      if (elm->invertMapping(point,xi))
        return elm.getReference();

#ifdef FFL_DEBUG
  std::cerr <<" *** FFlLinkHandler::findPoint: "<< point
            <<" does not match any element in group "<< group.getID()
            << std::endl;
#endif
  return NULL;
}


FFlElementBase* FFlLinkHandler::findPoint(const FaVec3& point, double* xi,
                                          int groupID) const
{
  if (groupID <= 0)
  {
    // The invertMapping method is implemented for shell elements only.
    // Filter out all other element types from the element search.
    CathegoryVec shell({ FFlTypeInfoSpec::SHELL_ELM });
    FFlElementBase* elm = this->findClosestElement(point,shell);
    if (elm && elm->invertMapping(point,xi))
      return elm;

    // Check all shell elements with its center closer to the target point
    // than the size of the smallest subscriping ball of the element
    for (FFlElementBase* elm : myElements)
      if (elm->getCathegory() == shell.front())
        if ((point - elm->getNodeCenter()).length() < elm->getSize())
          if (elm->invertMapping(point,xi))
            return elm;

#ifdef FFL_DEBUG
    std::cerr <<" *** FFlLinkHandler::findPoint: "<< point
              <<" does not match any element."<< std::endl;
#endif
    return NULL;
  }

  FFlGroup* group = this->getGroup(groupID);
  if (!group)
  {
#ifdef FFL_DEBUG
    std::cerr <<" *** FFlLinkHandler::findPoint: Non-existing element group "
              << groupID << std::endl;
#endif
    return NULL;
  }

  return this->findPoint(*group,point,xi);
}


#ifdef FT_USE_VERTEX
FFlVertex* FFlLinkHandler::getVertex(size_t i) const
{
  return i < myVertices.size() ? static_cast<FFlVertex*>(myVertices[i]) : NULL;
}


const FFlrVxToElmMap& FFlLinkHandler::getVxToElementMapping()
{
  if (myVxMapping.empty())
  {
    myVxMapping.resize(myVertices.size());
    int n, id;
    NodeCIter nit;
    for (FFlElementBase* elm : myElements)
      if (filterFiniteElements(elm,false,true)) // Skip result-less elements
        for (nit = elm->nodesBegin(), n = 1; nit != elm->nodesEnd(); ++nit, n++)
          if ((id = (*nit)->getVertexID()) >= 0)
            myVxMapping[id].push_back({elm,n});
  }

  return myVxMapping;
}


void FFlLinkHandler::findWindowedNodes(std::map<int,FaVec3>& nodes,
                                       const std::vector<int>& indices,
                                       const FaMat34& lCS, bool lFirst,
                                       WindowTester isInsideWindow) const
{
  static std::set<int> visited;
  if (lFirst) visited.clear();

  int nVertices = myVertices.size();
  for (int idx : indices)
    if (idx < nVertices && visited.insert(idx).second)
    {
      const FFlVertex* vrtx = static_cast<const FFlVertex*>(myVertices[idx]);
      if (vrtx && vrtx->getNode())
      {
        FaVec3 globalPos = lCS * (*vrtx);
        if (isInsideWindow(globalPos))
          nodes[vrtx->getNode()->getID()] = globalPos;
      }
    }
}

#endif


void FFlLinkHandler::deleteResults()
{
  delete myResults;
  myResults = NULL;

  myFElements.clear();
}


void FFlLinkHandler::removeElements(const ElementsVec& toBeErased)
{
  for (FFlElementBase* elm : toBeErased)
  {
    myElements.erase(std::find(myElements.begin(),myElements.end(),elm));
    for (GroupMap::value_type& g : myGroupMap)
      g.second->removeElement(elm->getID());
    delete elm;
  }
}


bool FFlLinkHandler::resolve(bool subdivParabolic, bool fromSESAM)
{
  if (isResolved) return true;

#ifdef FFL_TIMER
  myProfiler->startTimer("resolve");
#endif
  int nError = 0;
  int nNotes = 0;
  const int maxErr = 100;

  // make sure everything is sorted
  if (!areElementsSorted) this->sortElements(true);
  if (!areNodesSorted)    this->sortNodes(true);
  if (!areLoadsSorted)    this->sortLoads();
#ifdef FT_USE_VISUALS
  if (!areVisualsSorted)  this->sortVisuals();
#endif

  ElementsVec oldElements;
  if (subdivParabolic)
    oldElements = myElements;
  else
  {
    // check for parabolic beam elements
    for (FFlElementBase* elm : myElements)
      if (elm->getTypeName() == "BEAM3")
        oldElements.push_back(elm);
    // always split the parabolic beams
    subdivParabolic = !oldElements.empty();
  }
  if (subdivParabolic)
  {
    // Subdivide all parabolic shells into four linear elements each
    // and parabolic beams into two linear elements
    ElementsVec oldElements(myElements);
    for (FFlElementBase* elm : oldElements)
      if (this->splitElement(elm) < 0)
        return false;
  }

  // resolve nodes:
  AttributeTypeCIter at = myAttributes.find("PCOORDSYS");
  if (at != myAttributes.end())
    for (FFlNode* node : myNodes)
      if (!node->resolveLocalSystem(at->second, nError >= maxErr))
        if (nError++ < maxErr)
          ListUI <<"\n *** Error: Resolving node "
                 << node->getID() <<" failed\n";

  // resolve finite elements:
  for (FFlElementBase* elm : myElements)
    if (!elm->resolveNodeRefs(myNodes, nError >= maxErr) ||
        !elm->resolveElmRef(myElements, nError >= maxErr) ||
        !elm->resolve(myAttributes, nError >= maxErr) ||
#ifdef FT_USE_VISUALS
        !elm->resolveVisuals(myVisuals, nError < maxErr)
#else
        false
#endif
        )
      if (nError++ < maxErr)
        ListUI <<"\n *** Error: Resolving "<< elm->getTypeName() <<" element "
               << elm->getID() <<" failed\n";

  // resolve loads:
  for (FFlLoadBase* load : myLoads)
    if (!load->resolveNodeRef(myNodes, nError >= maxErr) ||
        !load->resolveElmRef(myElements, nError >= maxErr) ||
        !load->resolve(myAttributes, nError >= maxErr))
      if (nError++ < maxErr)
        ListUI <<"\n *** Error: Resolving "<< load->getTypeName() <<" load "
               << load->getID() <<" failed\n";

  // resolve groups:
  for (const GroupMap::value_type& g : myGroupMap)
    if (!g.second->resolveElemRefs(myElements, nError >= maxErr))
      if (nError++ < maxErr)
        ListUI <<"\n *** Error: Resolving element group "
               << g.second->getID() <<" failed\n";

  // resolve attributes:
  for (const AttributeTypeMap::value_type& am : myAttributes)
    for (const AttributeMap::value_type& attr : am.second)
      if (!attr.second->resolve(myAttributes, nError >= maxErr))
        if (nError++ < maxErr)
          ListUI <<"\n *** Error: Resolving "<< am.first <<" attribute "
                 << attr.second->getID() <<" failed\n";

  // Remove loose nodes from WAVGM and RGD elements
  int nNod = 0;
  ElementsVec toBeErased;
  for (FFlElementBase* elm : myElements)
    if ((nNod = elm->getNodeCount()) < 1)
      toBeErased.push_back(elm);
    else if (elm->getTypeName() == "WAVGM")
    {
      NodeCIter refn = elm->nodesBegin();
      if (fromSESAM && !(*refn)->hasDOFs())
      {
        // Loose reference node ==> remove element
        toBeErased.push_back(elm);
        continue;
      }

      // Find all independent nodes that are not connected to any other elements
      NodeCIter n = refn;
      std::vector<int> looseNodes;
      for (++n; n != elm->nodesEnd(); ++n)
        if (!(*n)->hasDOFs())
        {
          if ((*n)->isExternal())
            (*n)->pushDOFs(6); // Ensure external node is not flagged DOF-less
          else
            looseNodes.push_back((*n)->getID());
        }

      if (nNod-looseNodes.size() < 2)
      {
        // Only the reference node is coupled to other elements.
        // This element will have no effect and is therefore removed.
        nNotes++;
        ListUI <<"\n   * Note: WAVGM element "<< elm->getID()
               <<" (size="<< nNod <<") has no connections";
        toBeErased.push_back(elm);
        continue;
      }

      (*refn)->pushDOFs(6); // Ensure the reference node is not flagged DOF-less

      if (looseNodes.empty()) continue;

      nNotes++;
      ListUI <<"\n  ** Removing "<< looseNodes.size()
             <<" loose nodes from WAVGM element "<< elm->getID();
      if (!static_cast<FFlWAVGM*>(elm)->removeMasterNodes(looseNodes))
        ++nError;
      else
      {
        FFlAttributeBase* newAtt = elm->getAttribute("PWAVGM");
        if (newAtt)
        {
          int ID = 1;
          while (this->getAttribute("PWAVGM",ID)) ID++;
          newAtt = static_cast<FFlPWAVGM*>(newAtt)->removeWeights(looseNodes,nNod);
          newAtt->setID(ID);
          elm->clearAttribute("PWAVGM"); // Replace the old attribute with the new one
          elm->setAttribute(this->getAttribute("PWAVGM",this->addUniqueAttribute(newAtt)));
        }
      }
    }
    else if (elm->getTypeName() == "RGD")
    {
      // Find all dependent nodes that are connected to other elements
      int looseDn = 0;
      NodeCIter n = elm->nodesBegin();
      std::vector<FFlNode*> depNodes;
      depNodes.reserve(nNod-1);
      for (++n; n != elm->nodesEnd(); ++n)
        if ((*n)->hasDOFs())
          depNodes.push_back(n->getReference());
        else
          looseDn++;

      if (depNodes.empty())
      {
        // Only the reference node is coupled to other elements.
        // This element will have no effect and is therefore removed.
        nNotes++;
        ListUI <<"\n   * Note: RGD element "<< elm->getID()
               <<" (size="<< nNod <<") has no connections";
        toBeErased.push_back(elm);
      }
      else if (looseDn > 0)
      {
        nNotes++;
        ListUI <<"\n  ** Removing "<< looseDn <<" loose nodes ("<< nNod-looseDn
               <<" remaining) from RGD element "<< elm->getID();
        elm->setNodes(depNodes,1,true);
      }
    }

  if (!toBeErased.empty())
  {
    // Erase all elements without nodal connections
    ListUI <<"\n  ** Erasing the ";
    if (toBeErased.size() == 1)
      ListUI << toBeErased.front()->getTypeName() <<" element";
    else
      ListUI << toBeErased.size() <<" constraint elements (RGD/WAVGM)";
    ListUI <<" without nodal connections.";
    this->removeElements(toBeErased);
  }

  // Check for unused nodes:
  size_t nnodes = 0;
  size_t unused = 0;
  for (size_t i = 0; i < myNodes.size(); i++)
    if (!myNodes[i]->getRefCount())
    {
      if (unused++ < 50)
      {
        nNotes++;
        ListUI <<"\n   * Note: Unused node "<< myNodes[i]->getID()
               <<" (deleted)";
      }
      delete myNodes[i];
    }
    else if (++nnodes < i+1)
      myNodes[nnodes-1] = myNodes[i];

  if (unused > 50)
    ListUI <<" ...\n  ** A total of "<< unused
           <<" unused nodes were deleted.\n";

  if (nnodes < myNodes.size())
    myNodes.resize(nnodes);

  // Check for unused attributes
  for (AttributeTypeMap::value_type& am : myAttributes)
    for (AttributeMap::iterator it = am.second.begin(); it != am.second.end();)
      if (!it->second->getRefCount())
      {
        nNotes++;
        ListUI <<"\n   * Note: Unused attribute "<< am.first <<" "
               << it->second->getID() <<" (deleted)";
        delete it->second;
        it = am.second.erase(it);
      }
      else
        ++it;

  if (nError > maxErr)
    ListUI <<"\n *** A total of "<< nError <<" resolve errors were detected.\n"
           <<"     (Only the first "<< maxErr <<" are reported.)\n";
  if (nNotes > 0)
    ListUI <<"\n  ** A total of "<< nNotes <<" resolve notes were detected.\n"
           <<"     Simulation continues, but you ought to verify that"
           <<" the model is consistent.\n";

  isResolved = nError == 0;
  int nTotNodes = myNodes.size();
  hasLooseNodes = isResolved && this->getNodeCount(FFL_FEM) < nTotNodes;

#ifdef FFL_TIMER
  myProfiler->stopTimer("resolve");
#endif
  return isResolved;
}


bool FFlLinkHandler::verify(bool fixNegElms)
{
  bool status = true;
  ElementsVec okElements;
  okElements.reserve(myElements.size());
  for (FFlElementBase* elm : myElements)
    switch (elm->checkOrientation(fixNegElms))
      {
      case -1:
        ListUI <<"  ** "<< elm->getTypeName() <<" element "<< elm->getID()
               <<" has negative volume\n";
        status = fixNegElms;
      case  1:
        okElements.push_back(elm);
        break;
      default:
        ListUI <<"  ** "<< elm->getTypeName() <<" element "<< elm->getID()
               <<" has zero volume (deleted)\n";
        delete elm;
      }

  if (okElements.size() < myElements.size())
    myElements.swap(okElements);

  return status;
}


bool FFlLinkHandler::getExtents(FaVec3& max, FaVec3& min) const
{
  bool ok = false;

  // Lambda function for updating the bounding box values
  auto&& updateExt = [&max,&min,&ok](const FaVec3& v)
  {
    if (!ok)
    {
      min = max = v;
      ok = true;
      return;
    }

    for (int i = 0; i < 3; i++)
      if (v[i] < min[i])
        min[i] = v[i];
      else if (v[i] > max[i])
        max[i] = v[i];
  };

#ifdef FT_USE_VERTEX
  for (FaVec3* vtx : myVertices)
    updateExt(*vtx);
#else
  for (FFlNode* node : myNodes)
    updateExt(node->getPos());
#endif

  return ok;
}


double FFlLinkHandler::getMeanElementSize() const
{
  double meanSize = 0.0;

  for (FFlElementBase* elm : myElements)
    meanSize += elm->getSize();

  if (!myElements.empty())
    meanSize /= myElements.size();

  return meanSize;
}


/*!
  Computes the mass properties of the FE model, i.e., the total mass \a M,
  the center of gravity \a Xcg and the inertia tensor \a I are computed.
*/

void FFlLinkHandler::getMassProperties(double& M, FaVec3& Xcg,
                                       FFaTensor3& I) const
{
  // First estimate the center of gravity based on the bounding box
  FaVec3 max, min;
  this->getExtents(max,min);
  FaVec3 BBcg = 0.5*(max+min);

  // Loop over all element integrating mass and inertias about BBcg
  double     Me  = M   = 0.0;
  FaVec3     Xec = Xcg = FaVec3();
  FFaTensor3 Ie  = I   = FFaTensor3();
  for (FFlElementBase* elm : myElements)
    if (elm->getMassProperties(Me,Xec,Ie))
    {
      Xec -= BBcg;
      Xcg += Xec*Me;
      M   += Me;
      I   += Ie.translateInertia(Xec,Me);
    }

  // Compute the real center of gravity and adjust the inertias accordingly
  double eps = 100.0*DBL_EPSILON;
  if (fabs(M) < eps) return;

  Xcg = BBcg + Xcg/M;
  I.translateInertia(BBcg-Xcg,-M);
}


#ifdef FT_USE_VERTEX

/*!
  Adds a vertex to the vertex container.
  The running ID of the vertex is also set.
*/

void FFlLinkHandler::addVertex(FFlVertex* aVertex)
{
  if (!aVertex) return;

  aVertex->ref();
  aVertex->setRunningID(myVertices.size());

  myVertices.push_back(aVertex);
}
#endif


/*!
  Count the number of elements of each type
*/

void FFlLinkHandler::countElements() const
{
  ElmTypeCount::iterator eit;
  for (FFlElementBase* elm : myElements)
    if ((eit = myNumElements.find(elm->getTypeName())) == myNumElements.end())
      myNumElements[elm->getTypeName()] = 1;
    else
      eit->second++;
}


void FFlLinkHandler::dump() const
{
  if (myNumElements.empty())
    this->countElements();

  std::cout <<"\nFFlLinkHandler::dump()"
            <<"\n   Elements:        "<< myElements.size();
  for (const ElmTypeCount::value_type& ec : myNumElements)
    std::cout <<"\n\t"<< ec.first <<"\t"<< ec.second;

  int nFEnod = this->getNodeCount(FFL_FEM);
  std::cout <<"\n   Nodes (free):    "<< myNodes.size()
            <<" ("<< myNodes.size() - nFEnod <<")";
#ifdef FT_USE_VERTEX
  std::cout <<"\n   Vertices:        "<< myVertices.size();
#endif
  std::cout <<"\n   Loads:           "<< myLoads.size();
  std::cout <<"\n   Groups:          "<< myGroupMap.size();
  std::cout <<"\n   Attribute types: "<< myAttributes.size();
  for (const AttributeTypeMap::value_type& am : myAttributes)
    std::cout <<"\n\t" << am.first <<":\t"<< am.second.size();

  std::cout <<"\n----"
            <<"\nFFlNode:          "<< sizeof(FFlNode)
#ifdef FT_USE_VERTEX
            <<"\nFFlVertex:        "<< sizeof(FFlVertex)
#endif
            <<"\nFFlElement(Base): "<< sizeof(FFlElementBase)
            <<"\nFFlLinkHander:    "<< sizeof(FFlLinkHandler)
            <<"\n----"<< std::endl;
}


/*!
  Sorting of elements and nodes in the order of increasing external IDs.
  A warning is issued if duplicated elements or nodes are found.
  If \a deleteDuplicates is true, the latter duplicated entry is deleted.
  \note Deleting duplicated entries is not safe after resolving has been done.
*/

int FFlLinkHandler::sortElementsAndNodes(bool deleteDuplicates) const
{
  int nDup = 0;
  if (!areElementsSorted) nDup += this->sortElements(deleteDuplicates);
  if (!areNodesSorted)    nDup += this->sortNodes(deleteDuplicates);
  return nDup;
}


int FFlLinkHandler::sortElements(bool deleteDuplicates) const
{
  std::sort(myElements.begin(),myElements.end(),FFlFEPartBaseLess());
  areElementsSorted = true;

  // Check that no elements share the same external ID
  int ndupElem = 0;
  for (size_t e = 1; e < myElements.size(); e++)
    if (myElements[e]->getID() == myElements[e-1]->getID())
    {
      ndupElem++;
      FFlElementBase* e1 = myElements[e-1];
      FFlElementBase* e2 = myElements[e];
      ListUI <<"\n  ** Warning: Two elements with the same ID  "
             << e1->getID() <<"  were found.";

      NodeCIter nit;
      ListUI <<"\n     First element:  "<< e1->getTypeName();
      for (nit = e1->nodesBegin(); nit != e1->nodesEnd(); nit++)
        ListUI <<" "<< nit->getID();
      ListUI <<"\n     Second element: "<< e2->getTypeName();
      for (nit = e2->nodesBegin(); nit != e2->nodesEnd(); nit++)
        ListUI <<" "<< nit->getID();
      ListUI <<"\n";

      // This matters only if the elements are of the same cathegory
      if (e1->getCathegory() == e2->getCathegory())
      {
        ListUI <<"     This may result in inconsistency in the FE topology.\n";
        if (deleteDuplicates)
        {
          ListUI <<"     The latter element is therefore deleted.\n";
          ElementsIter eit = myElements.begin(); eit += e;
          delete *eit;
          myElements.erase(eit);
          myNumElements.clear();
          e--;
        }
      }
    }

  return ndupElem; // Number of duplicated elements found
}


int FFlLinkHandler::sortNodes(bool deleteDuplicates) const
{
  std::sort(myNodes.begin(),myNodes.end(),FFlFEPartBaseLess());
  areNodesSorted = true;

  // Check that no nodes share the same external ID
  int ndupNode = 0;
  for (size_t n = 1; n < myNodes.size(); n++)
    if (myNodes[n]->getID() == myNodes[n-1]->getID())
    {
      ndupNode++;
      ListUI <<"\n  ** Warning: Two nodes with the same ID  "
             << myNodes[n]->getID() <<"  were found."
             <<"\n     First node:  "<< myNodes[n-1]->getPos()
             <<"\n     Second node: "<< myNodes[n]->getPos() <<"\n";

      if (deleteDuplicates)
      {
        ListUI <<"     The latter node is deleted.\n";
        NodesIter nit = myNodes.begin(); nit += n;
        delete *nit;
        myNodes.erase(nit);
        n--;
      }
    }

  return ndupNode; // Number of duplicated nodes found
}


void FFlLinkHandler::sortLoads() const
{
  std::stable_sort(myLoads.begin(),myLoads.end(),FFlFEPartBaseLess());
  areLoadsSorted = true;
}


#ifdef FT_USE_VISUALS
void FFlLinkHandler::sortVisuals() const
{
  std::sort(myVisuals.begin(),myVisuals.end(),FFlFEPartBaseLess());
  areVisualsSorted = true;
}
#endif


#ifdef FT_USE_CONNECTORS

/*!
  Creates a connector based on input geometry and the nodal position.
  The elements and nodes are returned through \a cItems.
  The \a spiderType is either:
   2 - RiGiD element (RBE2)
   3 - Weighted AVerage Motion element (RBE3)
*/

int FFlLinkHandler::createConnector(const FFaCompoundGeometry& compound,
                                    const FaVec3& nodePos, int spiderType,
                                    FFlConnectorItems& cItems)
{
  if (spiderType < 2 || spiderType > 3)
    return -1; // Invalid spider type

  std::vector<FFlNode*> nodes;
  for (FFlNode* node : myNodes)
    if (node->hasDOFs() && node->getStatus() == FFlNode::INTERNAL)
      if (compound.isInside(node->getPos()))
        nodes.push_back(node);

  if (nodes.size() < 2)
    return 0; // No FE nodes on the provided geometry

  cItems.clear();

  FFlNode* refNode = new FFlNode(this->getNewNodeID(),nodePos);
  refNode->pushDOFs(6);
  this->addNode(refNode,true);
  cItems.addNode(refNode->getID());
  ListUI <<"  -> Creating FE node "<< refNode->getID() <<"\n";
  nodes.insert(nodes.begin(),refNode);

  const char* spider[2] = { "RGD", "WAVGM" };
  FFlElementBase* newEl = ElementFactory::instance()->create(spider[spiderType-2],
                                                             this->getNewElmID());
  newEl->setNodes(nodes);
  this->addElement(newEl,true);
  cItems.addElement(newEl->getID());
  ListUI <<"  -> Creating "<< (spiderType == 2 ? "RGD":"WAVGM")
         <<" element "<< newEl->getID()
         <<" with reference node "<< refNode->getID() <<"\n";

  if (spiderType == 3)
    // Create an attachable node co-located with the reference node
    refNode = this->createAttachableNode(refNode,refNode->getPos(),&cItems);

  // Ensure the attachable node has external status
  refNode->setExternal();

  return nodes.size();
}


/*!
  Deletes the elements and nodes in the connector.
*/

int FFlLinkHandler::deleteConnector(const FFlConnectorItems& cItems)
{
  int nDeleted = 0;
  ElementsIter eit;
  NodesIter nit;

  for (int elm : cItems.getElements())
    if ((eit = this->getElementIter(elm)) != myElements.end())
    {
      ListUI <<"  -> Deleting "<< (*eit)->getTypeName()
             <<" element "<< elm <<"\n";
      if ((*eit)->getTypeName() == "BUSH")
        myBushElements.erase(*eit);
      delete *eit;
      myElements.erase(eit);
      nDeleted++;
    }

  for (int node : cItems.getNodes())
    if ((nit = this->getNodeIter(node)) != myNodes.end())
    {
      ListUI <<"  -> Deleting FE node "<< node <<"\n";
      delete *nit;
      myNodes.erase(nit);
      nDeleted++;
    }

  this->sortElements();
  this->sortNodes();

  return nDeleted;
}
#endif


/*!
  This method will subdivide a 3-noded beam into two linear beam elements, and
  a 6-noded triangle and a 8-noded quadrilateral into four linear equivalents.
*/

int FFlLinkHandler::splitElement(FFlElementBase* elm)
{
  int nod_id = this->getNewNodeID();
  if (elm->getTypeName() == "QUAD8")
  {
    // Need a new node at the center of the 8-noded element
    FFlNode* newNode = NULL;
    FaVec3 X;
    for (int i = 1; i <= 8; i++)
      if ((newNode = this->getNode(elm->getNodeID(i))))
        X += ((i % 2) ? -0.25 : 0.5)*newNode->getPos();
    if (!this->addNode(newNode = new FFlNode(nod_id,X)))
      return -5;
  }

  ElementsVec newElems;
  if (!elm->split(newElems,this,nod_id))
    return 0;

  std::vector<int> newElmID;
  newElmID.reserve(newElems.size());
  for (FFlElementBase* newElm : newElems)
  {
    newElmID.push_back(newElm->getID());
    newElm->useAttributesFrom(elm);
#ifdef FT_USE_VISUALS
    newElm->useVisualsFrom(elm);
#endif
    if (!this->addElement(newElm))
      return -newElmID.size();
  }

  // Erase the parabolic element
  int oldElmID = elm->getID();
  myElements.erase(std::find(myElements.begin(),myElements.end(),elm));
  delete elm;

  // Update element groups
  for (GroupMap::value_type& g : myGroupMap)
    g.second->swapElement(oldElmID,newElmID);

  return newElems.size();
}

#ifdef FF_NAMESPACE
} // namespace
#endif
