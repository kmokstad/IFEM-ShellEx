// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlFEParts/FFlPCOORDSYS.H"
#include "FFlLib/FFlFEResultBase.H"
#ifdef FT_USE_VERTEX
#include "FFlLib/FFlVertex.H"
#endif
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"
#include "FFaLib/FFaAlgebra/FFaCheckSum.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"


/*!
  \class FFlNode FFlNode.H
  \brief FE node class. Stores node position and additional data for the node.
*/

FFlNode::FFlNode(int id) : FFlPartBase(id)
{
  status = 0; // Default is internal node
  myDOFCount = 0;

#ifdef FT_USE_VERTEX
  myVertex = NULL;
#endif
  myResults = NULL;
}


FFlNode::FFlNode(int id, double x, double y, double z, int s) : FFlPartBase(id)
{
  status = s;
  myDOFCount = 0;

#ifdef FT_USE_VERTEX
  myVertex = new FFlVertex(x,y,z);
  myVertex->ref();
  myVertex->setNode(this);
#else
  myPos = FaVec3(x,y,z);
#endif

  myResults = NULL;
}


FFlNode::FFlNode(int id, const FaVec3& pos, int s) : FFlPartBase(id)
{
  status = s;
  myDOFCount = 0;

#ifdef FT_USE_VERTEX
  myVertex = new FFlVertex(pos);
  myVertex->ref();
  myVertex->setNode(this);
#else
  myPos = pos;
#endif

  myResults = NULL;
}


FFlNode::FFlNode(const FFlNode& otherNode) : FFlPartBase(otherNode)
{
  status = otherNode.status;
  myDOFCount = otherNode.myDOFCount;

#ifdef FT_USE_VERTEX
  if (otherNode.myVertex)
  {
    myVertex = new FFlVertex(*otherNode.myVertex);
    myVertex->ref();
    myVertex->setNode(this);
  }
  else
    myVertex = NULL;
#else
  myPos = otherNode.myPos;
#endif

  myResults = NULL;
}


FFlNode::~FFlNode()
{
#ifdef FT_USE_VERTEX
  if (myVertex)
  {
    myVertex->setNode(NULL);
    myVertex->unRef();
  }
#endif
  this->deleteResults();
}


void FFlNode::init()
{
  FFlNodeTypeInfoSpec::instance()->setTypeName("Node");
  FFlNodeTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::NODE);
}


void FFlNode::calculateChecksum(FFaCheckSum* cs, int precision,
                                bool includeExtNodeInfo)
{
  FFlPartBase::checksum(cs);

#ifdef FT_USE_VERTEX
  if (myVertex)
    cs->add(*myVertex,precision);
#else
  cs->add(myPos,precision);
#endif

  cs->add(includeExtNodeInfo && status == 1);
  if (status < 0) cs->add((int)status);

  int localCS = myLocalSystem.getID();
  if (localCS > 0) cs->add(localCS);
}


void FFlNode::convertUnits(const FFaUnitCalculator* convCal)
{
#ifdef FT_USE_VERTEX
  if (myVertex) convCal->convert(*myVertex, "LENGTH");
#else
  convCal->convert(myPos, "LENGTH");
#endif
}


#ifdef FT_USE_VERTEX
int FFlNode::getVertexID() const
{
  if (myVertex) return myVertex->getRunningID();

  std::cerr <<"FFlNode::getVertexID(): No vertex set up for node "
            << this->getID() << std::endl;

  return -1;
}


const FaVec3& FFlNode::getPos() const
{
  if (myVertex) return *myVertex;

  std::cerr <<"FFlNode::getPos(): No vertex set up for node "
            << this->getID() << std::endl;

  static FaVec3 dummy;
  return dummy;
}


void FFlNode::setVertex(FFlVertex* aVertex)
{
  if (myVertex)
  {
    myVertex->setNode(NULL);
    myVertex->unRef();
  }
  myVertex = aVertex;
  myVertex->ref();
  myVertex->setNode(this);
}
#endif


bool FFlNode::setStatus(int newStat)
{
  if (status == newStat) return false;

  status = newStat;
  return true;
}


bool FFlNode::setExternal(bool ext)
{
  if (status > 1 || status == (int)ext) return false;

  status = ext;
  return true;
}


bool FFlNode::isFixed(int dof) const
{
  if (status < 0)
    switch (dof) {
    case 1: return -status & 1;
    case 2: return -status & 2;
    case 3: return -status & 4;
    case 4: return -status & 8;
    case 5: return -status & 16;
    case 6: return -status & 32;
    default: return true;
    }

  return false;
}


void FFlNode::setLocalSystem(const FFlPCOORDSYS* coordSys)
{
  myLocalSystem = coordSys;
}


bool FFlNode::resolveLocalSystem(const std::map<int,FFlAttributeBase*>& possibleCSs,
                                 bool suppressErrmsg)
{
  if (myLocalSystem.resolve(possibleCSs))
  {
    FFlAttributeBase* localCS = myLocalSystem.getReference();
    if (!localCS)
      return true; // no local coordinate system for this node
    else if (localCS->getTypeName() == "PCOORDSYS")
      return true; // we have local system

    // Logic error, indicates programming error
    myLocalSystem = (FFlAttributeBase*)NULL;
    std::cerr <<"FFlNode::resolveSystem(): Invalid attribute type provided "
              << localCS->getTypeName() << std::endl;
  }
  else if (!suppressErrmsg)
    ListUI <<"\n *** Error: Failed to resolve PCOORDSYS "
           << myLocalSystem.getID() <<"\n";

  return false;
}


int FFlNode::getLocalSystemID() const
{
  return myLocalSystem.getID();
}


FFlPCOORDSYS* FFlNode::getLocalSystem() const
{
  if (myLocalSystem.isResolved())
    return static_cast<FFlPCOORDSYS*>(myLocalSystem.getReference());
  else
    return NULL;
}


void FFlNode::deleteResults()
{
  delete myResults;
  myResults = NULL;
}
