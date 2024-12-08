// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlLoads.H"
#include "FFlLib/FFlFEParts/FFlPORIENT.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


FFlCLOAD::FFlCLOAD(int id, const char type) : FFlLoadBase(id), myType(type)
{
  this->addField(P);

  target_counter = 0;
}


FFlCLOAD::FFlCLOAD(const FFlCLOAD& obj) : FFlLoadBase(obj), myType(obj.myType)
{
  this->addField(P);

  P = obj.P;
  target_counter = 0;
  target.resize(obj.target.size());
  for (size_t i = 0; i < target.size(); i++)
    target[i] = obj.target[i].getID();
}


bool FFlCLOAD::resolveNodeRef(const std::vector<FFlNode*>& possibleNodes,
			      bool suppressErrmsg)
{
  if (possibleNodes.empty())
    return false;

  bool retVar = true;
  for (size_t i = 0; i < target.size(); i++)
    if (!target[i].resolve(possibleNodes))
    {
      if (!suppressErrmsg)
	ListUI <<"\n *** Error: Failed to resolve reference to node "
	       << target[i].getID() <<"\n";
      retVar = false;
    }

  return retVar;
}


void FFlCLOAD::setTarget(int nodeId, int)
{
  target.resize(target.size()+1,nodeId);
}


void FFlCLOAD::setTarget(const std::vector<int>& nodIds)
{
  size_t n = target.size();
  target.resize(n+nodIds.size());
  for (size_t i = 0; i < nodIds.size(); i++)
    target[n+i] = nodIds[i];
}


bool FFlCLOAD::getTarget(int& nodID, int& type) const
{
  if (type == 0) target_counter = 0;

  if (target_counter >= target.size())
  {
    target_counter = 0;
    return false;
  }

  nodID = target[target_counter++].getID();
  type = myType == 'F' ? -1 : -2;
  return true;
}


int FFlCLOAD::getLoad(std::vector<FaVec3>& F, int& type) const
{
  if (type == 0) target_counter = 0;

  if (target_counter >= target.size())
  {
    target_counter = 0;
    return 0;
  }

  F.clear();
  F.push_back(P.getValue());
  type = myType == 'F' ? -1 : -2;

  return target[target_counter++].getID();
}


void FFlCLOAD::calculateChecksum(FFaCheckSum* cs, int csMask) const
{
  FFlLoadBase::calculateChecksum(cs,csMask);

  cs->add(myType);
  for (size_t i = 0; i < target.size(); i++)
    cs->add(target[i].getID());
}


void FFlCFORCE::convertUnits(const FFaUnitCalculator* convCal)
{
  convCal->convert(P.data(),"FORCE");
}


void FFlCFORCE::init()
{
  using TypeInfoSpec = FFaSingelton<FFlTypeInfoSpec,FFlCFORCE>;

  TypeInfoSpec::instance()->setTypeName("CFORCE");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::LOAD);

  LoadFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlCFORCE::create,int,FFlLoadBase*&));
}


void FFlCMOMENT::convertUnits(const FFaUnitCalculator* convCal)
{
  convCal->convert(P.data(),"FORCE/LENGTH");
}


void FFlCMOMENT::init()
{
  using TypeInfoSpec = FFaSingelton<FFlTypeInfoSpec,FFlCMOMENT>;

  TypeInfoSpec::instance()->setTypeName("CMOMENT");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::LOAD);

  LoadFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlCMOMENT::create,int,FFlLoadBase*&));
}


FFlPLOAD::FFlPLOAD(int id) : FFlLoadBase(id)
{
  this->addField(P);
}


FFlPLOAD::FFlPLOAD(const FFlPLOAD& obj) : FFlLoadBase(obj)
{
  this->addField(P);
  P = obj.P;
}


void FFlPLOAD::convertUnits(const FFaUnitCalculator* convCal)
{
  std::vector<double>& p = P.data();
  for (size_t i = 0; i < p.size(); i++)
    convCal->convert(p[i],"FORCE/AREA",10);
}


/*!
  Returns the nodal face load intensities in global coordinates.
  One PLOAD object may act on several element faces, but the load on one face
  only is returned in each call. To get the first element/face, \a face must be
  set equal to zero on entry. The external ID of the element the load is acting
  on is returned, and \a face is updated to indicate which face for solids.
  When all load targets have been processed, zero is returned.
*/

int FFlPLOAD::getLoad(std::vector<FaVec3>& Pglb, int& face) const
{
  const std::vector<double>& p = P.getValue();
  if (p.empty()) return 0;

  // Get the next element and local face which this load acts on
  const FFlElementBase* elm = this->getTarget(face);
  if (!elm) return 0;

  // Determine the load intensity at the face nodes in local coordinate system
  size_t nFaceNodes = face > 0 ? elm->getFaceSize(face) : elm->getNodeCount();
  std::vector<double> Ploc(nFaceNodes,p.front());
  if (p.size() > 1)
  {
    if (nFaceNodes <= p.size() || nFaceNodes <= 4)
      for (size_t i = 1; i < nFaceNodes; i++)
	Ploc[i] = p[i < p.size() ? i : 0];
    else if (nFaceNodes <= 2*p.size()) // 2nd order face, interpolate mid-points
      for (size_t i = 1; i < nFaceNodes; i++)
	Ploc[i] = i%2 ? 0.5*(p[i/2]+p[i/2+1 < p.size() ? i/2+1 : 0]) : p[i/2];
    else
      return 0;
  }

  // Compute load intensities in global coordinate system
  Pglb.clear();
  Pglb.reserve(nFaceNodes);
  FFlPORIENT* ori = dynamic_cast<FFlPORIENT*>(this->getAttribute("PORIENT"));
  if (ori) // Load acts in a specified global direction
    for (size_t i = 0; i < nFaceNodes; i++)
      Pglb.push_back(ori->directionVector.getValue()*Ploc[i]);
  else if (elm->getFaceNormals(Pglb,face)) // Load acts normal to face
    for (size_t i = 0; i < nFaceNodes; i++)
      Pglb[i] *= (face > 0 ? -Ploc[i] : Ploc[i]);
  else
    return 0;

  if (face < 0) face = 1;
  return elm->getID();
}


FFlSURFLOAD::FFlSURFLOAD(const FFlSURFLOAD& obj) : FFlPLOAD(obj)
{
  target_counter = 0;
  target.resize(obj.target.size());
  for (size_t i = 0; i < target.size(); i++)
    target[i] = obj.target[i].getID();
}


bool FFlSURFLOAD::resolveElmRef(const std::vector<FFlElementBase*>& possibleElms,
				bool suppressErrmsg)
{
  if (possibleElms.empty())
    return false;

  bool retVar = true;
  for (size_t i = 0; i < target.size(); i++)
    if (!target[i].resolve(possibleElms))
    {
      if (!suppressErrmsg)
	ListUI <<"\n *** Error: Failed to resolve reference to element "
	       << target[i].getID() <<"\n";
      retVar = false;
    }
    else if (target[i]->getCathegory() != FFlTypeInfoSpec::SHELL_ELM)
    {
      if (!suppressErrmsg)
	ListUI <<"\n *** Error: Surface load is referring to non-shell element "
	       << target[i]->getTypeName() <<" "<< target[i].getID() <<"\n";
      retVar = false;
    }

  return retVar;
}


void FFlSURFLOAD::setTarget(const std::vector<int>& elmIds)
{
  size_t n = target.size();
  target.resize(n+elmIds.size());
  for (size_t i = 0; i < elmIds.size(); i++)
    target[n+i] = elmIds[i];
}


bool FFlSURFLOAD::getTarget(int& elmID, int& face) const
{
  if (face == 0) target_counter = 0;

  if (target_counter >= target.size())
  {
    target_counter = 0;
    return false;
  }

  elmID = target[target_counter++].getID();
  face = -1;
  return true;
}


const FFlElementBase* FFlSURFLOAD::getTarget(int& face) const
{
  if (face == 0) target_counter = 0;

  if (target_counter >= target.size())
  {
    target_counter = 0;
    return 0;
  }

  face = -1;
  if (target[target_counter].isResolved())
    return target[target_counter++].getReference();

  target_counter++;
  return 0;
}


void FFlSURFLOAD::calculateChecksum(FFaCheckSum* cs, int csMask) const
{
  FFlPLOAD::calculateChecksum(cs,csMask);

  for (size_t i = 0; i < target.size(); i++)
    cs->add(target[i].getID());
}


void FFlSURFLOAD::init()
{
  using TypeInfoSpec  = FFaSingelton<FFlTypeInfoSpec,FFlSURFLOAD>;
  using AttributeSpec = FFaSingelton<FFlFEAttributeSpec,FFlSURFLOAD>;

  TypeInfoSpec::instance()->setTypeName("SURFLOAD");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::LOAD);
  AttributeSpec::instance()->addLegalAttribute("PORIENT",false);

  LoadFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlSURFLOAD::create,int,FFlLoadBase*&));
}


FFlFACELOAD::FFlFACELOAD(const FFlFACELOAD& obj) : FFlPLOAD(obj)
{
  target_counter = 0;
  target.resize(obj.target.size());
  for (size_t i = 0; i < target.size(); i++)
  {
    target[i].first  = obj.target[i].first.getID();
    target[i].second = obj.target[i].second;
  }
}


bool FFlFACELOAD::resolveElmRef(const std::vector<FFlElementBase*>& possibleElms,
				bool suppressErrmsg)
{
  if (possibleElms.empty())
    return false;

  bool retVar = true;
  for (size_t i = 0; i < target.size(); i++)
    if (!target[i].first.resolve(possibleElms))
    {
      if (!suppressErrmsg)
	ListUI <<"\n *** Error: Failed to resolve reference to element "
	       << target[i].first.getID() <<"\n";
      retVar = false;
    }
    else if (target[i].first->getCathegory() != FFlTypeInfoSpec::SOLID_ELM)
    {
      if (!suppressErrmsg)
	ListUI <<"\n *** Error: Face load is referring to non-solid element "
	       << target[i].first->getTypeName() <<" "
	       << target[i].first->getID() <<"\n";
      retVar = false;
    }

  return retVar;
}


void FFlFACELOAD::setTarget(int elmID, int face)
{
  for (size_t i = 0; i < target.size(); i++)
    if (target[i].first.getID() == elmID)
      target[i].second = face;
}


void FFlFACELOAD::setTarget(const std::vector<int>& elmIds)
{
  size_t n = target.size();
  target.resize(n+(elmIds.size()+1)/2);
  for (size_t i = 0; i < elmIds.size(); i++)
    if (i%2 == 0)
      target[n+i/2].first = elmIds[i];
    else
      target[n+i/2].second = elmIds[i];
}


bool FFlFACELOAD::getTarget(int& elmID, int& face) const
{
  if (face == 0) target_counter = 0;

  if (target_counter >= target.size())
  {
    target_counter = 0;
    return false;
  }

  elmID = target[target_counter].first.getID();
  face = target[target_counter++].second;
  return true;
}


const FFlElementBase* FFlFACELOAD::getTarget(int& face) const
{
  if (face == 0) target_counter = 0;

  if (target_counter >= target.size())
  {
    target_counter = 0;
    return 0;
  }

  face = target[target_counter].second;
  if (target[target_counter].first.isResolved())
    return target[target_counter++].first.getReference();

  target_counter++;
  return 0;
}


void FFlFACELOAD::calculateChecksum(FFaCheckSum* cs, int csMask) const
{
  FFlPLOAD::calculateChecksum(cs,csMask);

  for (size_t i = 0; i < target.size(); i++)
  {
    cs->add(target[i].first.getID());
    cs->add(target[i].second);
  }
}


void FFlFACELOAD::init()
{
  using TypeInfoSpec  = FFaSingelton<FFlTypeInfoSpec,FFlFACELOAD>;
  using AttributeSpec = FFaSingelton<FFlFEAttributeSpec,FFlFACELOAD>;

  TypeInfoSpec::instance()->setTypeName("FACELOAD");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::LOAD);
  AttributeSpec::instance()->addLegalAttribute("PORIENT",false);

  LoadFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlFACELOAD::create,int,FFlLoadBase*&));
}

#ifdef FF_NAMESPACE
} // namespace
#endif
