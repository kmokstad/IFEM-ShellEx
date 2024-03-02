// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlCLOAD.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"


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
  FFlCFORCETypeInfoSpec::instance()->setTypeName("CFORCE");
  FFlCFORCETypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::LOAD);

  LoadFactory::instance()->registerCreator
    (FFlCFORCETypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlCFORCE::create,int,FFlLoadBase*&));
}


void FFlCMOMENT::convertUnits(const FFaUnitCalculator* convCal)
{
  convCal->convert(P.data(),"FORCE/LENGTH");
}


void FFlCMOMENT::init()
{
  FFlCMOMENTTypeInfoSpec::instance()->setTypeName("CMOMENT");
  FFlCMOMENTTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::LOAD);

  LoadFactory::instance()->registerCreator
    (FFlCMOMENTTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlCMOMENT::create,int,FFlLoadBase*&));
}
