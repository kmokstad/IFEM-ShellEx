// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <algorithm>

#include "FFlLib/FFlFEAttributeRefs.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFaLib/FFaAlgebra/FFaCheckSum.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


FFlFEAttributeRefs::FFlFEAttributeRefs(const FFlFEAttributeRefs& obj)
{
  myAttributes.reserve(obj.myAttributes.size());
  for (const AttribData& aref : obj.myAttributes)
    myAttributes.push_back(AttribData(aref.first,AttribRef(aref.second->getID())));
}


bool FFlFEAttributeRefs::useAttributesFrom(const FFlFEAttributeRefs* obj)
{
  bool ok = true;
  myAttributes.reserve(obj->myAttributes.size());
  for (const AttribData& aref : obj->myAttributes)
  {
    const std::string& typeName = obj->getAttributeName(aref.first);
    if (this->getAttributeTypeID(typeName))
      ok &= this->setAttribute(typeName,aref.second.getID());
  }

  return ok;
}


bool FFlFEAttributeRefs::setAttribute(FFlAttributeBase* attrObject)
{
  const std::string& type = attrObject->getTypeInfoSpec()->getTypeName();
  unsigned char typeID = this->getAttributeTypeID(type);
  if (!typeID)
  {
    std::cerr <<" *** Internal error: \""<< type
              <<"\" is not a legal attribute."<< std::endl;
    return false;
  }

  std::pair<AttribsCIter,AttribsCIter> p;
  p = std::equal_range(myAttributes.begin(), myAttributes.end(),
                       typeID, FFlFEAttribTypeLess());
  if (p.first == p.second ||
      this->getFEAttributeSpec()->multipleRefsAllowed(type))
  {
    myAttributes.push_back(AttribData(typeID,AttribRef(attrObject)));
    std::sort(myAttributes.begin(), myAttributes.end(), FFlFEAttribTypeLess());
  }
#ifdef FFL_DEBUG
  else
    ListUI <<"\n  ** Warning: Attribute reference \""<< type
           <<"\" is already set to "<< p.first->second.getID()
           <<", "<< attrObject->getID() <<" is ignored.\n";
#endif

  return true;
}


bool FFlFEAttributeRefs::setAttribute(const std::string& type, int ID)
{
  unsigned char typeID = this->getAttributeTypeID(type);
  if (!typeID)
  {
    if (FFlFEAttributeSpec::isObsolete(type))
      return true; // silently ignore the obsolete attributes

    std::cerr <<" *** Internal error: \""<< type
              <<"\" is not a legal attribute."<< std::endl;
    return false;
  }

  std::pair<AttribsCIter,AttribsCIter> p;
  p = std::equal_range(myAttributes.begin(), myAttributes.end(),
                       typeID, FFlFEAttribTypeLess());
  if (p.first == p.second ||
      this->getFEAttributeSpec()->multipleRefsAllowed(type))
  {
    myAttributes.push_back(AttribData(typeID,AttribRef(ID)));
    std::sort(myAttributes.begin(), myAttributes.end(), FFlFEAttribTypeLess());
  }
#ifdef FFL_DEBUG
  else
    ListUI <<"\n  ** Warning: Attribute reference \""<< type
           <<"\" is already set to "<< p.first->second.getID()
           <<", "<< ID <<" is ignored.\n";
#endif

  return true;
}


bool FFlFEAttributeRefs::clearAttribute(const std::string& type)
{
  unsigned char typeID = this->getAttributeTypeID(type);
  if (typeID)
  {
    std::pair<AttribsVec::iterator,AttribsVec::iterator> p;
    p = std::equal_range(myAttributes.begin(), myAttributes.end(),
                         typeID, FFlFEAttribTypeLess());
    if (p.first != p.second)
      myAttributes.erase(p.first,p.second);
    return true;
  }
  else if (!FFlFEAttributeSpec::isObsolete(type))
    std::cerr <<" *** Internal error: \""<< type
              <<"\" is not a legal attribute."<< std::endl;

  return false;
}


bool FFlFEAttributeRefs::hasAttribute(const FFlAttributeBase* attrib) const
{
  for (const AttribData& aref : myAttributes)
    if (aref.second.isResolved())
    {
      if (aref.second.getReference() == attrib)
        return true;
      else if (aref.second->hasAttribute(attrib))
        return true;
    }

  return false;
}


bool FFlFEAttributeRefs::hasAttribute(const std::vector<FFlAttributeBase*>& av) const
{
  for (FFlAttributeBase* attr : av)
    for (const AttribData& aref : myAttributes)
      if (aref.second.isResolved())
        if (aref.second.getReference() == attr)
          return true;

  // not found in first pass, running recursively:
  for (const AttribData& aref : myAttributes)
    if (aref.second.isResolved())
      if (aref.second->hasAttribute(av))
        return true;

  return false;
}


bool FFlFEAttributeRefs::resolve(const AttribTypMap& possibleRefs,
                                 bool suppressErrmsg)
{
  if (possibleRefs.empty() && !myAttributes.empty())
  {
    ListUI <<"\n *** Error: No attributes!\n";
    return false;
  }

  bool retVar = true;

  for (AttribData& attr : myAttributes)
  {
    // find correct sub-map in the input argument
    const std::string& attrName = this->getAttributeName(attr.first);
    AttribTypMap::const_iterator refIt = possibleRefs.find(attrName);

    // Workaround for conversion of obsolete field names into new name
    if (refIt == possibleRefs.end())
      if (attrName == "PBEAMORIENT" || attrName == "PBUSHORIENT")
        refIt = possibleRefs.find("PORIENT");

    if (refIt != possibleRefs.end())
      if (attr.second.resolve(refIt->second))
        continue;

    if (!suppressErrmsg)
      ListUI <<"\n *** Error: Failed to resolve "<< attrName
             <<" attribute "<< attr.second.getID() <<"\n";
    retVar = false;
  }

  return retVar;
}


void FFlFEAttributeRefs::checksum(FFaCheckSum* cs) const
{
  for (const AttribData& attr : myAttributes)
    cs->add(attr.second.getID());
}


FFlAttributeBase* FFlFEAttributeRefs::getAttribute(const std::string& atType) const
{
  unsigned char typeID = this->getAttributeTypeID(atType);
  if (typeID)
  {
    std::pair<AttribsCIter,AttribsCIter> p;
    p = std::equal_range(myAttributes.begin(), myAttributes.end(),
                         typeID, FFlFEAttribTypeLess());
    if (p.first != p.second && p.first->second.isResolved())
      return p.first->second.getReference();
  }
  else if (!FFlFEAttributeSpec::isObsolete(atType))
    std::cerr <<" *** Internal error: \""<< atType
              <<"\" is not a legal attribute."<< std::endl;

  return NULL;
}


std::vector<FFlAttributeBase*>
FFlFEAttributeRefs::getAttributes(const std::string& atType) const
{
  std::vector<FFlAttributeBase*> retVar;

  unsigned char typeID = this->getAttributeTypeID(atType);
  if (typeID)
  {
    std::pair<AttribsCIter,AttribsCIter> p;
    p = std::equal_range(myAttributes.begin(), myAttributes.end(),
                         typeID, FFlFEAttribTypeLess());
    for (AttribsCIter it = p.first; it < p.second; ++it)
      if (it->second.isResolved())
        retVar.push_back(it->second.getReference());
  }
  else if (!FFlFEAttributeSpec::isObsolete(atType))
    std::cerr <<" *** Internal error: \""<< atType
              <<"\" is not a legal attribute."<< std::endl;

  return retVar;
}


int FFlFEAttributeRefs::getAttributeID(const std::string& atType) const
{
  unsigned char typeID = this->getAttributeTypeID(atType);
  if (typeID)
  {
    std::pair<AttribsCIter,AttribsCIter> p;
    p = std::equal_range(myAttributes.begin(), myAttributes.end(),
                         typeID, FFlFEAttribTypeLess());
    if (p.first != p.second)
      return p.first->second.getID();
  }

  // No message if no attribute of given type or illegal attribute type
  return 0;
}


const std::string& FFlFEAttributeRefs::getAttributeName(unsigned char typeID) const
{
  FFlFEAttributeSpec* spec = this->getFEAttributeSpec();
  if (spec) return spec->getAttributeName(typeID);

  static std::string empty;
  return empty;
}


unsigned char FFlFEAttributeRefs::getAttributeTypeID(const std::string& name) const
{
  FFlFEAttributeSpec* spec = this->getFEAttributeSpec();
  return spec ? spec->getAttributeTypeID(name) : 0;
}

#ifdef FF_NAMESPACE
} // namespace
#endif
