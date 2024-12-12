// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_FE_ATTRIBUTE_REFS_H
#define FFL_FE_ATTRIBUTE_REFS_H

#include "FFlLib/FFlReference.H"
#include <vector>
#include <string>
#include <map>

class FFaCheckSum;
class FFlFEAttributeSpec;

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif

class FFlAttributeBase;

using AttribRef    = FFlReference<FFlAttributeBase>;
using AttribData   = std::pair<unsigned char,AttribRef>;
using AttribsVec   = std::vector<AttribData>;
using AttribsCIter = AttribsVec::const_iterator;
using AttributeMap = std::map<int,FFlAttributeBase*>;
using AttribTypMap = std::map<std::string,AttributeMap>;


class FFlFEAttributeRefs
{
protected:
  FFlFEAttributeRefs() {}
  FFlFEAttributeRefs(const FFlFEAttributeRefs& obj);
  FFlFEAttributeRefs& operator=(const FFlFEAttributeRefs&) = delete;

public:
  virtual ~FFlFEAttributeRefs() {}

  bool useAttributesFrom(const FFlFEAttributeRefs* obj);

  // Attributes are always sorted when inserted. If this becomes a performance
  // problem, we should see what to do with that without increasing object size.
  // We usually have only a few attributes, so it should not be a large problem.
  bool setAttribute(FFlAttributeBase* attribute);
  bool setAttribute(const std::string& type, int ID);
  bool clearAttribute(const std::string& type);

  bool resolve(const AttribTypMap& possibleRefs, bool suppressErrmsg = false);

  virtual FFlFEAttributeSpec* getFEAttributeSpec() const { return NULL; }

  AttribsCIter attributesBegin() const { return myAttributes.begin(); }
  AttribsCIter attributesEnd()   const { return myAttributes.end(); }

  FFlAttributeBase*              getAttribute (const std::string& type) const;
  std::vector<FFlAttributeBase*> getAttributes(const std::string& type) const;

  int getAttributeID(const std::string& type) const;

  bool hasAttribute(const FFlAttributeBase* attribute) const;
  bool hasAttribute(const std::vector<FFlAttributeBase*>& attributes) const;

private:
  const std::string& getAttributeName(unsigned char typeID) const;
  unsigned char getAttributeTypeID(const std::string& name) const;

protected:
  void checksum(FFaCheckSum* cs) const;

  AttribsVec myAttributes;
};


/*!
  Sorts attributes based on attribute type
*/

struct FFlFEAttribTypeLess
{
  bool operator() (const AttribData& lhs, const AttribData& rhs) const
  { return lhs.first < rhs.first; }

  bool operator() (const AttribData& lhs, unsigned char rhs) const
  { return lhs.first < rhs; }

  bool operator() (unsigned char lhs, const AttribData& rhs) const
  { return lhs < rhs.first; }
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
