// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_FE_ATTRIBUTE_REFS_H
#define FFL_FE_ATTRIBUTE_REFS_H

#include "FFlReference.H"
#include <vector>
#include <string>
#include <map>

class FFlAttributeBase;
class FFlFEAttributeSpec;
class FFaCheckSum;

typedef FFlReference<FFlAttributeBase>     AttribRef;
typedef std::pair<unsigned char,AttribRef> AttribData;
typedef std::vector<AttribData>            AttribsVec;
typedef AttribsVec::const_iterator         AttribsCIter;
typedef std::map<int,FFlAttributeBase*>    AttributeMap;


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

  bool resolve(const std::map<std::string,AttributeMap>& possibleReferences,
	       bool suppressErrmsg = false);

  virtual FFlFEAttributeSpec* getFEAttributeSpec() const = 0;

  AttribsCIter attributesBegin() const { return myAttributes.begin(); }
  AttribsCIter attributesEnd()   const { return myAttributes.end(); }

  FFlAttributeBase*              getAttribute (const std::string& type) const;
  std::vector<FFlAttributeBase*> getAttributes(const std::string& type) const;

  int getAttributeID(const std::string& type) const;

  bool hasAttribute(const FFlAttributeBase* attribute) const;
  bool hasAttribute(const std::vector<FFlAttributeBase*>& attributes) const;

protected:
  void checksum(FFaCheckSum* cs) const;

  AttribsVec myAttributes;
};


/*!
  Sorts attributes based on attribute type
*/

class FFlFEAttribTypeLess
{
public:
  bool operator() (const AttribData& lhs, const AttribData& rhs) const
  { return keyLess(lhs.first, rhs.first); }

  bool operator() (const AttribData& lhs, const AttribData::first_type& rhs) const
  { return keyLess(lhs.first, rhs); }

  bool operator() (const AttribData::first_type& lhs, const AttribData& rhs) const
  { return keyLess(lhs, rhs.first); }

private:
  bool keyLess (const AttribData::first_type& k1, const AttribData::first_type& k2) const
  { return k1 < k2; }
};

#endif
