// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_FE_ATTRIBUTE_SPEC_H
#define FFL_FE_ATTRIBUTE_SPEC_H

#include <string>
#include <set>

#define FFL_FE_ATTRIBUTE_SPEC(class) \
  virtual FFlFEAttributeSpec* getFEAttributeSpec() const \
  { return FFaSingelton<FFlFEAttributeSpec,class>::instance(); }


class FFlFEAttributeSpec
{
  class AttribType
  {
  public:
    AttribType() = delete; // No default constructor
    AttribType(const std::string& aName, bool req, bool multi)
      : myName(aName), isRequired(req), allowMultiple(multi)
    { typeID = ++typeIDCounter; }

    void static init() { typeIDCounter = 0; }

    std::string myName;
    bool isRequired;
    bool allowMultiple;
    unsigned char typeID;

  private:
    static unsigned char typeIDCounter;

    friend bool operator<(const AttribType& a, const AttribType& b)
    { return a.myName < b.myName; }
  };

public:
  void dump() const;

  bool addLegalAttribute(const std::string& name,
                         bool required = true, bool allowMulti = false);

  static bool isObsolete(const std::string& name);

  bool multipleRefsAllowed(const std::string& name) const;

  const std::string& getAttributeName(unsigned char typeID) const;
  unsigned char getAttributeTypeID(const std::string& name) const;

  static void initTypeID() { AttribType::init(); }

private:
  std::set<AttribType> myLegalAttributes;
};

#endif
