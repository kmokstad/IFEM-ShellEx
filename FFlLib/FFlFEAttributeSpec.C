// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "FFlLib/FFlFEAttributeSpec.H"


unsigned char FFlFEAttributeSpec::AttribType::typeIDCounter = 0;


/*!
  \brief Adds a legal attribute to the object.
  \param[in] name The name of the attribute
  \param[in] required If \e true, the attribute is required in the object
  \param[in] allowMulti If \e true, multiple instances are allowed in the object
*/

bool FFlFEAttributeSpec::addLegalAttribute(const std::string& name,
                                           bool required, bool allowMulti)
{
  return myLegalAttributes.insert(AttribType(name,required,allowMulti)).second;
}


/*!
  Checks if the specified attribute is obsolete and should be silently ignored.
*/

bool FFlFEAttributeSpec::isObsolete(const std::string& attr)
{
  static std::set<std::string> obsAtts;
  if (obsAtts.empty())
  {
    // Enter here all attributes that are no longer in use
    // but may exist in older ftl-files
    obsAtts.insert("PBEAMVISUAL");
  }

  return obsAtts.find(attr) != obsAtts.end();
}


/*!
  Checks if multiple references are allowed for the specified attribute.
*/

bool FFlFEAttributeSpec::multipleRefsAllowed(const std::string& name) const
{
  for (const AttribType& attr : myLegalAttributes)
    if (attr.myName == name)
      return attr.allowMultiple;;

  return false;
}


/*!
  Returns the corresponding name to the input type.
*/

const std::string& FFlFEAttributeSpec::getAttributeName(unsigned char typeID) const
{
  for (const AttribType& attr : myLegalAttributes)
    if (attr.typeID == typeID)
      return attr.myName;

  static std::string unknown("unknown");
  return unknown;
}


/*!
  Returns the corresponding type ID to the input name.
*/

unsigned char FFlFEAttributeSpec::getAttributeTypeID(const std::string& name) const
{
  for (const AttribType& attr : myLegalAttributes)
    if (attr.myName == name)
      return attr.typeID;

  return 0;
}


void FFlFEAttributeSpec::dump() const
{
  std::cout <<"Legal properties:";
  for (const AttribType& attr : myLegalAttributes)
    std::cout <<"\n\t"<< attr.myName <<"\t"<< attr.typeID;
  std::cout << std::endl;
}
