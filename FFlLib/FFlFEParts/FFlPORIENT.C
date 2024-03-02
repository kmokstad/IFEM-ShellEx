// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPORIENT.H"


FFlPORIENT::FFlPORIENT(int id) : FFlAttributeBase(id)
{
  this->addField(directionVector);

  directionVector.setValue(FaVec3(1.0,0.0,0.0));
}


FFlPORIENT::FFlPORIENT(const FFlPORIENT& obj) : FFlAttributeBase(obj)
{
  this->addField(directionVector);

  directionVector = obj.directionVector;
}


bool FFlPORIENT::isIdentic(const FFlAttributeBase* otherAttrib) const
{
  const FFlPORIENT* other = dynamic_cast<const FFlPORIENT*>(otherAttrib);
  return other ? (this->directionVector == other->directionVector) : false;
}


void FFlPORIENT::init()
{
  FFlPORIENTTypeInfoSpec::instance()->setTypeName("PORIENT");
  FFlPORIENTTypeInfoSpec::instance()->setDescription("Orientation vectors");
  FFlPORIENTTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPORIENTTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPORIENT::create,int,FFlAttributeBase*&));

  // Obsolete field names that should be converted into a PORIENT field
  // when reading old link data files
  AttributeFactory::instance()->registerCreator
    ("PBEAMORIENT", FFaDynCB2S(FFlPORIENT::create,int,FFlAttributeBase*&));
  AttributeFactory::instance()->registerCreator
    ("PBUSHORIENT", FFaDynCB2S(FFlPORIENT::create,int,FFlAttributeBase*&));
}


FFlPORIENT3::FFlPORIENT3(int id) : FFlAttributeBase(id)
{
  for (int i = 0; i < 3; i++)
  {
    this->addField(directionVector[i]);
    directionVector[i].setValue(FaVec3(1.0,0.0,0.0));
  }
}


FFlPORIENT3::FFlPORIENT3(const FFlPORIENT3& obj) : FFlAttributeBase(obj)
{
  for (int i = 0; i < 3; i++)
    this->addField(directionVector[i]);

  directionVector = obj.directionVector;
}


bool FFlPORIENT3::isIdentic(const FFlAttributeBase* otherAttrib) const
{
  const FFlPORIENT3* other = dynamic_cast<const FFlPORIENT3*>(otherAttrib);
  if (!other) return false;

  for (int i = 0; i < 3; i++)
    if (this->directionVector[i] != other->directionVector[i])
      return false;

  return true;
}


void FFlPORIENT3::init()
{
  FFlPORIENT3TypeInfoSpec::instance()->setTypeName("PORIENT3");
  FFlPORIENT3TypeInfoSpec::instance()->setDescription("Orientation vectors");
  FFlPORIENT3TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPORIENT3TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPORIENT3::create,int,FFlAttributeBase*&));
}
