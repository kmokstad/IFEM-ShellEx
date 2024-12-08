// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPORIENT.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


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
  using TypeInfoSpec = FFaSingelton<FFlTypeInfoSpec,FFlPORIENT>;

  TypeInfoSpec::instance()->setTypeName("PORIENT");
  TypeInfoSpec::instance()->setDescription("Orientation vectors");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
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
  using TypeInfoSpec = FFaSingelton<FFlTypeInfoSpec,FFlPORIENT3>;

  TypeInfoSpec::instance()->setTypeName("PORIENT3");
  TypeInfoSpec::instance()->setDescription("Orientation vectors");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPORIENT3::create,int,FFlAttributeBase*&));
}

#ifdef FF_NAMESPACE
} // namespace
#endif
