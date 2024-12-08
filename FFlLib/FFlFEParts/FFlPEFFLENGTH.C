// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPEFFLENGTH.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


FFlPEFFLENGTH::FFlPEFFLENGTH(int id) : FFlAttributeBase(id)
{
  this->addField(length);
}


FFlPEFFLENGTH::FFlPEFFLENGTH(const FFlPEFFLENGTH& obj) : FFlAttributeBase(obj)
{
  this->addField(length);
  length.setValue(obj.length.getValue());
}


void FFlPEFFLENGTH::convertUnits(const FFaUnitCalculator* convCal)
{
  convCal->convert(length.data(), "LENGTH");
}


bool FFlPEFFLENGTH::isIdentic(const FFlAttributeBase* otherAttrib) const
{
  const FFlPEFFLENGTH* other = dynamic_cast<const FFlPEFFLENGTH*>(otherAttrib);
  return other ? (this->length == other->length) : false;
}


void FFlPEFFLENGTH::init()
{
  using TypeInfoSpec = FFaSingelton<FFlTypeInfoSpec,FFlPEFFLENGTH>;

  TypeInfoSpec::instance()->setTypeName("PEFFLENGTH");
  TypeInfoSpec::instance()->setDescription("Effective beam lengths");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPEFFLENGTH::create,int,FFlAttributeBase*&));
}

#ifdef FF_NAMESPACE
} // namespace
#endif
