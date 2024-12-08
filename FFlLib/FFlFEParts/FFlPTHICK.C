// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPTHICK.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


FFlPTHICK::FFlPTHICK(int id) : FFlAttributeBase(id)
{
  this->addField(thickness);
  thickness = 0.1;
}


FFlPTHICK::FFlPTHICK(const FFlPTHICK& obj) : FFlAttributeBase(obj)
{
  this->addField(thickness);
  thickness = obj.thickness;
}


void FFlPTHICK::convertUnits(const FFaUnitCalculator* convCal)
{
  // Round to 10 significant digits
  convCal->convert(thickness.data(),"LENGTH",10);
}


void FFlPTHICK::init()
{
  using TypeInfoSpec = FFaSingelton<FFlTypeInfoSpec,FFlPTHICK>;

  TypeInfoSpec::instance()->setTypeName("PTHICK");
  TypeInfoSpec::instance()->setDescription("Shell thicknesses");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPTHICK::create,int,FFlAttributeBase*&));
}

#ifdef FF_NAMESPACE
} // namespace
#endif
