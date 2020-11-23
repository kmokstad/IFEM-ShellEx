// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPNSM.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"


FFlPNSM::FFlPNSM(int id) : FFlAttributeBase(id)
{
  this->addField(NSM);
  this->addField(isShell);
}


FFlPNSM::FFlPNSM(const FFlPNSM& obj) : FFlAttributeBase(obj)
{
  this->addField(NSM);
  this->addField(isShell);

  NSM.setValue(obj.NSM.getValue());
  isShell.setValue(obj.isShell.getValue());
}


void FFlPNSM::convertUnits(const FFaUnitCalculator* convCal)
{
  if (isShell.getValue() == 1)
    convCal->convert(NSM.data(),"MASS/AREA",10);
  else
    convCal->convert(NSM.data(),"MASS/LENGTH",10);
}


void FFlPNSM::init()
{
  FFlPNSMTypeInfoSpec::instance()->setTypeName("PNSM");
  FFlPNSMTypeInfoSpec::instance()->setDescription("Non-structural masses");
  FFlPNSMTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::MATERIAL_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPNSMTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPNSM::create,int,FFlAttributeBase*&));
}
