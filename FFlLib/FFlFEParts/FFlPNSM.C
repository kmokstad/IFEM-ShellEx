// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlPNSM.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"


FFlPNSM::FFlPNSM(int id) : FFlAttributeBase(id)
{
  this->addField(NSM);
  this->addField(isShell);

  NSM = 0.0;
  isShell = 0;
}


FFlPNSM::FFlPNSM(const FFlPNSM& obj) : FFlAttributeBase(obj)
{
  this->addField(NSM);
  this->addField(isShell);

  NSM = obj.NSM;
  isShell = obj.isShell;
}


void FFlPNSM::convertUnits(const FFaUnitCalculator* convCal)
{
  if (isShell.data() == 1)
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
