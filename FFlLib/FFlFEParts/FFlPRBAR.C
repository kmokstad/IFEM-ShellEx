// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlPRBAR.H"


FFlPRBAR::FFlPRBAR(int id) : FFlAttributeBase(id)
{
  this->addField(CNA);
  this->addField(CNB);
  this->addField(CMA);
  this->addField(CMB);

  CNA = 0;
  CNB = 0;
  CMA = 0;
  CMB = 0;
}


FFlPRBAR::FFlPRBAR(const FFlPRBAR& obj) : FFlAttributeBase(obj)
{
  this->addField(CNA);
  this->addField(CNB);
  this->addField(CMA);
  this->addField(CMB);

  CNA = obj.CNA;
  CNB = obj.CNB;
  CMA = obj.CMA;
  CMB = obj.CMB;
}


bool FFlPRBAR::isIdentic(const FFlAttributeBase* otherAttrib) const
{
  const FFlPRBAR* other = dynamic_cast<const FFlPRBAR*>(otherAttrib);
  if (!other) return false;

  if (this->CNA != other->CNA) return false;
  if (this->CNB != other->CNB) return false;
  if (this->CMA != other->CMA) return false;
  if (this->CMB != other->CMB) return false;

  return true;
}


void FFlPRBAR::init()
{
  FFlPRBARTypeInfoSpec::instance()->setTypeName("PRBAR");
  FFlPRBARTypeInfoSpec::instance()->setDescription("Rigid bar properties");
  FFlPRBARTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::OTHER_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPRBARTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPRBAR::create,int,FFlAttributeBase*&));
}
