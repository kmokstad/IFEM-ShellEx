// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlPFATIGUE.H"


FFlPFATIGUE::FFlPFATIGUE(int id) : FFlAttributeBase(id)
{
  this->addField(snCurveStd);
  this->addField(snCurveIndex);
  this->addField(stressConcentrationFactor);

  snCurveStd = 0;
  snCurveIndex = 0;
  stressConcentrationFactor = 1.0;
}


FFlPFATIGUE::FFlPFATIGUE(const FFlPFATIGUE& obj) : FFlAttributeBase(obj)
{
  this->addField(snCurveStd);
  this->addField(snCurveIndex);
  this->addField(stressConcentrationFactor);

  snCurveStd = obj.snCurveStd;
  snCurveIndex = obj.snCurveIndex;
  stressConcentrationFactor = obj.stressConcentrationFactor;
}


bool FFlPFATIGUE::isIdentic(const FFlAttributeBase* otherAttrib) const
{
  const FFlPFATIGUE* other = dynamic_cast<const FFlPFATIGUE*>(otherAttrib);
  if (!other) return false;

  if (this->snCurveStd != other->snCurveStd) return false;
  if (this->snCurveIndex != other->snCurveIndex) return false;
  if (this->stressConcentrationFactor != other->stressConcentrationFactor) return false;

  return true;
}


void FFlPFATIGUE::init()
{
  FFlPFATIGUETypeInfoSpec::instance()->setTypeName("PFATIGUE");
  FFlPFATIGUETypeInfoSpec::instance()->setDescription("Fatigue properties");
  FFlPFATIGUETypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::STRC_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPFATIGUETypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPFATIGUE::create,int,FFlAttributeBase*&));
}
