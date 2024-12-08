// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPMASS.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


FFlPMASS::FFlPMASS(int id) : FFlAttributeBase(id)
{
  this->addField(M);
}


FFlPMASS::FFlPMASS(const FFlPMASS& obj) : FFlAttributeBase(obj)
{
  this->addField(M);
  M = obj.M;
}


void FFlPMASS::convertUnits(const FFaUnitCalculator* convCal)
{
  int i, j;
  DoubleVec::iterator mit = M.data().begin();
  for (i = 0; i < 6; i++)
    for (j = 0; j <= i; j++)
      if (mit == M.data().end())
        return;
      else if (i < 3 && j < 3)
        convCal->convert(*(mit++),"MASS",10);
      else if (j < 3)
        convCal->convert(*(mit++),"MASS*LENGTH",10);
      else
        convCal->convert(*(mit++),"MASS*AREA",10);
}


void FFlPMASS::calculateChecksum(FFaCheckSum* cs, int csMask) const
{
  this->FFlAttributeBase::calculateChecksum(cs, csMask);

  // Add zeroes to fill up a full 6x6 matrix to match older files
  for (int i = M.getValue().size(); i < 21; i++)
    cs->add((double)0.0);
}


void FFlPMASS::init()
{
  using TypeInfoSpec = FFaSingelton<FFlTypeInfoSpec,FFlPMASS>;

  TypeInfoSpec::instance()->setTypeName("PMASS");
  TypeInfoSpec::instance()->setDescription("Concentrated mass properties");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::MASS_PROP);

  AttributeFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPMASS::create,int,FFlAttributeBase*&));
}

#ifdef FF_NAMESPACE
} // namespace
#endif
