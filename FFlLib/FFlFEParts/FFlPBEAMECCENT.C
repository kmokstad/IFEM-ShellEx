// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPBEAMECCENT.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


FFlPBEAMECCENT::FFlPBEAMECCENT(int id) : FFlAttributeBase(id)
{
  this->addField(node1Offset);
  this->addField(node2Offset);
}


FFlPBEAMECCENT::FFlPBEAMECCENT(const FFlPBEAMECCENT& ob) : FFlAttributeBase(ob)
{
  if (ob.size() > 0) this->addField(node1Offset);
  if (ob.size() > 1) this->addField(node2Offset);
  if (ob.size() > 2) this->addField(node3Offset);

  node1Offset.setValue(ob.node1Offset.getValue());
  node2Offset.setValue(ob.node2Offset.getValue());
  node3Offset.setValue(ob.node3Offset.getValue());
}


void FFlPBEAMECCENT::convertUnits(const FFaUnitCalculator* convCal)
{
  convCal->convert(node1Offset.data(),"LENGTH");
  convCal->convert(node2Offset.data(),"LENGTH");
  convCal->convert(node3Offset.data(),"LENGTH");
}


void FFlPBEAMECCENT::init()
{
  using TypeInfoSpec = FFaSingelton<FFlTypeInfoSpec,FFlPBEAMECCENT>;

  TypeInfoSpec::instance()->setTypeName("PBEAMECCENT");
  TypeInfoSpec::instance()->setDescription("Beam eccentricities");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPBEAMECCENT::create,int,FFlAttributeBase*&));
}


void FFlPBEAMECCENT::resize(size_t n)
{
  if (n == 9 && this->size() == 2)
    this->addField(node3Offset);
  else
    this->FFlAttributeBase::resize(n/3);
}

#ifdef FF_NAMESPACE
} // namespace
#endif
