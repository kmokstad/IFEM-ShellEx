// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlPCOMP.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


FFlPCOMP::FFlPCOMP(int id) : FFlAttributeBase(id)
{
  this->addField(Z0);
  this->addField(plySet);
}


FFlPCOMP::FFlPCOMP(const FFlPCOMP& obj) : FFlAttributeBase(obj)
{
  this->addField(Z0);
  Z0.setValue(obj.Z0.getValue());

  this->addField(plySet);
  plySet.setValue(obj.plySet.getValue());
}


void FFlPCOMP::convertUnits(const FFaUnitCalculator* convCal)
{
  // Round to 10 significant digits
  convCal->convert(Z0.data(),"LENGTH",10);

  for (FFlPly& ply : plySet.data())
    convCal->convert(ply.T,"LENGTH",10);
}


void FFlPCOMP::init()
{
  using TypeInfoSpec = FFaSingelton<FFlTypeInfoSpec,FFlPCOMP>;

  TypeInfoSpec::instance()->setTypeName("PCOMP");
  TypeInfoSpec::instance()->setDescription("Composite Shell properties");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (TypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPCOMP::create,int,FFlAttributeBase*&));
}


void FFlPCOMP::calculateChecksum(FFaCheckSum* cs, int csMask) const
{
  FFlAttributeBase::calculateChecksum(cs, csMask);

  Z0.calculateChecksum(cs);
  plySet.calculateChecksum(cs);

  return;
}

#ifdef FF_NAMESPACE
} // namespace
#endif


// Explicit instantiation
template FFlField<FFlPlyVec>::FFlField();

template<> inline bool
FFlField<FFlPlyVec>::parse(std::vector<std::string>::const_iterator& begin,
                           const std::vector<std::string>::const_iterator& end)
{
  while (begin != end)
  {
    FFlPly ply;
    if (!parseNumericField(ply.MID,*begin++))   return false;
    if ( begin == end ) return false;
    if (!parseNumericField(ply.T,*begin++))     return false;
    if ( begin == end ) return false;
    if (!parseNumericField(ply.theta,*begin++)) return false;

    myData.push_back(ply);
  }

  return true;
}

template<> inline void
FFlField<FFlPlyVec>::write(std::ostream& os) const
{
  for (const FFlPly& ply : myData)
    os <<"   "<< ply.MID <<" "<< ply.T <<" "<< ply.theta;
}


std::ostream& operator<<(std::ostream& os, const FFlPlyVec& val)
{
  for (const FFlPly& ply : val)
    os <<"\n"<< ply.MID <<" "<< ply.T <<" "<< ply.theta;

  return os;
}


template<> void FFaCheckSum::add(const FFlPlyVec& val)
{
  for (const FFlPly& ply : val)
  {
    this->add(ply.MID);
    this->add(ply.T);
    this->add(ply.theta);
  }
}
