// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlPCOMP.H"
#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"


template<> void FFaCheckSum::add(const FFlPCOMP::PlyVec& v)
{
  for (const FFlPCOMP::Ply& ply : v)
  {
    this->add(ply.MID);
    this->add(ply.T);
    this->add(ply.thetaInDeg);
  }
}

FFlPCOMP::FFlPCOMP(int id) : FFlAttributeBase(id)
{
  this->addField(Z0);
  this->addField(plySet);
  Z0 = 0.0;
}


FFlPCOMP::FFlPCOMP(const FFlPCOMP& obj) : FFlAttributeBase(obj)
{
  this->addField(Z0);
  this->addField(plySet);
  Z0 = obj.Z0;

  for (const Ply& ply : obj.plySet.getValue())
    plySet.data().push_back(ply);
}


void FFlPCOMP::convertUnits(const FFaUnitCalculator* convCal)
{
  // Round to 10 significant digits
  convCal->convert(Z0.data(),"LENGTH",10);

  for (Ply& ply : plySet.data())
    convCal->convert(ply.T,"LENGTH",10);
}


void FFlPCOMP::init()
{
  FFlPCOMPTypeInfoSpec::instance()->setTypeName("PCOMP");
  FFlPCOMPTypeInfoSpec::instance()->setDescription("Composite Shell properties");
  FFlPCOMPTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::GEOMETRY_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPCOMPTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPCOMP::create,int,FFlAttributeBase*&));
}

void FFlPCOMP::calculateChecksum(FFaCheckSum* cs, int csMask) const
{
  FFlAttributeBase::calculateChecksum(cs, csMask);

  Z0.calculateChecksum(cs);
  plySet.calculateChecksum(cs);

  return;
}


//Explicit instantiation
template FFlField<FFlPCOMP::PlyVec>::FFlField();

template<> inline bool
FFlField<FFlPCOMP::PlyVec>::parse(std::vector<std::string>::const_iterator& begin,
				  const std::vector<std::string>::const_iterator& end)
{
  while (begin != end)
  {
    FFlPCOMP::Ply ply;
    if (!parseNumericField(ply.MID,*begin++))   return false;
    if ( begin == end ) return false;
    if (!parseNumericField(ply.T,*begin++))     return false;
    if ( begin == end ) return false;
    if (!parseNumericField(ply.thetaInDeg,*begin++)) return false;

    myData.push_back(ply);
  }

  return true;
}

template<> inline void
FFlField<FFlPCOMP::PlyVec>::write(std::ostream& os) const
{
  for (const FFlPCOMP::Ply& ply : myData)
    os <<"   "<< ply.MID <<" "<< ply.T <<" "<< ply.thetaInDeg;
}

std::ostream& operator<< (std::ostream& os, const FFlPCOMP::PlyVec& val)
{
  for (const FFlPCOMP::Ply& ply : val)
    os <<"\n"<< ply.MID <<" "<< ply.T <<" "<< ply.thetaInDeg;

  return os;
}
