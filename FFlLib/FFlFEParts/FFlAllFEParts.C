// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlAllFEParts.H"
#include "FFlLib/FFlMemPool.H"

//! Set to \e true when initialized, to avoid initializing more than once.
static bool initialized = false;


void FFl::initAllElements()
{
  if (initialized) return;

  FFlFEAttributeSpec::initTypeID();

  FFlNode::init();
  FFlRGD::init();
  FFlRBAR::init();
  FFlWAVGM::init();
  FFlCMASS::init();
  FFlSPRING::init();
  FFlRSPRING::init();
  FFlBUSH::init();
  FFlBEAM2::init();
  FFlBEAM3::init();
  FFlTRI3::init();
  FFlTRI6::init();
  FFlQUAD4::init();
  FFlQUAD8::init();
  FFlTET4::init();
  FFlTET10::init();
  FFlWEDG6::init();
  FFlWEDG15::init();
  FFlHEX8::init();
  FFlHEX20::init();
#ifdef FT_USE_STRAINCOAT
  FFlSTRCT3::init();
  FFlSTRCT6::init();
  FFlSTRCQ4::init();
  FFlSTRCQ8::init();
#endif

  FFlCFORCE::init();
  FFlCMOMENT::init();
  FFlFACELOAD::init();
  FFlSURFLOAD::init();

  FFlPMASS::init();
  FFlPSPRING::init();
  FFlPBUSHCOEFF::init();
  FFlPBUSHECCENT::init();
  FFlPORIENT::init();
  FFlPORIENT3::init();
  FFlPBEAMECCENT::init();
  FFlPBEAMSECTION::init();
  FFlPBEAMPIN::init();
  FFlPEFFLENGTH::init();
  FFlPSPRING::init();
  FFlPTHICK::init();
  FFlPCOMP::init();
  FFlPNSM::init();
  FFlPRGD::init();
  FFlPWAVGM::init();
  FFlPRBAR::init();
  FFlPMAT::init();
  FFlPMAT2D::init();
  FFlPMAT3D::init();
  FFlPMATSHELL::init();
  FFlPCOORDSYS::init();
#ifdef FT_USE_STRAINCOAT
  FFlPSTRC::init();
  FFlPHEIGHT::init();
  FFlPTHICKREF::init();
  FFlPFATIGUE::init();
#endif
#ifdef FT_USE_VISUALS
  FFlVAppearance::init();
  FFlVDetail::init();
#endif

  FFlGroup::init();

  initialized = true;
}


template<class T> void releaseElement()
{
  FFaSingelton<FFlFEElementTopSpec,T>::removeInstance();
  FFaSingelton<FFlFEAttributeSpec,T>::removeInstance();
  FFaSingelton<FFlTypeInfoSpec,T>::removeInstance();
}


template<class T> void releaseAttribute()
{
  FFaSingelton<FFlFEAttributeSpec,T>::removeInstance();
  FFaSingelton<FFlTypeInfoSpec,T>::removeInstance();
}


void FFl::releaseAllElements()
{
  FFlNodeTypeInfoSpec::removeInstance();

  releaseElement<FFlRGD>();
  releaseElement<FFlRBAR>();
  releaseElement<FFlWAVGM>();
  releaseElement<FFlCMASS>();
  releaseElement<FFlSPRING>();
  releaseElement<FFlRSPRING>();
  releaseElement<FFlBUSH>();
  releaseElement<FFlBEAM2>();
  releaseElement<FFlBEAM3>();
  releaseElement<FFlTRI3>();
  releaseElement<FFlTRI6>();
  releaseElement<FFlQUAD4>();
  releaseElement<FFlQUAD8>();
  releaseElement<FFlTET4>();
  releaseElement<FFlTET10>();
  releaseElement<FFlWEDG6>();
  releaseElement<FFlWEDG15>();
  releaseElement<FFlHEX8>();
  releaseElement<FFlHEX20>();
#ifdef FT_USE_STRAINCOAT
  releaseElement<FFlSTRCT3>();
  releaseElement<FFlSTRCT6>();
  releaseElement<FFlSTRCQ4>();
  releaseElement<FFlSTRCQ8>();
#endif
  ElementFactory::removeInstance();

  releaseAttribute<FFlCFORCE>();
  releaseAttribute<FFlCMOMENT>();
  releaseAttribute<FFlFACELOAD>();
  releaseAttribute<FFlSURFLOAD>();
  LoadFactory::removeInstance();

  releaseAttribute<FFlPMASS>();
  releaseAttribute<FFlPSPRING>();
  releaseAttribute<FFlPBUSHCOEFF>();
  releaseAttribute<FFlPBUSHECCENT>();
  releaseAttribute<FFlPORIENT>();
  releaseAttribute<FFlPORIENT3>();
  releaseAttribute<FFlPBEAMECCENT>();
  releaseAttribute<FFlPBEAMSECTION>();
  releaseAttribute<FFlPBEAMPIN>();
  releaseAttribute<FFlPEFFLENGTH>();
  releaseAttribute<FFlPSPRING>();
  releaseAttribute<FFlPTHICK>();
  releaseAttribute<FFlPCOMP>();
  releaseAttribute<FFlPNSM>();
  releaseAttribute<FFlPRGD>();
  releaseAttribute<FFlPWAVGM>();
  releaseAttribute<FFlPRBAR>();
  releaseAttribute<FFlPMAT>();
  releaseAttribute<FFlPMAT2D>();
  releaseAttribute<FFlPMAT3D>();
  releaseAttribute<FFlPMATSHELL>();
  releaseAttribute<FFlPCOORDSYS>();
#ifdef FT_USE_STRAINCOAT
  releaseAttribute<FFlPSTRC>();
  releaseAttribute<FFlPHEIGHT>();
  releaseAttribute<FFlPTHICKREF>();
  releaseAttribute<FFlPFATIGUE>();
#endif
  AttributeFactory::removeInstance();

#ifdef FT_USE_VISUALS
  FFlVAppearanceTypeInfoSpec::removeInstance();
  FFlVDetailTypeInfoSpec::removeInstance();
  VisualFactory::removeInstance();
#endif

  FFlGroupTypeInfoSpec::removeInstance();

  FFlMemPool::deleteAllLinkMemPools();

  initialized = false;
}
