// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlSTRCoat.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFaLib/FFaAlgebra/FFaMat33.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"


void FFlSTRCT3::init()
{
  FFlSTRCT3TypeInfoSpec::instance()->setTypeName("STRCT3");
  FFlSTRCT3TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::STRC_ELM);

  ElementFactory::instance()->registerCreator(FFlSTRCT3TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlSTRCT3::create,int,FFlElementBase*&));

  FFlSTRCT3AttributeSpec::instance()->addLegalAttribute("PSTRC",true,true);
  FFlSTRCT3AttributeSpec::instance()->addLegalAttribute("PFATIGUE",false);

  FFlSTRCT3ElementTopSpec::instance()->setNodeCount(3);
  FFlSTRCT3ElementTopSpec::instance()->setNodeDOFs(0);
  FFlSTRCT3ElementTopSpec::instance()->setShellFaces(true);

  static int faces[] = { 1, 2, 3,-1 };
  FFlSTRCT3ElementTopSpec::instance()->setTopology(1,faces);
}


FaMat33 FFlSTRCT3::getGlobalizedElmCS() const
{
  FaMat33 T;
  return T.makeGlobalizedCS(this->getNode(1)->getPos(),
			    this->getNode(2)->getPos(),
			    this->getNode(3)->getPos());
}


void FFlSTRCQ4::init()
{
  FFlSTRCQ4TypeInfoSpec::instance()->setTypeName("STRCQ4");
  FFlSTRCQ4TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::STRC_ELM);

  ElementFactory::instance()->registerCreator(FFlSTRCQ4TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlSTRCQ4::create,int,FFlElementBase*&));

  FFlSTRCQ4AttributeSpec::instance()->addLegalAttribute("PSTRC",true,true);
  FFlSTRCQ4AttributeSpec::instance()->addLegalAttribute("PFATIGUE",false);

  FFlSTRCQ4ElementTopSpec::instance()->setNodeCount(4);
  FFlSTRCQ4ElementTopSpec::instance()->setNodeDOFs(0);
  FFlSTRCQ4ElementTopSpec::instance()->setShellFaces(true);

  static int faces[] = { 1, 2, 3, 4,-1 };
  FFlSTRCQ4ElementTopSpec::instance()->setTopology(1,faces);
}


FaMat33 FFlSTRCQ4::getGlobalizedElmCS() const
{
  FaMat33 T;
  return T.makeGlobalizedCS(this->getNode(1)->getPos(),
			    this->getNode(2)->getPos(),
			    this->getNode(3)->getPos(),
			    this->getNode(4)->getPos());
}


void FFlSTRCT6::init()
{
  FFlSTRCT6TypeInfoSpec::instance()->setTypeName("STRCT6");
  FFlSTRCT6TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::STRC_ELM);

  ElementFactory::instance()->registerCreator(FFlSTRCT6TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlSTRCT6::create,int,FFlElementBase*&));

  FFlSTRCT6AttributeSpec::instance()->addLegalAttribute("PSTRC",true,true);
  FFlSTRCT6AttributeSpec::instance()->addLegalAttribute("PFATIGUE",false);

  FFlSTRCT6ElementTopSpec::instance()->setNodeCount(6);
  FFlSTRCT6ElementTopSpec::instance()->setNodeDOFs(0);
  FFlSTRCT6ElementTopSpec::instance()->setShellFaces(true);

  static int faces[] = { 1, 2, 3, 4, 5, 6,-1 };
  FFlSTRCT6ElementTopSpec::instance()->setTopology(1,faces);
}


FaMat33 FFlSTRCT6::getGlobalizedElmCS() const
{
  FaMat33 T;
  return T.makeGlobalizedCS(this->getNode(1)->getPos(),
			    this->getNode(3)->getPos(),
			    this->getNode(5)->getPos());
}


void FFlSTRCQ8::init()
{
  FFlSTRCQ8TypeInfoSpec::instance()->setTypeName("STRCQ8");
  FFlSTRCQ8TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::STRC_ELM);

  ElementFactory::instance()->registerCreator(FFlSTRCQ8TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlSTRCQ8::create,int,FFlElementBase*&));

  FFlSTRCQ8AttributeSpec::instance()->addLegalAttribute("PSTRC",true,true);
  FFlSTRCQ8AttributeSpec::instance()->addLegalAttribute("PFATIGUE",false);

  FFlSTRCQ8ElementTopSpec::instance()->setNodeCount(8);
  FFlSTRCQ8ElementTopSpec::instance()->setNodeDOFs(0);
  FFlSTRCQ8ElementTopSpec::instance()->setShellFaces(true);

  static int faces[] = { 1, 2, 3, 4, 5, 6, 7, 8,-1 };
  FFlSTRCQ8ElementTopSpec::instance()->setTopology(1,faces);
}


FaMat33 FFlSTRCQ8::getGlobalizedElmCS() const
{
  FaMat33 T;
  return T.makeGlobalizedCS(this->getNode(1)->getPos(),
			    this->getNode(3)->getPos(),
			    this->getNode(5)->getPos(),
			    this->getNode(7)->getPos());
}


bool FFlStrainCoatBase::resolveElmRef(const std::vector<FFlElementBase*>& possibleElms,
				      bool suppressErrmsg)
{
  if (possibleElms.empty())
    return false;

  if (myFElm.resolve(possibleElms))
    return true;

  if (!suppressErrmsg)
    ListUI <<"\n *** Error: Failed to resolve reference to finite element "
	   << myFElm.getID() <<"\n";

  return false;
}


FFlElementBase* FFlStrainCoatBase::getFElement() const
{
  if (myFElm.isResolved())
    return myFElm.getReference();
  else
    return 0;
}
