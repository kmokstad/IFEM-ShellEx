// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlTRI6.H"
#include "FFlLib/FFlFEParts/FFlCurvedFace.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlLinkHandler.H"
#include "FFaLib/FFaAlgebra/FFaMat33.H"
#include "FFaLib/FFaAlgebra/FFaVolume.H"


void FFlTRI6::init()
{
  FFlTRI6TypeInfoSpec::instance()->setTypeName("TRI6");
  FFlTRI6TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SHELL_ELM);

  ElementFactory::instance()->registerCreator(FFlTRI6TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlTRI6::create,int,FFlElementBase*&));

  FFlTRI6AttributeSpec::instance()->addLegalAttribute("PTHICK",false);
  FFlTRI6AttributeSpec::instance()->addLegalAttribute("PMAT",false);
  FFlTRI6AttributeSpec::instance()->addLegalAttribute("PMATSHELL",false);
  FFlTRI6AttributeSpec::instance()->addLegalAttribute("PCOMP",false);
  FFlTRI6AttributeSpec::instance()->addLegalAttribute("PNSM",false);
  FFlTRI6AttributeSpec::instance()->addLegalAttribute("PCOORDSYS",false);

  FFlTRI6ElementTopSpec::instance()->setNodeCount(6);
  FFlTRI6ElementTopSpec::instance()->setNodeDOFs(6);
  FFlTRI6ElementTopSpec::instance()->setShellFaces(true);

  static int faces[] = { 1,2,3,4,5,6,-1 };
  FFlTRI6ElementTopSpec::instance()->setTopology(1,faces);
}


bool FFlTRI6::split(Elements& newElem, FFlLinkHandler* owner, int)
{
  newElem.clear();
  newElem.reserve(4);

  int elm_id = owner->getNewElmID();
  for (int i = 0; i < 4; i++)
    newElem.push_back(ElementFactory::instance()->create("TRI3",elm_id++));

  newElem[0]->setNode(1,this->getNodeID(1));
  newElem[0]->setNode(2,this->getNodeID(2));
  newElem[0]->setNode(3,this->getNodeID(6));
  newElem[1]->setNode(1,this->getNodeID(2));
  newElem[1]->setNode(2,this->getNodeID(3));
  newElem[1]->setNode(3,this->getNodeID(4));
  newElem[2]->setNode(1,this->getNodeID(6));
  newElem[2]->setNode(2,this->getNodeID(4));
  newElem[2]->setNode(3,this->getNodeID(5));
  newElem[3]->setNode(1,this->getNodeID(4));
  newElem[3]->setNode(2,this->getNodeID(6));
  newElem[3]->setNode(3,this->getNodeID(2));

  return true;
}


FaMat33 FFlTRI6::getGlobalizedElmCS() const
{
  //TODO,bh: use the possible PCOORDSYS here
  FaMat33 T;
  return T.makeGlobalizedCS(this->getNode(1)->getPos(),
			    this->getNode(3)->getPos(),
			    this->getNode(5)->getPos());
}


bool FFlTRI6::getFaceNormals(std::vector<FaVec3>& normals, short int,
			     bool switchNormal) const
{
  std::vector<FFlNode*> nodes;
  if (!this->getFaceNodes(nodes,1,switchNormal)) return false;

  return FFlCurvedFace::faceNormals(nodes,normals);
}


bool FFlTRI6::getVolumeAndInertia(double& volume, FaVec3& cog,
				  FFaTensor3& inertia) const
{
  double th = this->getThickness();
  if (th <= 0.0) return false; // Should not happen

  //TODO,kmo: Account for curved edges and face (by numerical integration)
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(3)->getPos());
  FaVec3 v3(this->getNode(5)->getPos());

  FaVec3 normal = (v2-v1) ^ (v3-v1);
  double length = normal.length();
  volume = 0.5 * length * th;
  cog = (v1 + v2 + v3) / 3.0;
  if (length < 1.0e-16) return false;

  // Compute inertias by expanding the shell into a solid

  normal *= th/length;
  v1 -= cog + 0.5*normal;
  v2 -= cog + 0.5*normal;
  v3 -= cog + 0.5*normal;

  FaVec3 v4(v1+normal);
  FaVec3 v5(v2+normal);
  FaVec3 v6(v3+normal);

  FFaVolume::wedMoment(v1,v2,v3,v4,v5,v6,inertia);
  return true;
}
