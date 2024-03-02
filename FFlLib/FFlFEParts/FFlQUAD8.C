// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlQUAD8.H"
#include "FFlLib/FFlFEParts/FFlCurvedFace.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlLinkHandler.H"
#include "FFaLib/FFaAlgebra/FFaMat33.H"
#include "FFaLib/FFaAlgebra/FFaVolume.H"


void FFlQUAD8::init()
{
  FFlQUAD8TypeInfoSpec::instance()->setTypeName("QUAD8");
  FFlQUAD8TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SHELL_ELM);

  ElementFactory::instance()->registerCreator(FFlQUAD8TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlQUAD8::create,int,FFlElementBase*&));

  FFlQUAD8AttributeSpec::instance()->addLegalAttribute("PTHICK",false);
  FFlQUAD8AttributeSpec::instance()->addLegalAttribute("PMAT",false);
  FFlQUAD8AttributeSpec::instance()->addLegalAttribute("PMATSHELL",false);
  FFlQUAD8AttributeSpec::instance()->addLegalAttribute("PCOMP",false);
  FFlQUAD8AttributeSpec::instance()->addLegalAttribute("PNSM",false);
  FFlQUAD8AttributeSpec::instance()->addLegalAttribute("PCOORDSYS",false);

  FFlQUAD8ElementTopSpec::instance()->setNodeCount(8);
  FFlQUAD8ElementTopSpec::instance()->setNodeDOFs(6);
  FFlQUAD8ElementTopSpec::instance()->setShellFaces(true);

  static int faces[] = { 1,2,3,4,5,6,7,8,-1 };
  FFlQUAD8ElementTopSpec::instance()->setTopology(1,faces);
}


bool FFlQUAD8::split(Elements& newElem, FFlLinkHandler* owner, int centerNode)
{
  newElem.clear();
  newElem.reserve(4);

  int elm_id = owner->getNewElmID();
  for (int i = 0; i < 4; i++)
    newElem.push_back(ElementFactory::instance()->create("QUAD4",elm_id++));

  newElem[0]->setNode(1,this->getNodeID(1));
  newElem[0]->setNode(2,this->getNodeID(2));
  newElem[0]->setNode(3,centerNode);
  newElem[0]->setNode(4,this->getNodeID(8));
  newElem[1]->setNode(1,this->getNodeID(2));
  newElem[1]->setNode(2,this->getNodeID(3));
  newElem[1]->setNode(3,this->getNodeID(4));
  newElem[1]->setNode(4,centerNode);
  newElem[2]->setNode(1,centerNode);
  newElem[2]->setNode(2,this->getNodeID(4));
  newElem[2]->setNode(3,this->getNodeID(5));
  newElem[2]->setNode(4,this->getNodeID(6));
  newElem[3]->setNode(1,this->getNodeID(8));
  newElem[3]->setNode(2,centerNode);
  newElem[3]->setNode(3,this->getNodeID(6));
  newElem[3]->setNode(4,this->getNodeID(7));

  return true;
}


FaMat33 FFlQUAD8::getGlobalizedElmCS() const
{
  //TODO,bh: use the possible PCOORDSYS here
  FaMat33 T;
  return T.makeGlobalizedCS(this->getNode(1)->getPos(),
			    this->getNode(3)->getPos(),
			    this->getNode(5)->getPos(),
			    this->getNode(7)->getPos());
}


bool FFlQUAD8::getFaceNormals(std::vector<FaVec3>& normals, short int,
			      bool switchNormal) const
{
  std::vector<FFlNode*> nodes;
  if (!this->getFaceNodes(nodes,1,switchNormal)) return false;

  return FFlCurvedFace::faceNormals(nodes,normals);
}


bool FFlQUAD8::getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const
{
  double th = this->getThickness();
  if (th <= 0.0) return false; // Should not happen

  //TODO,kmo: Account for curved edges and face (by numerical integration)
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(3)->getPos());
  FaVec3 v3(this->getNode(5)->getPos());
  FaVec3 v4(this->getNode(7)->getPos());

  double a1 = ((v2-v1) ^ (v3-v1)).length();
  double a2 = ((v3-v1) ^ (v4-v1)).length();

  volume = 0.5*(a1+a2)*th;
  cog = ((v1+v3)*(a1+a2) + v2*a1 + v4*a2) / ((a1+a2)*3.0);

  // Compute inertias by expanding the shell into a solid

  FaVec3 normal = (v3-v1) ^ (v4-v2);
  double length = normal.length();
  if (length < 1.0e-16) return false;

  normal *= th/length;
  v1 -= cog + 0.5*normal;
  v2 -= cog + 0.5*normal;
  v3 -= cog + 0.5*normal;
  v4 -= cog + 0.5*normal;

  FaVec3 v5(v1+normal);
  FaVec3 v6(v2+normal);
  FaVec3 v7(v3+normal);
  FaVec3 v8(v4+normal);

  FFaVolume::hexMoment(v1,v2,v3,v4,v5,v6,v7,v8,inertia);
  return true;
}
