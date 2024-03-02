// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlTET10.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlFEParts/FFlCurvedFace.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include "FFaLib/FFaAlgebra/FFaVolume.H"


void FFlTET10::init()
{
  FFlTET10TypeInfoSpec::instance()->setTypeName("TET10");
  FFlTET10TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SOLID_ELM);

  ElementFactory::instance()->registerCreator(FFlTET10TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlTET10::create,int,FFlElementBase*&));

  FFlTET10AttributeSpec::instance()->addLegalAttribute("PMAT");

  FFlTET10ElementTopSpec::instance()->setNodeCount(10);
  FFlTET10ElementTopSpec::instance()->setNodeDOFs(3);

  static int faces[] = { 1, 6, 5, 4, 3, 2,-1,
			 3, 4, 5, 9,10, 8,-1,
			 1, 7,10, 9, 5, 6,-1,
    			 1, 2, 3, 8,10, 7,-1 };

  FFlTET10ElementTopSpec::instance()->setTopology(4,faces);
}


bool FFlTET10::getFaceNormals(std::vector<FaVec3>& normals, short int face,
			      bool switchNormal) const
{
  std::vector<FFlNode*> nodes;
  if (!this->getFaceNodes(nodes,face,switchNormal)) return false;

  return FFlCurvedFace::faceNormals(nodes,normals);
}


bool FFlTET10::getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const
{
  //TODO: Account for possibly curved edges
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(3)->getPos());
  FaVec3 v3(this->getNode(5)->getPos());
  FaVec3 v4(this->getNode(10)->getPos());

  FFaVolume::tetVolume(v1,v2,v3,v4,volume);
  FFaVolume::tetCenter(v1,v2,v3,v4,cog);
  FFaVolume::tetMoment(v1,v2,v3,v4,inertia);

  return true;
}


int FFlTET10::checkOrientation(bool fixIt)
{
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(3)->getPos());
  FaVec3 v3(this->getNode(5)->getPos());
  FaVec3 v4(this->getNode(10)->getPos());

  double volume;
  FFaVolume::tetVolume(v1,v2,v3,v4,volume);
  if (volume >= 1.0e-16) return 1; // ok orientation

  // Try swap nodes 2 and 3
  FFaVolume::tetVolume(v1,v3,v2,v4,volume);
  if (volume < 1.0e-16) return 0; // flat element

  if (fixIt)
  {
    FFlNode* node3 = this->getNode(3);
    FFlNode* node5 = this->getNode(5);
    this->setNode(3,node5);
    this->setNode(5,node3);
    FFlNode* node2 = this->getNode(2);
    FFlNode* node6 = this->getNode(6);
    this->setNode(2,node6);
    this->setNode(6,node2);
    FFlNode* node8 = this->getNode(8);
    FFlNode* node9 = this->getNode(9);
    this->setNode(8,node9);
    this->setNode(9,node8);
  }
  return -1;
}
