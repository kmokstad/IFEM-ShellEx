// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlTET4.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include "FFaLib/FFaAlgebra/FFaVolume.H"


void FFlTET4::init()
{
  FFlTET4TypeInfoSpec::instance()->setTypeName("TET4");
  FFlTET4TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SOLID_ELM);

  ElementFactory::instance()->registerCreator(FFlTET4TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlTET4::create,int,FFlElementBase*&));

  FFlTET4AttributeSpec::instance()->addLegalAttribute("PMAT");

  FFlTET4ElementTopSpec::instance()->setNodeCount(4);
  FFlTET4ElementTopSpec::instance()->setNodeDOFs(3);

  static int faces[] = { 1, 3, 2,-1,
			 1, 2, 4,-1,
			 2, 3, 4,-1,
			 1, 4, 3,-1 };

  FFlTET4ElementTopSpec::instance()->setTopology(4,faces);
}


bool FFlTET4::getFaceNormals(std::vector<FaVec3>& normals, short int face,
			     bool switchNormal) const
{
  std::vector<FFlNode*> nodes;
  if (!this->getFaceNodes(nodes,face,switchNormal)) return false;

  FaVec3 v0 =  nodes[0]->getPos();
  FaVec3 vn = (nodes[1]->getPos() - v0) ^ (nodes[2]->getPos() - v0);

  normals.clear();
  normals.push_back(vn.normalize());
  normals.push_back(vn);
  normals.push_back(vn);

  return true;
}


bool FFlTET4::getVolumeAndInertia(double& volume, FaVec3& cog,
				  FFaTensor3& inertia) const
{
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(2)->getPos());
  FaVec3 v3(this->getNode(3)->getPos());
  FaVec3 v4(this->getNode(4)->getPos());

  FFaVolume::tetVolume(v1,v2,v3,v4,volume);
  FFaVolume::tetCenter(v1,v2,v3,v4,cog);
  FFaVolume::tetMoment(v1,v2,v3,v4,inertia);

  return true;
}


int FFlTET4::checkOrientation(bool fixIt)
{
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(2)->getPos());
  FaVec3 v3(this->getNode(3)->getPos());
  FaVec3 v4(this->getNode(4)->getPos());

  double volume;
  FFaVolume::tetVolume(v1,v2,v3,v4,volume);
  if (volume >= 1.0e-16) return 1; // ok orientation

  // Try swap nodes 2 and 3
  FFaVolume::tetVolume(v1,v3,v2,v4,volume);
  if (volume < 1.0e-16) return 0; // flat element

  if (fixIt)
  {
    FFlNode* node2 = this->getNode(2);
    FFlNode* node3 = this->getNode(3);
    this->setNode(2,node3);
    this->setNode(3,node2);
  }
  return -1;
}
