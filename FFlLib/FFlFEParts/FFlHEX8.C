// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlHEX8.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include "FFaLib/FFaAlgebra/FFaVolume.H"


void FFlHEX8::init()
{
  FFlHEX8TypeInfoSpec::instance()->setTypeName("HEX8");
  FFlHEX8TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SOLID_ELM);

  ElementFactory::instance()->registerCreator(FFlHEX8TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlHEX8::create,int,FFlElementBase*&));

  FFlHEX8AttributeSpec::instance()->addLegalAttribute("PMAT");

  FFlHEX8ElementTopSpec::instance()->setNodeCount(8);
  FFlHEX8ElementTopSpec::instance()->setNodeDOFs(3);

  static int faces[] = { 2, 3, 7, 6,-1,
			 3, 4, 8, 7,-1,
			 1, 5, 8, 4,-1,
			 1, 2, 6, 5,-1,
			 5, 6, 7, 8,-1,
			 1, 4, 3, 2,-1 };

  FFlHEX8ElementTopSpec::instance()->setTopology(6,faces);
}


bool FFlHEX8::getFaceNormals(std::vector<FaVec3>& normals, short int face,
			     bool switchNormal) const
{
  std::vector<FFlNode*> nodes;
  if (!this->getFaceNodes(nodes,face,switchNormal)) return false;

  FaVec3 vn = (nodes[2]->getPos() - nodes[0]->getPos()) ^
	      (nodes[3]->getPos() - nodes[1]->getPos());

  normals.clear();
  normals.push_back(vn.normalize());
  normals.push_back(vn);
  normals.push_back(vn);
  normals.push_back(vn);

  return true;
}


bool FFlHEX8::getVolumeAndInertia(double& volume, FaVec3& cog,
				  FFaTensor3& inertia) const
{
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(2)->getPos());
  FaVec3 v3(this->getNode(3)->getPos());
  FaVec3 v4(this->getNode(4)->getPos());
  FaVec3 v5(this->getNode(5)->getPos());
  FaVec3 v6(this->getNode(6)->getPos());
  FaVec3 v7(this->getNode(7)->getPos());
  FaVec3 v8(this->getNode(8)->getPos());

  FFaVolume::hexVolume(v1,v2,v3,v4,v5,v6,v7,v8,volume);
  FFaVolume::hexCenter(v1,v2,v3,v4,v5,v6,v7,v8,cog);
  FFaVolume::hexMoment(v1,v2,v3,v4,v5,v6,v7,v8,inertia);

  return true;
}
