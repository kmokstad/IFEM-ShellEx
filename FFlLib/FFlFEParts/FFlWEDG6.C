// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlWEDG6.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include "FFaLib/FFaAlgebra/FFaVolume.H"


void FFlWEDG6::init()
{
  FFlWEDG6TypeInfoSpec::instance()->setTypeName("WEDG6");
  FFlWEDG6TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SOLID_ELM);

  ElementFactory::instance()->registerCreator(FFlWEDG6TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlWEDG6::create,int,FFlElementBase*&));

  FFlWEDG6AttributeSpec::instance()->addLegalAttribute("PMAT");

  FFlWEDG6ElementTopSpec::instance()->setNodeCount(6);
  FFlWEDG6ElementTopSpec::instance()->setNodeDOFs(3);

  static int faces[] = { 1, 2, 5, 4,-1,
			 2, 3, 6, 5,-1,
			 1, 4, 6, 3,-1,
			 1, 3, 2,   -1,
			 4, 5, 6,   -1 };

  FFlWEDG6ElementTopSpec::instance()->setTopology(5,faces);
}


bool FFlWEDG6::getFaceNormals(std::vector<FaVec3>& normals, short int face,
			      bool switchNormal) const
{
  std::vector<FFlNode*> nodes;
  if (!this->getFaceNodes(nodes,face,switchNormal)) return false;

  FaVec3 v0 = nodes[0]->getPos();
  FaVec3 vn = nodes.size() == 4 ?
    (nodes[2]->getPos() - v0) ^ (nodes[3]->getPos() - nodes[1]->getPos()) :
    (nodes[1]->getPos() - v0) ^ (nodes[2]->getPos() - v0);

  normals.clear();
  normals.push_back(vn.normalize());
  for (size_t i = 1; i < nodes.size(); i++)
    normals.push_back(vn);

  return true;
}


bool FFlWEDG6::getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const
{
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(2)->getPos());
  FaVec3 v3(this->getNode(3)->getPos());
  FaVec3 v4(this->getNode(4)->getPos());
  FaVec3 v5(this->getNode(5)->getPos());
  FaVec3 v6(this->getNode(6)->getPos());

  FFaVolume::wedVolume(v1,v2,v3,v4,v5,v6,volume);
  FFaVolume::wedCenter(v1,v2,v3,v4,v5,v6,cog);
  FFaVolume::wedMoment(v1,v2,v3,v4,v5,v6,inertia);

  return true;
}
