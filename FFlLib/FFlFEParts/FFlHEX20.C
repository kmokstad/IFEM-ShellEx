// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlHEX20.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlFEParts/FFlCurvedFace.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include "FFaLib/FFaAlgebra/FFaVolume.H"


void FFlHEX20::init()
{
  FFlHEX20TypeInfoSpec::instance()->setTypeName("HEX20");
  FFlHEX20TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SOLID_ELM);
  ElementFactory::instance()->registerCreator(FFlHEX20TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlHEX20::create,int,FFlElementBase*&));

  FFlHEX20AttributeSpec::instance()->addLegalAttribute("PMAT");

  FFlHEX20ElementTopSpec::instance()->setNodeCount(20);
  FFlHEX20ElementTopSpec::instance()->setNodeDOFs(3);

  static int faces[] = { 1, 2, 3,10,15,14,13, 9,-1,
			 3, 4, 5,11,17,16,15,10,-1,
			 5, 6, 7,12,19,18,17,11,-1,
			 1, 9,13,20,19,12, 7, 8,-1,
			13,14,15,16,17,18,19,20,-1,
			 1, 8, 7, 6, 5, 4, 3, 2,-1 };

  FFlHEX20ElementTopSpec::instance()->setTopology(6,faces);
}


bool FFlHEX20::getFaceNormals(std::vector<FaVec3>& normals, short int face,
			      bool switchNormal) const
{
  std::vector<FFlNode*> nodes;
  if (!this->getFaceNodes(nodes,face,switchNormal)) return false;

  return FFlCurvedFace::faceNormals(nodes,normals);
}


bool FFlHEX20::getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const
{
  //TODO: Account for possibly curved egdges
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(3)->getPos());
  FaVec3 v3(this->getNode(5)->getPos());
  FaVec3 v4(this->getNode(7)->getPos());
  FaVec3 v5(this->getNode(13)->getPos());
  FaVec3 v6(this->getNode(15)->getPos());
  FaVec3 v7(this->getNode(17)->getPos());
  FaVec3 v8(this->getNode(19)->getPos());

  FFaVolume::hexVolume(v1,v2,v3,v4,v5,v6,v7,v8,volume);
  FFaVolume::hexCenter(v1,v2,v3,v4,v5,v6,v7,v8,cog);
  FFaVolume::hexMoment(v1,v2,v3,v4,v5,v6,v7,v8,inertia);

  return true;
}
