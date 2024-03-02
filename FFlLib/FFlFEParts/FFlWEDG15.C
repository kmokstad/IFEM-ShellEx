// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlWEDG15.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlFEParts/FFlCurvedFace.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include "FFaLib/FFaAlgebra/FFaVolume.H"


void FFlWEDG15::init()
{
  FFlWEDG15TypeInfoSpec::instance()->setTypeName("WEDG15");
  FFlWEDG15TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SOLID_ELM);

  ElementFactory::instance()->registerCreator(FFlWEDG15TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlWEDG15::create,int,FFlElementBase*&));

  FFlWEDG15AttributeSpec::instance()->addLegalAttribute("PMAT");

  FFlWEDG15ElementTopSpec::instance()->setNodeCount(15);
  FFlWEDG15ElementTopSpec::instance()->setNodeDOFs(3);

  static int faces[] = { 1, 2, 3, 8,12,11,10, 7,-1,
			 3, 4, 5, 9,14,13,12, 8,-1,
			 1, 7,10,15,14, 9, 5, 6,-1,
			 1, 6, 5, 4, 3, 2,      -1,
			10,11,12,13,14,15,      -1 };

  FFlWEDG15ElementTopSpec::instance()->setTopology(5,faces);
}


bool FFlWEDG15::getFaceNormals(std::vector<FaVec3>& normals, short int face,
			       bool switchNormal) const
{
  std::vector<FFlNode*> nodes;
  if (!this->getFaceNodes(nodes,face,switchNormal)) return false;

  return FFlCurvedFace::faceNormals(nodes,normals);
}


bool FFlWEDG15::getVolumeAndInertia(double& volume, FaVec3& cog,
				    FFaTensor3& inertia) const
{
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(3)->getPos());
  FaVec3 v3(this->getNode(5)->getPos());
  FaVec3 v4(this->getNode(10)->getPos());
  FaVec3 v5(this->getNode(12)->getPos());
  FaVec3 v6(this->getNode(14)->getPos());

  FFaVolume::wedVolume(v1,v2,v3,v4,v5,v6,volume);
  FFaVolume::wedCenter(v1,v2,v3,v4,v5,v6,cog);
  FFaVolume::wedMoment(v1,v2,v3,v4,v5,v6,inertia);

  return true;
}
