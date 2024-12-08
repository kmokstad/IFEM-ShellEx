// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlSolids.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlFEParts/FFlCurvedFace.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include "FFaLib/FFaAlgebra/FFaVolume.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


void FFlTET4::init()
{
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlTET4>;
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlTET4>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlTET4>;

  TypeInfoSpec::instance()->setTypeName("TET4");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SOLID_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlTET4::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PMAT");

  ElementTopSpec::instance()->setNodeCount(4);
  ElementTopSpec::instance()->setNodeDOFs(3);

  int faces[] = { 1, 3, 2,-1,
                  1, 2, 4,-1,
                  2, 3, 4,-1,
                  1, 4, 3,-1 };
  ElementTopSpec::instance()->setTopology(4,faces);
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


void FFlWEDG6::init()
{
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlWEDG6>;
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlWEDG6>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlWEDG6>;

  TypeInfoSpec::instance()->setTypeName("WEDG6");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SOLID_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlWEDG6::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PMAT");

  ElementTopSpec::instance()->setNodeCount(6);
  ElementTopSpec::instance()->setNodeDOFs(3);

  int faces[] = { 1, 2, 5, 4,-1,
                  2, 3, 6, 5,-1,
                  1, 4, 6, 3,-1,
                  1, 3, 2,   -1,
                  4, 5, 6,   -1 };
  ElementTopSpec::instance()->setTopology(5,faces);
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


void FFlHEX8::init()
{
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlHEX8>;
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlHEX8>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlHEX8>;

  TypeInfoSpec::instance()->setTypeName("HEX8");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SOLID_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlHEX8::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PMAT");

  ElementTopSpec::instance()->setNodeCount(8);
  ElementTopSpec::instance()->setNodeDOFs(3);

  int faces[] = { 2, 3, 7, 6,-1,
                  3, 4, 8, 7,-1,
                  1, 5, 8, 4,-1,
                  1, 2, 6, 5,-1,
                  5, 6, 7, 8,-1,
                  1, 4, 3, 2,-1 };
  ElementTopSpec::instance()->setTopology(6,faces);
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


void FFlTET10::init()
{
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlTET10>;
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlTET10>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlTET10>;

  TypeInfoSpec::instance()->setTypeName("TET10");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SOLID_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlTET10::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PMAT");

  ElementTopSpec::instance()->setNodeCount(10);
  ElementTopSpec::instance()->setNodeDOFs(3);

  int faces[] = { 1, 6, 5, 4, 3, 2,-1,
                  3, 4, 5, 9,10, 8,-1,
                  1, 7,10, 9, 5, 6,-1,
                  1, 2, 3, 8,10, 7,-1 };
  ElementTopSpec::instance()->setTopology(4,faces);
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


void FFlWEDG15::init()
{
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlWEDG15>;
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlWEDG15>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlWEDG15>;

  TypeInfoSpec::instance()->setTypeName("WEDG15");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SOLID_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlWEDG15::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PMAT");

  ElementTopSpec::instance()->setNodeCount(15);
  ElementTopSpec::instance()->setNodeDOFs(3);

  int faces[] = { 1, 2, 3, 8,12,11,10, 7,-1,
                  3, 4, 5, 9,14,13,12, 8,-1,
                  1, 7,10,15,14, 9, 5, 6,-1,
                  1, 6, 5, 4, 3, 2,      -1,
                 10,11,12,13,14,15,      -1 };
  ElementTopSpec::instance()->setTopology(5,faces);
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


void FFlHEX20::init()
{
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlHEX20>;
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlHEX20>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlHEX20>;

  TypeInfoSpec::instance()->setTypeName("HEX20");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SOLID_ELM);
  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlHEX20::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PMAT");

  ElementTopSpec::instance()->setNodeCount(20);
  ElementTopSpec::instance()->setNodeDOFs(3);

  int faces[] = { 1, 2, 3,10,15,14,13, 9,-1,
                  3, 4, 5,11,17,16,15,10,-1,
                  5, 6, 7,12,19,18,17,11,-1,
                  1, 9,13,20,19,12, 7, 8,-1,
                 13,14,15,16,17,18,19,20,-1,
                  1, 8, 7, 6, 5, 4, 3, 2,-1 };
  ElementTopSpec::instance()->setTopology(6,faces);
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

#ifdef FF_NAMESPACE
} // namespace
#endif
