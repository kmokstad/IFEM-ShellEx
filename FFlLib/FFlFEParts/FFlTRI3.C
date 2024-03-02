// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlTRI3.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFaLib/FFaAlgebra/FFaMat33.H"
#include "FFaLib/FFaAlgebra/FFaVolume.H"
#include "FFaLib/FFaAlgebra/FFaMath.H"


void FFlTRI3::init()
{
  FFlTRI3TypeInfoSpec::instance()->setTypeName("TRI3");
  FFlTRI3TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SHELL_ELM);

  ElementFactory::instance()->registerCreator(FFlTRI3TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlTRI3::create,int,FFlElementBase*&));

  FFlTRI3AttributeSpec::instance()->addLegalAttribute("PTHICK",false);
  FFlTRI3AttributeSpec::instance()->addLegalAttribute("PMAT",false);
  FFlTRI3AttributeSpec::instance()->addLegalAttribute("PMATSHELL",false);
  FFlTRI3AttributeSpec::instance()->addLegalAttribute("PCOMP",false);
  FFlTRI3AttributeSpec::instance()->addLegalAttribute("PNSM",false);
  FFlTRI3AttributeSpec::instance()->addLegalAttribute("PCOORDSYS",false);

  FFlTRI3ElementTopSpec::instance()->setNodeCount(3);
  FFlTRI3ElementTopSpec::instance()->setNodeDOFs(6);
  FFlTRI3ElementTopSpec::instance()->setShellFaces(true);

  static int faces[] = { 1,2,3,-1 };
  FFlTRI3ElementTopSpec::instance()->setTopology(1,faces);
}


FaMat33 FFlTRI3::getGlobalizedElmCS() const
{
  //TODO,bh: use the possible PCOORDSYS here
  FaMat33 T;
  return T.makeGlobalizedCS(this->getNode(1)->getPos(),
			    this->getNode(2)->getPos(),
			    this->getNode(3)->getPos());
}


bool FFlTRI3::getFaceNormals(std::vector<FaVec3>& normals, short int,
			     bool switchNormal) const
{
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(2)->getPos());
  FaVec3 v3(this->getNode(3)->getPos());
  FaVec3 vn = (v2-v1) ^ (switchNormal ? v1-v3 : v3-v1);

  normals.clear();
  normals.push_back(vn.normalize());
  normals.push_back(vn);
  normals.push_back(vn);

  return true;
}


bool FFlTRI3::getVolumeAndInertia(double& volume, FaVec3& cog,
				  FFaTensor3& inertia) const
{
  double th = this->getThickness();
  if (th <= 0.0) return false; // Should not happen

  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(2)->getPos());
  FaVec3 v3(this->getNode(3)->getPos());

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


FaVec3 FFlTRI3::interpolate(const double* Xi,
                            const std::vector<FaVec3>& v) const
{
  return Xi[0]*v[0] + Xi[1]*v[1] + (1.0-Xi[0]-Xi[1])*v[2];
}

FaVec3 FFlTRI3::mapping(double xi, double eta, double) const
{
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(2)->getPos());
  FaVec3 v3(this->getNode(3)->getPos());

  return xi*v1 + eta*v2 + (1.0-xi-eta)*v3;
}


bool FFlTRI3::invertMapping(const FaVec3& X, double* Xi) const
{
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(2)->getPos());
  FaVec3 v3(this->getNode(3)->getPos());

  FaVec3 normal = (v2-v1) ^ (v3-v1);
  double elArea = normal.length();
  if (elArea < 1.0e-16)
  {
    std::cout <<" *** FFlTRI3::invertMapping: Degenerated element "
              << this->getID() << std::endl;
    return false;
  }

  normal /= elArea;
  FaVec3 aux = X - v1;
  double dis = normal*aux;
  FaVec3 v1X = aux - dis*normal;

  aux = (v2-v1) ^ v1X;
  double A3 = copysign(aux.length(),aux*normal);
  aux = v1X ^ (v3-v1);
  double A2 = copysign(aux.length(),aux*normal);

  Xi[0] = (elArea-A2-A3) / elArea;
  Xi[1] = A2 / elArea;

  const double epsO = 0.001;
  const double tolO = 1.0 + epsO;

  if (Xi[0] < -epsO || Xi[1] > tolO ||
      Xi[1] < -epsO || Xi[0] > tolO ||
      Xi[0]+Xi[1] < -epsO || Xi[0]+Xi[1] > tolO)
  {
#if FFL_DEBUG > 1
    std::cout <<"  ** FFlTRI3::invertMapping: Point "<< X
              <<" is not inside element "<< this->getID()
              <<", Xi1 = "<< Xi[0] <<", Xi2 = "<< Xi[1] << std::endl;
#endif
    return false;
  }

  double elDia2 = elArea*2.0/M_PI;
  if (dis*dis > offPlaneTol*offPlaneTol*elDia2)
  {
    // The point should be closer to the element surface than
    // (offPlaneTol*100)% of the equivalent element diameter
#ifdef FFL_DEBUG
    std::cout <<"  ** FFlTRI3::invertMapping: Point "<< X
              <<" is inside element "<< this->getID()
              <<", Xi1 = "<< Xi[0] <<", Xi2 = "<< Xi[1]
              <<"\n     but is too far ("<< dis
              <<") from the element plane (D="<< sqrt(elDia2)
              <<")."<< std::endl;
#endif
    return false;
  }

  return true;
}
