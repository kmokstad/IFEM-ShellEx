// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlQUAD4.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFaLib/FFaAlgebra/FFaMat33.H"
#include "FFaLib/FFaAlgebra/FFaVolume.H"
#include "FFaLib/FFaAlgebra/FFaMath.H"


void FFlQUAD4::init()
{
  FFlQUAD4TypeInfoSpec::instance()->setTypeName("QUAD4");
  FFlQUAD4TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SHELL_ELM);

  ElementFactory::instance()->registerCreator(FFlQUAD4TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlQUAD4::create,int,FFlElementBase*&) );

  FFlQUAD4AttributeSpec::instance()->addLegalAttribute("PTHICK",false);
  FFlQUAD4AttributeSpec::instance()->addLegalAttribute("PMAT",false);
  FFlQUAD4AttributeSpec::instance()->addLegalAttribute("PMATSHELL",false);
  FFlQUAD4AttributeSpec::instance()->addLegalAttribute("PCOMP",false);
  FFlQUAD4AttributeSpec::instance()->addLegalAttribute("PNSM",false);
  FFlQUAD4AttributeSpec::instance()->addLegalAttribute("PCOORDSYS",false);

  FFlQUAD4ElementTopSpec::instance()->setNodeCount(4);
  FFlQUAD4ElementTopSpec::instance()->setNodeDOFs(6);
  FFlQUAD4ElementTopSpec::instance()->setShellFaces(true);

  static int faces[] = { 1,2,3,4,-1 };
  FFlQUAD4ElementTopSpec::instance()->setTopology(1,faces);
}


FaMat33 FFlQUAD4::getGlobalizedElmCS() const
{
  //TODO,bh: use the possible PCOORDSYS here
  FaMat33 T;
  return T.makeGlobalizedCS(this->getNode(1)->getPos(),
			    this->getNode(2)->getPos(),
			    this->getNode(3)->getPos(),
			    this->getNode(4)->getPos());
}


bool FFlQUAD4::getFaceNormals(std::vector<FaVec3>& normals, short int,
			      bool switchNormal) const
{
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(2)->getPos());
  FaVec3 v3(this->getNode(3)->getPos());
  FaVec3 v4(this->getNode(3)->getPos());
  FaVec3 vn = (v3-v1) ^ (switchNormal ? v2-v4 : v4-v2);

  normals.clear();
  normals.push_back(vn.normalize());
  normals.push_back(vn);
  normals.push_back(vn);
  normals.push_back(vn);

  return true;
}


bool FFlQUAD4::getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const
{
  double th = this->getThickness();
  if (th <= 0.0) return false; // Should not happen

  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(2)->getPos());
  FaVec3 v3(this->getNode(3)->getPos());
  FaVec3 v4(this->getNode(4)->getPos());

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


#define N0(x) 0.5*(1.0-x)
#define N1(x) 0.5*(x+1.0)

FaVec3 FFlQUAD4::interpolate(const double* Xi,
                             const std::vector<FaVec3>& v) const
{
  return N0(Xi[0])*N0(Xi[1]) * v[0]
    +    N1(Xi[0])*N0(Xi[1]) * v[1]
    +    N1(Xi[0])*N1(Xi[1]) * v[2]
    +    N0(Xi[0])*N1(Xi[1]) * v[3];
}

FaVec3 FFlQUAD4::mapping(double xi, double eta, double) const
{
  return N0(xi)*N0(eta) * this->getNode(1)->getPos()
    +    N1(xi)*N0(eta) * this->getNode(2)->getPos()
    +    N1(xi)*N1(eta) * this->getNode(3)->getPos()
    +    N0(xi)*N1(eta) * this->getNode(4)->getPos();
}


bool FFlQUAD4::invertMapping(const FaVec3& X, double* Xi) const
{
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(2)->getPos());
  FaVec3 v3(this->getNode(3)->getPos());
  FaVec3 v4(this->getNode(4)->getPos());
  FaMat33 T = T.makeGlobalizedCS(v1,v2,v3,v4).transpose();
  FaVec3 Xp = T*X;
  FaVec3 X1 = T*v1;
  FaVec3 X2 = T*v2;
  FaVec3 X3 = T*v3;
  FaVec3 X4 = T*v4;

  // Initialize auxilliary variables
  double A[2][4];
  for (int j = 0; j < 2; j++)
  {
    A[j][0] =         0.25*( X1[j] - X2[j] + X3[j] - X4[j]);
    A[j][1] =         0.25*(-X1[j] + X2[j] + X3[j] - X4[j]);
    A[j][2] =         0.25*(-X1[j] - X2[j] + X3[j] + X4[j]);
    A[j][3] = Xp[j] - 0.25*( X1[j] + X2[j] + X3[j] + X4[j]);
  }

  /* We have to solve the following set of equations (j=0,1):

     Aj0*XI*ETA + Aj1*XI + Aj2*ETA = Aj3

     The way that we may solve this nonlinear (in XI and ETA)
     set of equations depends on the coefficients Aij.
     The solution is unique for proper input.
     See FFa::bilinearSolve for details.
  */

  double S1[4], S2[4];
  int NSOL = FFa::bilinearSolve(A[0],A[1],S1,S2);
  if (NSOL < 1)
  {
    std::cerr <<" *** FFlQUAD4::invertMapping: Failure for element "
              << this->getID() <<", NSOL = "<< NSOL;
    std::cerr <<"\n     A0 =";
    for (int i = 0; i < 4; i++) std::cerr <<" "<< A[0][i];
    std::cerr <<"\n     A1 =";
    for (int j = 0; j < 4; j++) std::cerr <<" "<< A[1][j];
    std::cerr <<"\n     X1 = "<< X1;
    std::cerr <<"\n     X2 = "<< X2;
    std::cerr <<"\n     X3 = "<< X3;
    std::cerr <<"\n     X4 = "<< X4;
    std::cerr <<"\n     Xp = "<< Xp << std::endl;
    return false;
  }

  const double epsO = 0.001;
  const double epsZ = 1.0e-8;
  const double tolO = 1.0+epsO;
  const double tolM = 2.0*epsZ*tolO;

  // Check that the solution(s) is(are) "inside" the element
  int numSol = 0;
  for (int i = 0; i < NSOL; i++)
    if (fabs(S1[i]) < tolO && fabs(S2[i]) < tolO)
    {
      if (numSol == 0)
      {
        numSol = 1;
        Xi[0] = S1[i];
        Xi[1] = S2[i];
      }
      else if (fabs(Xi[0]-S1[i]) < tolM && fabs(Xi[1]-S2[i]) < tolM)
      {
        // Multiple solutions, but they are "almost equal", use the average
        Xi[0] = 0.5*(Xi[0]+S1[i]);
        Xi[1] = 0.5*(Xi[1]+S2[i]);
      }
      else
      {
        ++numSol; // Choose the solution that is closest to the point X
        double d1 = (X - this->mapping(Xi[0],Xi[1])).sqrLength();
        double d2 = (X - this->mapping(S1[i],S2[i])).sqrLength();
        if (d2 < d1)
        {
          Xi[0] = S1[i];
          Xi[1] = S2[i];
        }
      }
    }

  if (numSol < 1)
  {
#if FFL_DEBUG > 1
    std::cout <<"  ** FFlQUAD4::invertMapping: Point "<< X
              <<" is not inside element "<< this->getID() << std::endl;
#endif
    return false;
  }

  double elArea = 0.5*(((v2-v1)^(v3-v1)).length() + ((v3-v1)^(v4-v1)).length());
  double elDia2 = elArea*4.0/M_PI;
  double distXp = Xp.z() - 0.25*(X1.z()+X2.z()+X3.z()+X4.z());
  if (distXp*distXp > offPlaneTol*offPlaneTol*elDia2)
  {
    // The point should be closer to the element surface than
    // (offPlaneTol*100)% of the equivalent element diameter
#ifdef FFL_DEBUG
    std::cout <<"  ** FFlQUAD::invertMapping: Point "<< X
              <<" is inside element "<< this->getID()
              <<", Xi1 = "<< Xi[0] <<", Xi2 = "<< Xi[1]
              <<"\n     but is too far ("<< distXp
              <<") from the element plane (D="<< sqrt(elDia2)
              <<")."<< std::endl;
#endif
    return false;
  }
#if FFL_DEBUG > 1
  else if (numSol > 1)
    std::cout <<"  ** FFlQUAD4::invertMapping: Point "<< X
              <<", multiple solutions for element "<< this->getID()
              << std::endl;
#endif

  return true;
}
