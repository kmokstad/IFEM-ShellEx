// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlShells.H"
#include "FFlLib/FFlFEParts/FFlPTHICK.H"
#include "FFlLib/FFlFEParts/FFlPCOMP.H"
#include "FFlLib/FFlFEParts/FFlPMAT.H"
#include "FFlLib/FFlFEParts/FFlPNSM.H"
#include "FFlLib/FFlFEParts/FFlCurvedFace.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlLinkHandler.H"
#include "FFaLib/FFaAlgebra/FFaMat33.H"
#include "FFaLib/FFaAlgebra/FFaVolume.H"
#include "FFaLib/FFaAlgebra/FFaMath.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif

double FFlShellElementBase::offPlaneTol = 0.1;


double FFlShellElementBase::getThickness() const
{
  FFlPTHICK* pthk = dynamic_cast<FFlPTHICK*>(this->getAttribute("PTHICK"));
  if (pthk)
    return pthk->thickness.getValue();
  FFlPCOMP* pcomp = dynamic_cast<FFlPCOMP*>(this->getAttribute("PCOMP"));
  if (pcomp)
    return -2.0*pcomp->Z0.getValue(); //TODO,bh: This is quick and dirty
  return 0.0; // Should not happen
}


double FFlShellElementBase::getMassDensity() const
{
  FFlPMAT* pmat = dynamic_cast<FFlPMAT*>(this->getAttribute("PMAT"));
  double rho = pmat ? pmat->materialDensity.getValue() : 0.0;

  // Check for non-structural mass
  FFlPNSM* pnsm = dynamic_cast<FFlPNSM*>(this->getAttribute("PNSM"));
  if (!pnsm) return rho;

  // Shell thickness
  double th = this->getThickness();
  if (th < 1.0e-16) return rho; // avoid division by zero

  // Modify the mass density to account for the non-structural mass
  return rho + pnsm->NSM.getValue() / th;
}


void FFlTRI3::init()
{
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlTRI3>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlTRI3>;
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlTRI3>;

  TypeInfoSpec::instance()->setTypeName("TRI3");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SHELL_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlTRI3::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PTHICK",false);
  AttributeSpec::instance()->addLegalAttribute("PMAT",false);
  AttributeSpec::instance()->addLegalAttribute("PMATSHELL",false);
  AttributeSpec::instance()->addLegalAttribute("PCOMP",false);
  AttributeSpec::instance()->addLegalAttribute("PNSM",false);
  AttributeSpec::instance()->addLegalAttribute("PCOORDSYS",false);

  ElementTopSpec::instance()->setNodeCount(3);
  ElementTopSpec::instance()->setNodeDOFs(6);
  ElementTopSpec::instance()->setShellFaces(true);

  static int faces[] = { 1,2,3,-1 };
  ElementTopSpec::instance()->setTopology(1,faces);
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


void FFlTRI6::init()
{
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlTRI6>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlTRI6>;
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlTRI6>;

  TypeInfoSpec::instance()->setTypeName("TRI6");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SHELL_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlTRI6::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PTHICK",false);
  AttributeSpec::instance()->addLegalAttribute("PMAT",false);
  AttributeSpec::instance()->addLegalAttribute("PMATSHELL",false);
  AttributeSpec::instance()->addLegalAttribute("PCOMP",false);
  AttributeSpec::instance()->addLegalAttribute("PNSM",false);
  AttributeSpec::instance()->addLegalAttribute("PCOORDSYS",false);

  ElementTopSpec::instance()->setNodeCount(6);
  ElementTopSpec::instance()->setNodeDOFs(6);
  ElementTopSpec::instance()->setShellFaces(true);

  static int faces[] = { 1,2,3,4,5,6,-1 };
  ElementTopSpec::instance()->setTopology(1,faces);
}


bool FFlTRI6::split(Elements& newElem, FFlLinkHandler* owner, int)
{
  newElem.clear();
  newElem.reserve(4);

  int elm_id = owner->getNewElmID();
  for (int i = 0; i < 4; i++)
    newElem.push_back(ElementFactory::instance()->create("TRI3",elm_id++));

  newElem[0]->setNode(1,this->getNodeID(1));
  newElem[0]->setNode(2,this->getNodeID(2));
  newElem[0]->setNode(3,this->getNodeID(6));
  newElem[1]->setNode(1,this->getNodeID(2));
  newElem[1]->setNode(2,this->getNodeID(3));
  newElem[1]->setNode(3,this->getNodeID(4));
  newElem[2]->setNode(1,this->getNodeID(6));
  newElem[2]->setNode(2,this->getNodeID(4));
  newElem[2]->setNode(3,this->getNodeID(5));
  newElem[3]->setNode(1,this->getNodeID(4));
  newElem[3]->setNode(2,this->getNodeID(6));
  newElem[3]->setNode(3,this->getNodeID(2));

  return true;
}


FaMat33 FFlTRI6::getGlobalizedElmCS() const
{
  //TODO,bh: use the possible PCOORDSYS here
  FaMat33 T;
  return T.makeGlobalizedCS(this->getNode(1)->getPos(),
                            this->getNode(3)->getPos(),
                            this->getNode(5)->getPos());
}


bool FFlTRI6::getFaceNormals(std::vector<FaVec3>& normals, short int,
                             bool switchNormal) const
{
  std::vector<FFlNode*> nodes;
  if (!this->getFaceNodes(nodes,1,switchNormal)) return false;

  return FFlCurvedFace::faceNormals(nodes,normals);
}


bool FFlTRI6::getVolumeAndInertia(double& volume, FaVec3& cog,
                                  FFaTensor3& inertia) const
{
  double th = this->getThickness();
  if (th <= 0.0) return false; // Should not happen

  //TODO,kmo: Account for curved edges and face (by numerical integration)
  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(3)->getPos());
  FaVec3 v3(this->getNode(5)->getPos());

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


int FFlTRI6::getNodalCoor(double* X, double* Y, double* Z) const
{
  int ierr = this->FFlElementBase::getNodalCoor(X,Y,Z);
  if (ierr < 0) return ierr;

  // Reorder the nodes such that the 3 mid-side nodes are ordered last
  // 1-2-3-4-5-6 --> 1-3-5-2-4-6
  std::swap(X[1],X[2]); std::swap(X[2],X[4]); std::swap(X[3],X[4]);
  std::swap(Y[1],Y[2]); std::swap(Y[2],Y[4]); std::swap(Y[3],Y[4]);
  std::swap(Z[1],Z[2]); std::swap(Z[2],Z[4]); std::swap(Z[3],Z[4]);

  return ierr;
}


void FFlQUAD4::init()
{
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlQUAD4>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlQUAD4>;
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlQUAD4>;

  TypeInfoSpec::instance()->setTypeName("QUAD4");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SHELL_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlQUAD4::create,int,FFlElementBase*&) );

  AttributeSpec::instance()->addLegalAttribute("PTHICK",false);
  AttributeSpec::instance()->addLegalAttribute("PMAT",false);
  AttributeSpec::instance()->addLegalAttribute("PMATSHELL",false);
  AttributeSpec::instance()->addLegalAttribute("PCOMP",false);
  AttributeSpec::instance()->addLegalAttribute("PNSM",false);
  AttributeSpec::instance()->addLegalAttribute("PCOORDSYS",false);

  ElementTopSpec::instance()->setNodeCount(4);
  ElementTopSpec::instance()->setNodeDOFs(6);
  ElementTopSpec::instance()->setShellFaces(true);

  static int faces[] = { 1,2,3,4,-1 };
  ElementTopSpec::instance()->setTopology(1,faces);
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


void FFlQUAD8::init()
{
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlQUAD8>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlQUAD8>;
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlQUAD8>;

  TypeInfoSpec::instance()->setTypeName("QUAD8");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::SHELL_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlQUAD8::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PTHICK",false);
  AttributeSpec::instance()->addLegalAttribute("PMAT",false);
  AttributeSpec::instance()->addLegalAttribute("PMATSHELL",false);
  AttributeSpec::instance()->addLegalAttribute("PCOMP",false);
  AttributeSpec::instance()->addLegalAttribute("PNSM",false);
  AttributeSpec::instance()->addLegalAttribute("PCOORDSYS",false);

  ElementTopSpec::instance()->setNodeCount(8);
  ElementTopSpec::instance()->setNodeDOFs(6);
  ElementTopSpec::instance()->setShellFaces(true);

  static int faces[] = { 1,2,3,4,5,6,7,8,-1 };
  ElementTopSpec::instance()->setTopology(1,faces);
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

#ifdef FF_NAMESPACE
} // namespace
#endif
