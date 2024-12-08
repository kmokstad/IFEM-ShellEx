// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include <cstring>

#include "FFlLib/FFlFEParts/FFlCMASS.H"
#include "FFlLib/FFlFEParts/FFlPMASS.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFaLib/FFaAlgebra/FFaTensor3.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


void FFlCMASS::init()
{
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlCMASS>;
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlCMASS>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlCMASS>;

  TypeInfoSpec::instance()->setTypeName("CMASS");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::OTHER_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlCMASS::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PMASS",false);

  ElementTopSpec::instance()->setNodeCount(1);
  ElementTopSpec::instance()->setNodeDOFs(3);
}


bool FFlCMASS::getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const
{
  volume = 0.0;
  cog = this->getNode(1)->getPos();

  FFlPMASS* pmass = dynamic_cast<FFlPMASS*>(this->getAttribute("PMASS"));
  if (!pmass)
  {
    // Property-less mass element, ignore
    inertia = FFaTensor3(0.0);
    return true;
  }

  const std::vector<double>& M = pmass->M.getValue();
  if (M.size() > 0 && M.size() <= 6)
  {
    // No inertia for this element, only mass
    volume = M[0];
    inertia = FFaTensor3(0.0);
    return true;
  }
  else if (M.size() != 21)
  {
    ListUI <<"\n  ** Warning: Mass matrix for CMASS element "<< this->getID()
           <<" has invalid size "<< M.size() <<"\n";
    return false;
  }

  // Compute the inertia tensor about origo
  double R[6][6], II[6][6];
  memset(R,0,36*sizeof(double));
  R[0][0] = R[1][1] = R[2][2] = 1.0;
  R[1][3] = -cog[2];
  R[2][3] =  cog[1];
  R[3][3] = 1.0;
  R[0][4] =  cog[2];
  R[2][4] = -cog[0];
  R[4][4] = 1.0;
  R[0][5] = -cog[1];
  R[1][5] =  cog[0];
  R[5][5] = 1.0;
  bool haveInertia = RtMR(M,R,II);

  // Compute the mass = volume
  volume = (II[0][0]+II[1][1]+II[2][2]) / 3.0;
  if (!haveInertia)
  {
    inertia = FFaTensor3(0.0);
    return true;
  }
  else if (volume < 1.0e-16)
  {
    ListUI <<"\n  ** Warning: CMASS element "<< this->getID()
           <<" has zero mass, but non-zero inertia. This is non-physical.\n";
    return false;
  }

  // Adjust the center of gravity for possible offset
  cog = FaVec3(II[1][5],II[2][3],II[0][4]) / volume;

  // Compute the inertia tensor about the center of gravity
  R[1][3] += cog[2];
  R[2][3] -= cog[1];
  R[0][4] -= cog[2];
  R[2][4] += cog[0];
  R[0][5] += cog[1];
  R[1][5] -= cog[0];
  RtMR(M,R,II);

  inertia = FFaTensor3(II[3][3],II[4][4],II[5][5],
                       II[3][4],II[3][5],II[4][5]);
  return true;
}


/*!
  Evaluate the triple matrix-product  II = R^t * M * R
*/

bool FFlCMASS::RtMR(const std::vector<double>& M,
                    double R[6][6], double II[6][6])
{
  bool haveNonZeroInertia = false;
  unsigned int i, j, k = 0;

  // Set up the full quadratic mass matrix in II
  // Check if there actually are inertia terms in it
  for (j = 0; j < 6; j++)
  {
    for (i = 0; i < j; i++)
      II[i][j] = II[j][i] = M[k++];
    II[j][j] = M[k++];
    for (i = 0; i <= j && !haveNonZeroInertia; i++)
      if (i < j || j > 2)
        if (II[i][j] != 0.0)
          haveNonZeroInertia = true;
  }

  // Multiply  II * R ==> II
  double tmp[6];
  for (i = 0; i < 6; i++)
  {
    for (j = 0; j < 6; j++)
    {
      tmp[j] = 0.0;
      for (k = 0; k < 6; k++)
        tmp[j] += II[i][k]*R[k][j];
    }
    for (j = 0; j < 6; j++)
      II[i][j] = tmp[j];
  }

  // Multiply  R^t * II ==> II
  for (j = 0; j < 6; j++)
  {
    for (i = 0; i < 6; i++)
    {
      tmp[i] = 0.0;
      for (k = 0; k < 6; k++)
        tmp[i] += R[k][i]*II[k][j];
    }
    for (i = 0; i < 6; i++)
      II[i][j] = tmp[i];
  }

  return haveNonZeroInertia;
}

#ifdef FF_NAMESPACE
} // namespace
#endif
