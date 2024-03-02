// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlCurvedFace.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"


/*!
  This method computes the outward-pointing normal vector at each nodal point
  of a possibly curved second-order element face.
*/

bool FFlCurvedFace::faceNormals(const std::vector<FFlNode*>& nodes,
				std::vector<FaVec3>& normals)
{
  /*
    We need the shape function derivatives at the nodal points of the face
    to account for the curvature. The local ordering is assumed as follows:

    |              4                      (5)
    |              *            6 *--------*--------* 4
    |             /|              |                 |
    |            / |              |        X2       |
    |           /  |              |        |        |
    |      (5) *   * (3)      (7) *        o---X1   * (3)
    |         /    |              |                 |
    |        /     |              |                 |
    |       /      |              |                 |
    |    0 *---*---* 2          0 *--------*--------* 2
    |         (1)                         (1)

  */

  //                       0    1    2    3    4    5    6    7
  static double A1[6] = { 1.0, 0.5, 0.0, 0.0, 0.0, 0.5 };
  static double A2[6] = { 0.0, 0.5, 1.0, 0.5, 0.0, 0.0 };

  static double X1[8] = {-1.0, 0.0, 1.0, 1.0, 1.0, 0.0,-1.0,-1.0 };
  static double X2[8] = {-1.0,-1.0,-1.0, 0.0, 1.0, 1.0, 1.0, 0.0 };

  normals.clear();
  normals.reserve(nodes.size());

  double N1[8], N2[8];
  for (size_t i = 0; i < nodes.size(); i++)
  {
    // Find shape function derivatives at this node
    switch (nodes.size())
      {
      case 6: shapeDerivs6(N1,N2,A1[i],A2[i]); break;
      case 8: shapeDerivs8(N1,N2,X1[i],X2[i]); break;
      default: return false;
      }

    // Evaluate the two tangent vectors at this node
    FaVec3 v1, v2;
    for (size_t j = 0; j < nodes.size(); j++)
    {
      if (N1[j] != 0.0) v1 += N1[j]*nodes[j]->getPos();
      if (N2[j] != 0.0) v2 += N2[j]*nodes[j]->getPos();
    }

    // Cross product to get normal vector
    normals.push_back(v1^v2);
    normals.back().normalize();
  }

  return true;
}


/*!
  Shape function derivatives for a 6-noded triangle
*/

void FFlCurvedFace::shapeDerivs6(double* N1, double* N2, double X1, double X2)
{
  double X3 = 1.0 - X1 - X2;

  N1[0] =  4.0*X1 - 1.0;
  N1[1] =  4.0*X2;
  N1[2] =  0.0;
  N1[3] = -4.0*X2;
  N1[4] =  1.0 - 4.0*X3;
  N1[5] =  4.0*(X3-X1);

  N2[0] =  0.0;
  N2[1] =  4.0*X1;
  N2[2] =  4.0*X2 - 1.0;
  N2[3] =  4.0*(X3-X2);
  N2[4] =  1.0 - 4.0*X3;
  N2[5] = -4.0*X1;
}


/*!
  Shape function derivatives for a 8-noded quadrilateral
*/

void FFlCurvedFace::shapeDerivs8(double* N1, double* N2, double X1, double X2)
{
  N1[1] = -X1*(1.0-X2);
  N1[3] =  0.5 - 0.5*X2*X2;
  N1[5] = -X1*(1.0+X2);
  N1[7] = -0.5 + 0.5*X2*X2;

  N1[0] = -0.25*(1.0-X2) - 0.5*(N1[7]+N1[1]);
  N1[2] =  0.25*(1.0-X2) - 0.5*(N1[1]+N1[3]);
  N1[4] =  0.25*(1.0+X2) - 0.5*(N1[3]+N1[5]);
  N1[6] = -0.25*(1.0+X2) - 0.5*(N1[5]+N1[7]);

  N2[1] = -0.5 + 0.5*X1*X1;
  N2[3] = -X2*(1.0+X1);
  N2[5] =  0.5 - 0.5*X1*X1;
  N2[7] = -X2*(1.0-X1);

  N2[0] = -0.25*(1.0-X1) - 0.5*(N2[7]+N2[1]);
  N2[2] = -0.25*(1.0+X1) - 0.5*(N2[1]+N2[3]);
  N2[4] =  0.25*(1.0+X1) - 0.5*(N2[3]+N2[5]);
  N2[6] =  0.25*(1.0-X1) - 0.5*(N2[5]+N2[7]);
}
