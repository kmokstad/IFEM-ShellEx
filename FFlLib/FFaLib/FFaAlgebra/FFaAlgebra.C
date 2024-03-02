// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFaLib/FFaAlgebra/FFaMat33.H"
#include "FFaLib/FFaAlgebra/FFaAlgebra.H"


/*!
  Performs an eccentricity transformation of a 6x6 element matrix.
*/

void FFaAlgebra::eccTransform6 (double mat[6][6], const FaVec3& X)
{
  // The following code is based on the fortran routine TRIX30
  // from the file vpmBase/Femlib/beamaux.f
  // Note: The vector X is here assumed to point FROM the nodal point location
  // TO the actual element location (this is opposite to the case in beamaux.f)

  int i;

  for (i = 0; i < 6; i++)
  {
    mat[3][i] -= X[2]*mat[1][i] - X[1]*mat[2][i];
    mat[4][i] -= X[0]*mat[2][i] - X[2]*mat[0][i];
    mat[5][i] -= X[1]*mat[0][i] - X[0]*mat[1][i];
  }

  for (int i = 0; i < 6; i++)
  {
    mat[i][3] -= X[2]*mat[i][1] - X[1]*mat[i][2];
    mat[i][4] -= X[0]*mat[i][2] - X[2]*mat[i][0];
    mat[i][5] -= X[1]*mat[i][0] - X[0]*mat[i][1];
  }
}


/*!
  Performs a congruence transformation of a symmetric (3*N)x(3*N) matrix.
  The transformation matrix consists of a 3*3 submatrix T, which is
  repeated along the diagonal (when node == 0).
  If node > 0, the transformation matrix equals the identity matrix, but with
  the submatrix T inserted on the diagonal at position 3*(node-1)+1 to 3*node.
*/

bool FFaAlgebra::congruenceTransform (double** mat, const FaMat33& T,
				      int N, int node)
{
  if (N < 1 || node > N) return false;

  // The following code is based on the fortran routine MPRO30
  // from the file vpmBase/Femlib/beamaux.f

  int i, j, iA, jA, kA;
  double tmp, B[3][3];

  for (i = 0; i < N; i++)
    for (j = i; j < N; j++)
    {
      for (iA = 0; iA < 3; iA++)
        for (jA = 0; jA < 3; jA++)
          B[iA][jA] = 0.0;

      if (node < 1 || i == node-1)
	for (iA = 0; iA < 3; iA++)
	  for (kA = 0; kA < 3; kA++)
	  {
	    tmp = mat[3*j+kA][3*i+iA];
	    for (jA = 0; jA < 3; jA++)
	      B[jA][iA] += T[kA][jA] * tmp;
	  }

      else if (j == node-1)
	for (iA = 0; iA < 3; iA++)
	  for (jA = 0; jA < 3; jA++)
	    B[jA][iA] = mat[3*j+jA][3*i+iA];

      if (node < 1 || j == node-1)
	for (iA = 0; iA < 3; iA++)
	  for (kA = 0; kA < 3; kA++)
	  {
	    tmp = 0.0;
	    for (jA = 0; jA < 3; jA++)
	      tmp += B[kA][jA] * T[jA][iA];
	    mat[3*j+kA][3*i+iA] = tmp;
	  }

      else if (i == node-1)
	for (iA = 0; iA < 3; iA++)
	  for (jA = 0; jA < 3; jA++)
	    mat[3*j+jA][3*i+iA] = B[jA][iA];
    }

  for (i = 0; i < 3*N; i++)
    for (j = i+1; j < 3*N; j++)
      mat[i][j] = mat[j][i];

  return true;
}
