// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFaLib/FFaAlgebra/FFaTensorTransforms.H"
#include "FFaLib/FFaAlgebra/FFaTensor3.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include "FFaLib/FFaOS/FFaFortran.H"

/*!
  von Mises calculation for stress and strain tensors.
*/

DOUBLE_FUNCTION(vonmises,VONMISES) (const int& N, double* S)
{
  return FFaTensorTransforms::vonMises(N,S);
}


/*!
  Principal values calculation.
*/

SUBROUTINE(princval,PRINCVAL) (const int& N, double* S, double* Pv)
{
  FFaTensorTransforms::principalValues(N,S,Pv);
}


/*!
  Maximum shear value calculation.
*/

SUBROUTINE(maxshearvalue,MAXSHEARVALUE) (const int& N, double* Pv, double& S)
{
  S = FFaTensorTransforms::maxShearValue(Pv[0],Pv[N-1]);
}


/*!
  Maximum shear value and associated direction.
*/

SUBROUTINE(maxshear,MAXSHEAR) (const int& N, double* Pv, double* Pd,
                               double& S, double* Sd)
{
  S = FFaTensorTransforms::maxShearValue(Pv[0],Pv[N-1]);
  FFaTensorTransforms::maxShearDir(N,Pd,Pd+N*N-N,Sd);
}


/*!
  Congruence transformation of 2D and 3D symmetric tensors.
*/

SUBROUTINE(tratensor,TRATENSOR) (const int& N, double* S, const double* T)
{
  if (N == 2)
    FFaTensorTransforms::rotate2D(S,T,S);
  else if (N == 3)
    FFaTensorTransforms::rotate3D(S,T,S);
}


/*!
  Inertia tensor transformation based on the parallel-axis theorem.
*/

SUBROUTINE(trainertia,TRAINERTIA) (const int& N, double* inertia,
                                   const double* x, const double& mass)
{
  if (N == 3)
  {
    FFaTensor3 I(inertia);
    FaVec3 X(x[0],x[1],x[2]);
    I.translateInertia(X,mass);
    for (int i = 0; i < 6; i++) inertia[i] = I[i];
  }
}
