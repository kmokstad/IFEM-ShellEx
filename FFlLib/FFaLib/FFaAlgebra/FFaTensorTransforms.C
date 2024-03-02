// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFaTensorTransforms.C
  \brief Tensor transformation utilities.
*/

#include <iostream>

#include "FFaLib/FFaAlgebra/FFaTensorTransforms.H"
#include "FFaLib/FFaAlgebra/FFaMath.H"

#ifdef FT_HAS_LAPACK
#include "FFaLib/FFaOS/FFaFortran.H"

SUBROUTINE(dsyev,DSYEV)
#ifdef _NCHAR_AFTER_CHARARG
  (const char*, const int, const char*, const int,
   const int&, double*, const int&, double*, double*, const int&, int&);
#else
  (const char*, const char*,
   const int&, double*, const int&, double*, double*, const int&, int&,
   const int, const int);
#endif
#endif


double FFaTensorTransforms::vonMises (double s11, double s22, double s12)
{
  return sqrt(s11*s11 + s22*s22 - s11*s22 + 3.0*s12*s12);
}

double FFaTensorTransforms::vonMises (double s11, double s22, double s33,
                                      double s12, double s13, double s23)
{
  return sqrt(s11*s11 + s22*s22 + s33*s33 - s11*s22 - s22*s33 - s33*s11
	      + 3.0*(s12*s12 + s13*s13 + s23*s23));
}

/*!
  \param[in] N Tensor dimension (1, 2 or 3)
  \param[in] S Tensor values, upper triangle of the symmetric tensor
  \details The tensor values are assumed laid out as follows
  - N = 1: S = {s11}
  - N = 2: S = {s11, s22, s12}
  - N = 3: S = {s11, s22, s33, s12, s13, s23}
*/

double FFaTensorTransforms::vonMises (int N, const double* S)
{
  switch (N)
    {
    case 1:
      return S[0];
    case 2:
      return vonMises(S[0], S[1], S[2]);
    case 3:
      return vonMises(S[0], S[1], S[2], S[3], S[4], S[5]);
    default:
      return -HUGE_VAL;
    }
}


/*!
  \param[in] N Tensor dimension (1, 2 or 3)
  \param[in] S Tensor values, upper triangle of the symmetric tensor
  \param[out] pVal Principal values in descending order, i.e., {P1,P2,P3}
  \param[out] pDir Associated principal direction vectors, e.g., in 3D:
              - dir1 = pDir[0:2]
              - dir2 = pDir[3:5]
              - dir3 = pDir[6:8]
  \return Negative value on error, otherwise zero

  \details The tensor values are assumed laid out as follows
  - N = 1: S = {s11}
  - N = 2: S = {s11, s22, s12}
  - N = 3: S = {s11, s22, s33, s12, s13, s23}

  This function uses the LAPACK eigenvalue solver DSYEV.
*/

int FFaTensorTransforms::principalDirs (int N, const double* S,
                                        double* pVal, double* pDir)
{
  if (N == 1)
  {
    pVal[0] = S[0];
    pDir[0] = 1.0;
    return 0;
  }

#ifdef FT_HAS_LAPACK

  // Set up the upper triangle of the tensor in pDir

  if (N == 2)
  {
    pDir[0] = S[0];
    pDir[3] = S[1];
    pDir[2] = S[2];
  }
  else if (N == 3)
  {
    pDir[0] = S[0];
    pDir[4] = S[1];
    pDir[8] = S[2];
    pDir[3] = S[3];
    pDir[6] = S[4];
    pDir[7] = S[5];
  }
  else
    return -1;

  // Solve the symmetric eigenvalue problem

  int       info = 0;
  const int Lwork = 12;
  double    Work[Lwork];
  F90_NAME(dsyev,DSYEV) ("V","U",N,pDir,N,pVal,Work,Lwork,info,1,1);
  if (info) return info;

  // DSYEV returns the eigenvalues in ascending order, but we
  // want them in descending order (largest principal value first)

  double tp = pVal[0];
  pVal[0]   = pVal[N-1];
  pVal[N-1] = tp;

  int i, j;
  for (i = 0; i < N; i++)
  {
    j = N*(N-1) + i;
    tp = pDir[i];
    pDir[i] = pDir[j];
    pDir[j] = tp;
  }

#else
  int i;
  for (i=0; i < N;   i++)    pVal[i] = HUGE_VAL;
  for (i=0; i < N*N; i++)    pDir[i] = 0.0;
  for (i=0; i < N*N; i+=N+1) pDir[i] = 1.0;
#endif
  return 0;
}

/*!
  \param[in] S Tensor values, {s11, s22, s33, s12, s13, s23}
  \param[out] pVal Principal values, {pMax,pMid,pMin}
  \param[out] p1Dir Direction vector for the first principal value
  \param[out] p2Dir Direction vector for the second principal value
  \param[out] p3Dir Direction vector for the third principal value
  \return Negative value on error, otherwise zero
*/

int FFaTensorTransforms::principalDirs (const double* S,
                                        double* pVal,
                                        double* p1Dir,
                                        double* p2Dir,
                                        double* p3Dir)
{
  double pDir[9];

  int retval = FFaTensorTransforms::principalDirs(3,S,pVal,pDir);

  p1Dir[0] = pDir[0];
  p1Dir[1] = pDir[1];
  p1Dir[2] = pDir[2];

  p2Dir[0] = pDir[3];
  p2Dir[1] = pDir[4];
  p2Dir[2] = pDir[5];

  p3Dir[0] = pDir[6];
  p3Dir[1] = pDir[7];
  p3Dir[2] = pDir[8];

  return retval;
}

/*!
  \param[in] S Tensor values, {s11, s22, s12}
  \param[out] pVal Principal values, {pMax,pMin}
  \param[out] p1Dir Direction vector for the first principal value
  \param[out] p2Dir Direction vector for the second principal value
  \return Negative value on error, otherwise zero
*/

int FFaTensorTransforms::principalDirs (const double* S,
                                        double* pVal,
                                        double* p1Dir,
                                        double* p2Dir)
{
  double pDir[4];

  int retval = FFaTensorTransforms::principalDirs(2,S,pVal,pDir);

  p1Dir[0] = pDir[0];
  p1Dir[1] = pDir[1];

  p2Dir[0] = pDir[2];
  p2Dir[1] = pDir[3];

  return retval;
}


/*!
  \param[in] N Tensor dimension (1, 2 or 3)
  \param[in] S Tensor values, upper triangle of the symmetric tensor
  \param[out] P Principal values in descending order, i.e., {P1,P2,P3}

  \details The tensor values are assumed laid out as follows
  - N = 1: S = {s11}
  - N = 2: S = {s11, s22, s12}
  - N = 3: S = {s11, s22, s33, s12, s13, s23}
*/

bool FFaTensorTransforms::principalValues (int N, const double* S, double* P)
{
  switch (N)
  {
    case 1:
      P[0] = S[0];
      return true;

    case 2:
      return FFaTensorTransforms::principalVals2D(S[0],S[1],S[2],P);
    case 3:
      return FFaTensorTransforms::principalVals3D(S[0],S[1],S[2],
                                                  S[3],S[4],S[5],P);
  }
  return false;
}

bool FFaTensorTransforms::principalVals2D (double s11, double s22, double s12,
                                           double* pVals)
{
  double C = -(s11+s22);
  double D = s11*s22 - s12*s12;

  if (FFa::cubicSolve(0.0,1.0,C,D,pVals) != 2)
  {
    std::cerr <<" *** FFaTensorTransforms::principalVals2D:"
              <<" Could not solve the quadratic equation  "
              <<" x^2 + ("<< C <<")*x + ("<< D <<") = 0"<< std::endl;
    return false;
  }

  // Sort in decreasing order
  if (pVals[0] < pVals[1]) std::swap(pVals[0],pVals[1]);

  return true;
}

bool FFaTensorTransforms::principalVals3D (double s11, double s22, double s33,
                                           double s12, double s13, double s23,
                                           double* pVals)
{
  double B = -(s11+s22+s33);
  double C = s11*s22 + s11*s33 + s22*s33 - s12*s12 - s13*s13 - s23*s23;
  double D = s11*s23*s23 + s22*s13*s13 + s33*s12*s12 - s11*s22*s33
             - 2.0*s12*s13*s23;

  if (FFa::cubicSolve(1.0,B,C,D,pVals) != 3)
  {
    std::cerr <<" *** FFaTensorTransforms::principalVals3D:"
              <<" Could not solve the cubic equation  "
              <<" x^3 + ("<< B <<")*x^2 + ("<< C << ")*x + ("<< D
              <<") =  0"<< std::endl;
    return false;
  }

  // Sort in decreasing order
  if (pVals[0] < pVals[1]) std::swap(pVals[0],pVals[1]);
  if (pVals[1] < pVals[2]) std::swap(pVals[1],pVals[2]);
  if (pVals[0] < pVals[1]) std::swap(pVals[0],pVals[1]);

  return true;
}


/*!
  \param[in] pMax Max principal value
  \param[in] pMin Min principal value
  \return Max shear value
*/

double FFaTensorTransforms::maxShearValue (double pMax, double pMin)
{
  return 0.5*(pMax-pMin);
}

/*!
  \param[in] N Vector dimension (1, 2 or 3)
  \param[in] pMaxDir Direction of the max principal value
  \param[in] pMinDir Direction of the min principal value
  \param[out] maxShearDir Direction of the max shear value,
              this is 45 degrees on the max and min directions
	      pointing from min to max
*/

void FFaTensorTransforms::maxShearDir (int N,
                                       const double* pMaxDir,
                                       const double* pMinDir,
                                       double* maxShearDir)
{
  int    i;
  double length = 0;

  for (i = 0; i < N; i++)
  {
    maxShearDir[i] = pMaxDir[i] - pMinDir[i];
    length += maxShearDir[i] * maxShearDir[i];
  }

  length = sqrt(length);
  for (i = 0; i < N; i++)
    maxShearDir[i] /= length;
}


/*!
  \param[in] inTensor Symmetric input tensor, {s11, s22, s12}
  \param[in] rotMx Transformation matrix, [eX, eY]
  \param[out] outTensor The transformed tensor
*/

void FFaTensorTransforms::rotate2D (const double* inTensor, const double* rotMx,
                                    double* outTensor)
{
  FFaTensorTransforms::rotate(inTensor, &rotMx[0], &rotMx[2], outTensor);
}

/*!
  \param[in] S Symmetric input tensor, S = {s11, s22, s12}
  \param[in] eX Unit X-direction vector
  \param[in] eY Unit Y-direction vector
  \param[out] S_transformed The transformed tensor
*/

void FFaTensorTransforms::rotate (const double* S,
                                  const double* eX,
                                  const double* eY,
                                  double* S_transformed)
{
  double TS11 = eX[0]*S[0] + eY[0]*S[2];
  double TS12 = eX[0]*S[2] + eY[0]*S[1];
  double TS21 = eX[1]*S[0] + eY[1]*S[2];
  double TS22 = eX[1]*S[2] + eY[1]*S[1];

  S_transformed[0] = TS11*eX[0] + TS12*eY[0];
  S_transformed[1] = TS21*eX[1] + TS22*eY[1];
  S_transformed[2] = TS11*eX[1] + TS12*eY[1];
}

/*!
  \param[in] inTensor Symmetric input tensor, {s11, s22, s33, s12, s13, s23}
  \param[in] rotMx Transformation matrix, [eX, eY, eZ]
  \param[out] outTensor The transformed tensor
*/

void FFaTensorTransforms::rotate3D (const double* inTensor, const double* rotMx,
                                    double* outTensor)
{
  FFaTensorTransforms::rotate (inTensor,&rotMx[0],&rotMx[3],&rotMx[6],outTensor);
}

/*!
  \param[in] S Symmetric input tensor, S = {s11, s22, s33, s12, s13, s23}
  \param[in] eX Unit X-direction vector
  \param[in] eY Unit Y-direction vector
  \param[in] eZ Unit Z-direction vector
  \param[out] S_transformed The transformed tensor
*/

void FFaTensorTransforms::rotate (const double* S,
                                  const double* eX,
                                  const double* eY,
                                  const double* eZ,
                                  double* S_transformed)
{
  double TS11 = eX[0]*S[0] + eY[0]*S[3] + eZ[0]*S[4];
  double TS12 = eX[0]*S[3] + eY[0]*S[1] + eZ[0]*S[5];
  double TS13 = eX[0]*S[4] + eY[0]*S[5] + eZ[0]*S[2];
  double TS21 = eX[1]*S[0] + eY[1]*S[3] + eZ[1]*S[4];
  double TS22 = eX[1]*S[3] + eY[1]*S[1] + eZ[1]*S[5];
  double TS23 = eX[1]*S[4] + eY[1]*S[5] + eZ[1]*S[2];
  double TS31 = eX[2]*S[0] + eY[2]*S[3] + eZ[2]*S[4];
  double TS32 = eX[2]*S[3] + eY[2]*S[1] + eZ[2]*S[5];
  double TS33 = eX[2]*S[4] + eY[2]*S[5] + eZ[2]*S[2];

  S_transformed[0] = TS11*eX[0] + TS12*eY[0] + TS13*eZ[0];
  S_transformed[1] = TS21*eX[1] + TS22*eY[1] + TS23*eZ[1];
  S_transformed[2] = TS31*eX[2] + TS32*eY[2] + TS33*eZ[2];
  S_transformed[3] = TS11*eX[1] + TS12*eY[1] + TS13*eZ[1];
  S_transformed[4] = TS11*eX[2] + TS12*eY[2] + TS13*eZ[2];
  S_transformed[5] = TS21*eX[2] + TS22*eY[2] + TS23*eZ[2];
}


/*!
  Padding with zeros to fit the 3D tensor.
*/

void FFaTensorTransforms::from2Dto3D (const double* S2D, double* S3D)
{
  S3D[0] = S2D[0];
  S3D[1] = S2D[1];
  S3D[2] = 0.0;
  S3D[3] = S2D[2];
  S3D[4] = 0.0;
  S3D[5] = 0.0;
}

/*!
  Just cutting off whats needed to fit the 2D tensor.
*/

void FFaTensorTransforms::from3Dto2D (const double* S3D, double* S2D)
{
  S2D[0] = S3D[0];
  S2D[1] = S3D[1];
  S2D[2] = S3D[3];
}
