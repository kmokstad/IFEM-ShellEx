// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFaMath.C
  \brief Various math utility functions.
*/

#include "FFaLib/FFaAlgebra/FFaMath.H"
#include <iostream>


#ifdef FFA_DEBUG
double atan3 (double y, double x, const char* func)
{
  if (fabs(y) > EPS_ZERO || fabs(x) > EPS_ZERO)
    return atan2(y,x);

  if (func)
    std::cerr << func <<": Singular rotation (x,y)="<< x <<","<< y << std::endl;

  return 0.0;
}
#endif


#ifndef FFA_NO_ROUND
double round (double value, int precision)
{
  if (value == 0.0) return 0.0; // Avoid log10(0)
  double aval  = std::fabs(value);
  int exponent = (int)std::log10(aval) - (aval < 1.0 ? precision : precision-1);
  double denom = std::pow(10.0,exponent);
#ifdef FFA_USE_ROUND
  return std::round(value/denom)*denom;
#else
  double delta = value > 0.0 ? 0.5*denom : -0.5*denom;
  double remnd = fmod(value+delta,denom);
  // fmod seems to return <denom> (or close to) if <value> is already rounded
  return fabs(remnd) < 0.99*denom ? value+delta - remnd : value+delta;
#endif
}
#endif


/*!
  See K. Rottmann, "Matematische Formelsammlung" (1960), pp. 13-16 for details.

  \param[in] A Coefficient of the cubic term
  \param[in] B Coefficient of the quadratic term
  \param[in] C Coefficient of the linear term
  \param[in] D The constant term
  \param[out] X Solution(s)
  \return Number of solutions found
  - -3 : No real solution, three complex roots detected
  - -2 : No real solution, two complex roots detected
  -  0 : No solutions at all (A, B and C were all zero)
  -  1 : One solution found (A and B were zero)
  -  2 : Two unique solution found (A was zero)
  -  3 : Three unique solutions found
*/

int FFa::cubicSolve (double A, double B, double C, double D, double* X)
{
  const double epsilon = 1.0e-16;

  if (fabs(A) > epsilon) // Cubic equation
  {
    const double epsmall = pow(epsilon,6.0);
    double P =  (C - B*B/(3.0*A)) / (3.0*A);
    double Q = ((2.0*B*B/(27.0*A) - C/3.0) * B/A + D) / (A+A);
    double W = Q*Q + P*P*P;

    if (W <= -epsmall && P < 0.0) // Casus irreducibilis, case delta
    {
      double FI = acos(-Q/sqrt(-P*P*P));
      X[0] =  2.0*sqrt(-P) * cos( FI      /3.0);
      X[1] = -2.0*sqrt(-P) * cos((FI+M_PI)/3.0);
      X[2] = -2.0*sqrt(-P) * cos((FI-M_PI)/3.0);
    }
    else if (fabs(W) < epsmall && Q <= 0.0) // Kardanische formel, case alpha
    {
      X[0] =  2.0 * pow(-Q,1.0/3.0);
      X[1] = -0.5 * X[0];
      X[2] =  X[1];
    }
    else if (W > -epsmall && Q+sqrt(W) <= 0.0 && Q-sqrt(W) <= 0.0) // Case alpha
    {
      X[0] = pow(-Q+sqrt(W),1.0/3.0) + pow(-Q-sqrt(W),1.0/3.0);
      X[1] = -0.5 * X[0]; // Ignore imaginary part of the complex solution
      X[2] =  X[1];
    }
    else if (W >= epsmall && fabs(Q) > epsmall && P > 0.0) // Case beta
    {
      // TODO: Try to simplify this.
      // TODO: The fabs(Q) and copysign-part was added 17.03.2020, verify it
      double FI = atan(sqrt(P*P*P)/fabs(Q));
      double KI = atan(copysign(pow(tan(0.5*FI),1.0/3.0),Q));
      X[0] = -2.0 * sqrt(P) / tan(KI+KI);
      X[1] = -0.5 * X[0]; // Ignore imaginary part of the complex solution
      X[2] =  X[1];
    }
    else if (W > -epsmall && fabs(Q) > epsmall && P < 0.0) // Case gamma
    {
      std::cerr <<"FFa::cubicSolve: Case gamma not implemented."<< std::endl;
      return -3;
    }
    else
      return -3;

    W = B/(3.0*A);
    for (int i = 0; i < 3; i++)
      X[i] -= W;

    return 3;
  }
  else if (fabs(B) > epsilon) // Quadratic equation
  {
    const double epsmall = pow(epsilon,4.0);
    double P = C*C - 4.0*B*D;
    if (P > 0.0)
    {
      double Q = sqrt(P);
      X[0] = (-C+Q)/(B+B);
      X[1] = (-C-Q)/(B+B);
    }
    else if (P > -epsmall)
    {
      X[0] = -C/(B+B);
      X[1] = X[0];
    }
    else
      return -2;

    return 2;
  }
  else if (fabs(C) > epsilon) // Linear equation
  {
    X[0] = -D/C;
    return 1;
  }

  return 0;
}


/*!
  The following set of equations are solved for the unknowns \a X and \a Y :
  \code
           A0 * x*y  +  A1 * x  +  A2 * y  =  A3
           B0 * x*y  +  B1 * x  +  B2 * y  =  B3
  \endcode
  \param[in] A Coefficients of the first equation
  \param[in] B Coefficients of the second equation
  \param[out] X Solution value(s) of the first unknown
  \param[out] Y Solution value(s) of the second unknown
  \return Number of unique solutions found, negative on error
*/

int FFa::bilinearSolve (const double* A, const double* B, double* X, double* Y)
{
  const double epsilon = 1.0e-16;

  // Initialize the geometric tolerance
  double tol = 0.0;
  for (int i = 0; i < 4; i++)
  {
    double a = epsilon*fabs(A[i]);
    double b = epsilon*fabs(B[i]);
    if (a > tol) tol = a;
    if (b > tol) tol = b;
  }

  int nsol = 0;
  double D = A[1]*B[2] - B[1]*A[2];
  if (fabs(A[0]) < tol && fabs(B[0]) < tol)
  {
    // Linear set of equations:
    // | A1 A2 | * | X | = | A3 |
    // | B1 B2 |   | Y | = | B3 |
    if (fabs(D) > tol*tol)
    {
      X[0] = ( B[2]*A[3] - A[2]*B[3]) / D;
      Y[0] = (-B[1]*A[3] + A[1]*B[3]) / D;
      nsol = 1;
    }
  }
  else
  {
    // Solve the second-order equation  Q2*x^2 + Q1*x + Q0 = 0

    double Q0 = B[2]*A[3] - A[2]*B[3];
    double Q1 = B[0]*A[3] - A[0]*B[3] - D;
    double Q2 = A[0]*B[1] - B[0]*A[1];
    double Z[3] = { 0.0, 0.0, 0.0 };
    int nx = cubicSolve(0.0,Q2,Q1,Q0,Z);

    // Find y = (A3-A1*x) / (A0*x+A2)
    for (int i = 0; i < nx; i++)
    {
      Q0 = A[0]*Z[i] + A[2];
      if (fabs(Q0) > tol)
      {
        X[nsol] = Z[i];
        Y[nsol] = (A[3] - A[1]*Z[i]) / Q0;
        ++nsol;
      }
    }

    // Solve the second-order equation  Q2*y^2 + Q1*y + Q0 = 0

    Q0 = B[1]*A[3] - A[1]*B[3];
    Q1 = B[0]*A[3] - A[0]*B[3] + D;
    Q2 = A[0]*B[2] - B[0]*A[2];
    int ny = cubicSolve(0.0,Q2,Q1,Q0,Z);

    // Find x = (A3-A2*y) / (A0*y+A1)
    int n1 = nsol;
    for (int j = 0; j < ny; j++)
    {
      Q0 = A[0]*Z[j] + A[1];
      bool newSol = fabs(Q0) > tol;
      for (int i = 0; i < n1 && newSol; i++)
      {
        // Avoid solutions equal to those found by solving for x first
        double aY = fabs(Y[i]);
        double aZ = fabs(Z[j]);
        newSol = fabs(Y[i]-Z[j]) > epsilon*(aY > aZ ? aY : aZ);
      }
      if (newSol)
      {
        Y[nsol] = Z[j];
        X[nsol] = (A[3] - A[2]*Z[j]) / Q0;
        ++nsol;
      }
    }

    if (nx < 0)
      nsol = nx;
    else if (ny < 0)
      nsol = ny;
  }

  return nsol;
}
