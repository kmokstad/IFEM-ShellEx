// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFaMath.H
  \brief Various math utility functions (extension of math.h).
*/

#ifndef FFA_MATH_H
#define FFA_MATH_H

#include <cmath>

#if _MSC_VER > 1300 && _MSC_VER < 1500
#define hypot _hypot
#endif

#ifndef EPS_ZERO
//! \brief Division by zero tolerance.
#define EPS_ZERO 1.0e-15
#endif

#ifndef M_PI
//! \brief The value of the numerical constant &pi;.
#define M_PI 3.14159265358979323846
#endif

#ifndef RAD
//! \brief Degrees-to-Radians conversion.
#define RAD(x) (x)*M_PI/180.0
#endif

//! \brief Rounds off \a value down to \a precision significant digits.
extern double round (double value, int precision);

#ifdef FFA_DEBUG
extern double atan3 (double y, double x, const char* func = 0);
#else
//! \brief Checks for singular rotations before using atan2 library function.
inline double atan3 (double y, double x, const char* = 0)
{
  return fabs(y) > EPS_ZERO || fabs(x) > EPS_ZERO ? atan2(y,x) : 0.0;
}
#endif


namespace FFa //! Various math utilities
{
  //! \brief Solves the cubic equation
  //! A*x<sup>3</sup> + B*x<sup>2</sup> + C*x + D = 0.
  int cubicSolve (double A, double B, double C, double D, double* X);
  //! \brief Solves a bi-linear set of equations in \a X and \a Y.
  int bilinearSolve (const double* A, const double* B, double* X, double* Y);
}

#endif
