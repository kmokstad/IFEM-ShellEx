// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFaVec3.C
  \brief Point vectors in 3D space.
*/

#include <cctype>

#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include "FFaLib/FFaAlgebra/FFaMath.H"


////////////////////////////////////////////////////////////////////////////////
//
// Local operators
//
////////////////////////////////////////////////////////////////////////////////

FaVec3& FaVec3::operator+= (const FaVec3& v)
{
  n[VX] += v.n[VX];
  n[VY] += v.n[VY];
  n[VZ] += v.n[VZ];
  return *this;
}


FaVec3& FaVec3::operator-= (const FaVec3& v)
{
  n[VX] -= v.n[VX];
  n[VY] -= v.n[VY];
  n[VZ] -= v.n[VZ];
  return *this;
}


FaVec3& FaVec3::operator*= (double d)
{
  n[VX] *= d;
  n[VY] *= d;
  n[VZ] *= d;
  return *this;
}


FaVec3& FaVec3::operator/= (double d)
{
  if (fabs(d) < EPS_ZERO)
  {
#ifdef FFA_DEBUG
    std::cerr <<"FaVec3::operator/=(double): Division by zero, d="<< d
              << std::endl;
#endif
    n[VX] = n[VY] = n[VZ] = HUGE_VAL;
  }
  else for (short int i = 0; i < 3; i++)
    n[i] /= d;

  return *this;
}


////////////////////////////////////////////////////////////////////////////////
//
// Special functions
//
////////////////////////////////////////////////////////////////////////////////

/*!
  Returns 0 if not parallel, 1 if same direction and -1 if opposite direction.
  Tolerance = 1 - cos(maxAngle)
*/

int FaVec3::isParallell (const FaVec3& otherVec, double tolerance) const
{
  double len = sqrt(this->sqrLength() * otherVec.sqrLength());
  if (fabs(len) < EPS_ZERO) return 0;

  double res = (*this * otherVec) / len;
  if (fabs(res+1.0) <= tolerance) return -1;
  if (fabs(res-1.0) <= tolerance) return 1;

  return 0;
}


/*!
  The angle is returned in radians in the range [0,&pi;].
*/

double FaVec3::angle (const FaVec3& otherVec) const
{
  double len = sqrt(this->sqrLength() * otherVec.sqrLength());
  if (fabs(len) < EPS_ZERO)
  {
#ifdef FFA_DEBUG
    std::cerr <<"FaVec3::angle(FaVec3&): Zero vector length(s), |*this|="
              << this->length() <<" |v|="<< otherVec.length() << std::endl;
#endif
    return 0.0;
  }

  double res = (*this * otherVec) / len;
  if (res >  1.0) return 0.0;
  if (res < -1.0) return acos(-1.0);

  return acos(res);
}


/*!
  Checks if the given vector is equal to this one within the given tolerance.
*/

bool FaVec3::equals (const FaVec3& otherVec, double tolerance) const
{
  double xd = n[VX] - otherVec.n[VX];
  if (fabs(xd) > tolerance) return false;

  double yd = n[VY] - otherVec.n[VY];
  if (fabs(yd) > tolerance) return false;

   double zd = n[VZ] - otherVec.n[VZ];
  if (fabs(zd) > tolerance) return false;

  return xd*xd+yd*yd+zd*zd > tolerance*tolerance ? false : true;
}


/*!
  Checks if this vector is zero within the given tolerance.
*/

bool FaVec3::isZero (double tolerance) const
{
  if (fabs(n[VX]) > tolerance) return false;
  if (fabs(n[VY]) > tolerance) return false;
  if (fabs(n[VZ]) > tolerance) return false;
  return true;
}


/*!
  Truncates components of this vector to zero if less than the tolerance.
*/

FaVec3& FaVec3::truncate (double tolerance)
{
  double len = this->length();
  if (len > EPS_ZERO && tolerance > 0.0)
  {
    if (len < 1.0) tolerance *= len;
    for (short int i = 0; i < 3; i++)
      if (fabs(n[i]) < tolerance)
        n[i] = 0.0;
  }
  return *this;
}


/*!
  If the vector has zero length it will end up as [1,0,0].
  If a vector component is less than \a truncTol in absolute value
  after the normalization, it is set to zero.
*/

FaVec3& FaVec3::normalize (double truncTol)
{
  double len = this->length();
  if (len < EPS_ZERO)
  {
    // Default direction if zero length
    n[VX] = 1.0;
    n[VY] = 0.0;
    n[VZ] = 0.0;
#ifdef FFA_DEBUG
    std::cerr <<"FaVec3::normalize(): Zero vector length, len="<< len
              <<"\n                     Replacing by default [1,0,0]"
              << std::endl;
#endif
  }
  else for (short int i = 0; i < 3; i++)
  {
    n[i] /= len;
    if (truncTol > 0.0 && fabs(n[i]) < truncTol)
      n[i] = 0.0;
  }
  return *this;
}


/*!
  This method is used when parsing Nastran bulk data files, to avoid
  (or at least reduce the probability of) checksum issues.
*/

FaVec3& FaVec3::round (int precision)
{
  for (short int i = 0; i < 3; i++)
    ::round(n[i],precision);

  return *this;
}


////////////////////////////////////////////////////////////////////////////////
//
// Cylindrical and spherical coordinate system operations
//
////////////////////////////////////////////////////////////////////////////////

/*!
  Assuming \a cylCoords = [radius, angleAboutAxis, lengthAlongRotAxis];
  The \a axis argument defines which axis of rotation to use.
*/

FaVec3& FaVec3::setByCylCoords (const FaVec3& cylCoords, FFaVec3IdxEnum axis)
{
  if (axis < 0 || axis > 2) axis = VZ;

  // Cyclic permutation:  indX -> indY -> axis
  //                         \____________/
  short int indX = (axis+1) % 3;
  short int indY = (axis+2) % 3;

  n[indX] = cylCoords[VX] * cos(cylCoords[VY]);
  n[indY] = cylCoords[VX] * sin(cylCoords[VY]);
  n[axis] = cylCoords[VZ];

  return *this;
}


/*!
  Returns v = [radius, angleAboutAxis, lengthAlongRotAxis].
  The \a axis parameter defines which axis of rotation to use.
*/

FaVec3 FaVec3::getAsCylCoords (FFaVec3IdxEnum axis) const
{
  if (axis < 0 || axis > 2) axis = VZ;

  // Cyclic permutation:  indX -> indY -> axis
  //                         \____________/
  short int indX = (axis+1) % 3;
  short int indY = (axis+2) % 3;

  double r     = hypot(n[indX],n[indY]);
  double theta = atan3(n[indY],n[indX]);

  return FaVec3(r,theta,n[axis]);
}


/*!
  Assuming \a sphCoords = [radius, angleAboutAxis, asimuthAngle];
  The \a axis parameter defines which axis of rotation to use.
*/

FaVec3& FaVec3::setBySphCoords (const FaVec3& sphCoords, FFaVec3IdxEnum axis)
{
  if (axis < 0 || axis > 2) axis = VZ;

  // Cyclic permutation:  indX -> indY -> axis
  //                         \____________/
  short int indX = (axis+1) % 3;
  short int indY = (axis+2) % 3;

  n[indX] = sphCoords[VX] * cos(sphCoords[VY]) * sin(sphCoords[VZ]);
  n[indY] = sphCoords[VX] * sin(sphCoords[VY]) * sin(sphCoords[VZ]);
  n[axis] = sphCoords[VX] * cos(sphCoords[VZ]);

  return *this;
}


/*!
  Returns v = [radius, angleAboutAxis, asimuthAngle].
  The \a axis parameter defines which axis of rotation to use.
*/

FaVec3 FaVec3::getAsSphCoords (FFaVec3IdxEnum axis) const
{
  double r = this->length();
  if (r < EPS_ZERO) return FaVec3();

  if (axis < 0 || axis > 2) axis = VZ;

  // Cyclic permutation:  indX -> indY -> axis
  //                         \____________/
  short int indX = (axis+1) % 3;
  short int indY = (axis+2) % 3;

  double theta = atan3(n[indY],n[indX]);
  double phi   = acos(n[axis]/r);

  return FaVec3(r,theta,phi);
}


////////////////////////////////////////////////////////////////////////////////
//
// Global operators
//
////////////////////////////////////////////////////////////////////////////////

FaVec3 operator- (const FaVec3& a)
{
  return FaVec3(-a.n[VX],-a.n[VY],-a.n[VZ]);
}

FaVec3 operator+ (const FaVec3& a, const FaVec3& b)
{
  return FaVec3(a.n[VX]+b.n[VX], a.n[VY]+b.n[VY], a.n[VZ]+b.n[VZ]);
}

FaVec3 operator- (const FaVec3& a, const FaVec3& b)
{
  return FaVec3(a.n[VX]-b.n[VX], a.n[VY]-b.n[VY], a.n[VZ]-b.n[VZ]);
}


FaVec3 operator* (const FaVec3& a, double d)
{
  return FaVec3(d*a.n[VX], d*a.n[VY], d*a.n[VZ]);
}

FaVec3 operator* (double d, const FaVec3& a)
{
  return FaVec3(d*a.n[VX], d*a.n[VY], d*a.n[VZ]);
}

FaVec3 operator/ (const FaVec3& a, double d)
{
  if (fabs(d) < EPS_ZERO)
  {
#if FFA_DEBUG
    std::cerr <<"FaVec3 operator/(FaVec3&,double): Division by zero, d="<< d
              << std::endl;
#endif
    return FaVec3(HUGE_VAL, HUGE_VAL, HUGE_VAL);
  }

  return FaVec3(a.n[VX]/d, a.n[VY]/d, a.n[VZ]/d);
}


double operator* (const FaVec3& a, const FaVec3& b)
{
  return a.n[VX]*b.n[VX] + a.n[VY]*b.n[VY] + a.n[VZ]*b.n[VZ];
}

FaVec3 operator^ (const FaVec3& a, const FaVec3& b)
{
  return FaVec3(a.n[VY]*b.n[VZ] - a.n[VZ]*b.n[VY],
                a.n[VZ]*b.n[VX] - a.n[VX]*b.n[VZ],
                a.n[VX]*b.n[VY] - a.n[VY]*b.n[VX]);
}


bool operator== (const FaVec3& a, const FaVec3& b)
{
  return (a.n[VX] == b.n[VX]) &&
         (a.n[VY] == b.n[VY]) &&
         (a.n[VZ] == b.n[VZ]);
}

bool operator!= (const FaVec3& a, const FaVec3& b)
{
  return (a.n[VX] != b.n[VX]) ||
         (a.n[VY] != b.n[VY]) ||
         (a.n[VZ] != b.n[VZ]);
}


std::ostream& operator<< (std::ostream& s, const FaVec3& v)
{
  s << v.n[VX] <<' '<< v.n[VY] <<' '<< v.n[VZ];
  return s;
}

std::istream& operator>> (std::istream& s, FaVec3& v)
{
  FaVec3 v_tmp;
  s >> v_tmp[VX] >> v_tmp[VY] >> v_tmp[VZ];
  if (s) v = v_tmp;
  return s;
}
