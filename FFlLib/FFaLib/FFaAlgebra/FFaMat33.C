// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFaLib/FFaAlgebra/FFaMat33.H"
#include "FFaLib/FFaAlgebra/FFaMath.H"
#include "FFaLib/FFaOS/FFaFortran.H"


/*!
  Matrix layout :

    [0][0]  [1][0]  [2][0]

    [0][1]  [1][1]  [2][1]

    [0][2]  [1][2]  [2][2]
*/


////////////////////////////////////////////////////////////////////////////////
//
// Constructors
//
////////////////////////////////////////////////////////////////////////////////

FaMat33::FaMat33 (const float* mat)
{
  v[0] = FaVec3(mat);
  v[1] = FaVec3(mat+3);
  v[2] = FaVec3(mat+6);
}


FaMat33::FaMat33 (const double* mat)
{
  v[0] = FaVec3(mat);
  v[1] = FaVec3(mat+3);
  v[2] = FaVec3(mat+6);
}


FaMat33::FaMat33 (const FaVec3& v0, const FaVec3& v1, const FaVec3& v2)
{
  v[0] = v0;
  v[1] = v1;
  v[2] = v2;
}


////////////////////////////////////////////////////////////////////////////////
//
// Local operators
//
////////////////////////////////////////////////////////////////////////////////

FaMat33& FaMat33::operator= (const FaMat33& m)
{
  v[0] = m.v[0];
  v[1] = m.v[1];
  v[2] = m.v[2];
  return *this;
}


FaMat33& FaMat33::operator+= (const FaMat33& m)
{
  v[0] += m.v[0];
  v[1] += m.v[1];
  v[2] += m.v[2];
  return *this;
}


FaMat33& FaMat33::operator-= (const FaMat33& m)
{
  v[0] -= m.v[0];
  v[1] -= m.v[1];
  v[2] -= m.v[2];
  return *this;
}


FaMat33& FaMat33::operator*= (double d)
{
  v[0] *= d;
  v[1] *= d;
  v[2] *= d;
  return *this;
}


FaMat33& FaMat33::operator/= (double d)
{
  if (fabs(d) < EPS_ZERO)
  {
#ifdef FFA_DEBUG
    std::cerr <<"FaMat33::operator/=(double): Division by zero "<< d
	      << std::endl;
#endif
    v[0] = v[1] = v[2] = FaVec3(HUGE_VAL, HUGE_VAL, HUGE_VAL);
  }
  else
  {
    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
  }

  return *this;
}


////////////////////////////////////////////////////////////////////////////////
//
// Special functions
//
////////////////////////////////////////////////////////////////////////////////

#define THIS(i,j) this->operator()(i,j)

FaMat33 FaMat33::inverse (double eps) const
{
  static FaMat33 b;

  double det
    = v[0][0]*(v[1][1]*v[2][2] - v[2][1]*v[1][2])
    - v[0][1]*(v[1][0]*v[2][2] - v[2][0]*v[1][2])
    + v[0][2]*(v[1][0]*v[2][1] - v[2][0]*v[1][1]);

  if (fabs(det) >= eps)
  {
    b.v[0][0] =  (v[1][1]*v[2][2] - v[2][1]*v[1][2]) / det;
    b.v[0][1] = -(v[0][1]*v[2][2] - v[2][1]*v[0][2]) / det;
    b.v[0][2] =  (v[0][1]*v[1][2] - v[1][1]*v[0][2]) / det;
    b.v[1][0] = -(v[1][0]*v[2][2] - v[2][0]*v[1][2]) / det;
    b.v[1][1] =  (v[0][0]*v[2][2] - v[2][0]*v[0][2]) / det;
    b.v[1][2] = -(v[0][0]*v[1][2] - v[1][0]*v[0][2]) / det;
    b.v[2][0] =  (v[1][0]*v[2][1] - v[2][0]*v[1][1]) / det;
    b.v[2][1] = -(v[0][0]*v[2][1] - v[2][0]*v[0][1]) / det;
    b.v[2][2] =  (v[0][0]*v[1][1] - v[1][0]*v[0][1]) / det;
  }
#ifdef FFA_DEBUG
  else
    std::cerr <<"FaMat33::inverse(): Singular matrix, det = "<< std::endl;
#endif

  return b;
}


FaMat33& FaMat33::setIdentity ()
{
  v[0] = FaVec3(1.0, 0.0, 0.0);
  v[1] = FaVec3(0.0, 1.0, 0.0);
  v[2] = FaVec3(0.0, 0.0, 1.0);
  return *this;
}


FaMat33 FaMat33::transpose () const
{
  return FaMat33(FaVec3(v[0][0], v[1][0], v[2][0]),
                 FaVec3(v[0][1], v[1][1], v[2][1]),
                 FaVec3(v[0][2], v[1][2], v[2][2]));
}


FaMat33& FaMat33::shift (int delta)
{
  if (delta < -2 || delta%3 == 0) return *this;

  // Perform a cyclic permutation of the matrix columns
  FaVec3 v1 = v[0];
  FaVec3 v2 = v[1];
  FaVec3 v3 = v[2];

  v[(3+delta)%3] = v1;
  v[(4+delta)%3] = v2;
  v[(5+delta)%3] = v3;

  return *this;
}


bool FaMat33::isCoincident (const FaMat33& oMat, double tolerance) const
{
  if (v[0].isParallell(oMat[0],tolerance) != 1) return false;
  if (v[1].isParallell(oMat[1],tolerance) != 1) return false;
  if (v[2].isParallell(oMat[2],tolerance) != 1) return false;
  return true;
}


/*!
  Compute a globalized coordinate system where the X-axis is parallel to the
  given vector v1, and the two other axes are as close as possible to the
  corresponding global coordinate axes.
*/

FaMat33& FaMat33::makeGlobalizedCS (const FaVec3& v1)
{
  FaVec3& eX = v[0];
  FaVec3& eY = v[1];
  FaVec3& eZ = v[2];

  eX = v1;
  eX.normalize();

  if (fabs(eX(3)) > fabs(eX(2)))
  {
    // Define eY by projecting the global Y-axis onto the plane defined by eX
    // eY = v1 x (Y x v1)
    eY(1) = -eX(2)*eX(1);
    eY(2) =  eX(1)*eX(1) + eX(3)*eX(3);
    eY(3) = -eX(2)*eX(3);

    eY.normalize();
    eZ = eX ^ eY;
  }
  else
  {
    // Define eZ by projecting the global Z-axis onto the plane defined by eX
    // eZ = v1 x (Z x v1)
    eZ(1) = -eX(3)*eX(1);
    eZ(2) = -eX(3)*eX(2);
    eZ(3) =  eX(1)*eX(1) + eX(2)*eX(2);

    eZ.normalize();
    eY = eZ ^ eX;
  }

  return *this;
}


/*!
  Compute a globalized coordinate system in the plane defined by the two given
  vectors v1 and v2, such that the local Z-axis is parallell to the normal
  vector of the plane, and the two other axes are as close as possible to the
  corresponding global coordinate axes.
*/

FaMat33& FaMat33::makeGlobalizedCS (const FaVec3& v1, const FaVec3& v2)
{
  FaVec3& eX = v[0];
  FaVec3& eY = v[1];
  FaVec3& eZ = v[2];

  eZ = v1 ^ v2;
  eZ.normalize();

  if (fabs(eZ(1)) < fabs(eZ(2)))
  {
    // Define eX by projecting the global X-axis onto the plane defined by eZ
    // eX = eZ x (X x eZ)
    eX(1) =  eZ(2)*eZ(2) + eZ(3)*eZ(3);
    eX(2) = -eZ(1)*eZ(2);
    eX(3) = -eZ(1)*eZ(3);

    eX.normalize();
    eY = eZ ^ eX;
  }
  else
  {
    // Define eY by projecting the global Y-axis onto the plane defined by eZ
    // eY = eZ x (Y x eZ)
    eY(1) = -eZ(2)*eZ(1);
    eY(2) =  eZ(1)*eZ(1) + eZ(3)*eZ(3);
    eY(3) = -eZ(2)*eZ(3);

    eY.normalize();
    eX = eY ^ eZ;
  }

  return *this;
}

FaMat33& FaMat33::makeGlobalizedCS (const FaVec3& v0,
                                    const FaVec3& v1, const FaVec3& v2)
{
  return makeGlobalizedCS(v1-v0,v2-v0);
}

FaMat33& FaMat33::makeGlobalizedCS (const FaVec3& v1, const FaVec3& v2,
                                    const FaVec3& v3, const FaVec3& v4)
{
  return makeGlobalizedCS(v3-v1,v4-v2);
}


/*!
  Compute an incremental rotation tensor from the given Euler angles.
*/

FaMat33& FaMat33::eulerRotateZYX (const FaVec3& angles)
{
  double ca = cos(angles(3));
  double cb = cos(angles(2));
  double cy = cos(angles(1));

  double sa = sin(angles(3));
  double sb = sin(angles(2));
  double sy = sin(angles(1));

  THIS(1,1) = ca*cb;
  THIS(1,2) = ca*sb*sy - sa*cy;
  THIS(1,3) = ca*sb*cy + sa*sy;
  THIS(2,1) = sa*cb;
  THIS(2,2) = sa*sb*sy + ca*cy;
  THIS(2,3) = sa*sb*cy - ca*sy;
  THIS(3,1) = -sb;
  THIS(3,2) = cb*sy;
  THIS(3,3) = cb*cy;

  return *this;
}


/*!
  Return the Euler angles corresponding to an incremental rotation.
*/

FaVec3 FaMat33::getEulerZYX () const
{
  // Principles from Computer-Aided modelling, dynamic simulation and
  // dimensioning of Mechanisms, Ole Ivar Sivertsen, NTH 1995, pp. 49-51.
  /*
  double  aZ =  atan3(THIS(2,1),THIS(1,1));
  FaMat33 Tb = *this * makeZrotation(aZ);
  double  aY = -atan3(  Tb(3,1),  Tb(1,1));
  FaMat33 Tc =   Tb  * makeYrotation(aY);
  double  aX =  atan3(  Tc(3,2),  Tc(2,2));
  */
#ifdef FFA_DEBUG
  const char* func = "FaMat33::getEulerZYX";
#else
  const char* func = 0;
#endif
  // New calculation, based on http://www.udel.edu/HNES/HESC427/EULER.DOC
  // See also http://en.wikipedia.org/wiki/Euler_angles
  double aZ =  atan3(THIS(2,1),THIS(1,1),func);
  double aY = -atan3(THIS(3,1),hypot(THIS(1,1),THIS(2,1)),func);
  double aX =  atan3(THIS(3,2),THIS(3,3),func);
  return FaVec3(aX,aY,aZ);
}


/*!
  Compute an incremental rotation tensor from the given rotation angles via a
  quaternion representation of the rotation. This function is equivalent to
  subroutine vec_to_mat in the Fortran module rotationModule (vpmUtilitiesF90).
  The angles provided are those related to a Rodrigues parameterization.
  Rotation axis, with length equal to the angle to rotate about that axis.
*/

FaMat33& FaMat33::incRotate (const FaVec3& angles)
{
  const double eps = EPS_ZERO;

  double theta = angles.length();
  double quat0 = cos(0.5*theta);
  FaVec3 quatr = angles * (theta < eps ? 0.5 : sin(0.5*theta)/theta);
  double quatl = sqrt(quat0*quat0 + quatr.sqrLength());
  quat0 /= quatl;
  quatr /= quatl;

  THIS(1,1) = 2.0*(quatr(1)*quatr(1) + quat0*quat0) - 1.0;
  THIS(2,2) = 2.0*(quatr(2)*quatr(2) + quat0*quat0) - 1.0;
  THIS(3,3) = 2.0*(quatr(3)*quatr(3) + quat0*quat0) - 1.0;

  THIS(1,2) = 2.0*(quatr(1)*quatr(2) - quatr(3)*quat0);
  THIS(1,3) = 2.0*(quatr(1)*quatr(3) + quatr(2)*quat0);
  THIS(2,3) = 2.0*(quatr(2)*quatr(3) - quatr(1)*quat0);

  THIS(2,1) = 2.0*(quatr(2)*quatr(1) + quatr(3)*quat0);
  THIS(3,1) = 2.0*(quatr(3)*quatr(1) - quatr(2)*quat0);
  THIS(3,2) = 2.0*(quatr(3)*quatr(2) + quatr(1)*quat0);

  return *this;
}


/*!
  Return the rotation angles corresponding to an incremental rotation.
  This function is equivalent to the Fortran subroutine mat_to_vec
  in module rotationModule (vpmUtilitiesF90).
*/

FaVec3 FaMat33::getRotation () const
{
  const double eps = EPS_ZERO;

  int i = 1;
  if (THIS(2,2) > THIS(i,i)) i = 2;
  if (THIS(3,3) > THIS(i,i)) i = 3;

  double quat0;
  FaVec3 quatr;
  double trace = THIS(1,1) + THIS(2,2) + THIS(3,3);
  if (trace > THIS(i,i))
  {
    quat0    = 0.5*sqrt(1.0+trace);
    quatr(1) = (THIS(3,2) - THIS(2,3)) * 0.25/quat0;
    quatr(2) = (THIS(1,3) - THIS(3,1)) * 0.25/quat0;
    quatr(3) = (THIS(2,1) - THIS(1,2)) * 0.25/quat0;
  }
  else
  {
    int j = 1 + i%3;
    int k = 1 + j%3;
    quatr(i) = sqrt(0.5*THIS(i,i) + 0.25*(1.0-trace));
    quat0    = (THIS(k,j) - THIS(j,k)) * 0.25/quatr(i);
    quatr(j) = (THIS(j,i) + THIS(i,j)) * 0.25/quatr(i);
    quatr(k) = (THIS(k,i) + THIS(i,k)) * 0.25/quatr(i);
  }

  double quatl = sqrt(quat0*quat0 + quatr.sqrLength());
  double sthh  = quatr.length() / quatl;
  double cthh  = quat0 / quatl;
  double theta = sthh < 0.7 ? 2.0*asin(sthh) : 2.0*acos(cthh);
  if (theta < eps)
    return quatr * 2.0;
  else if (sthh < 1.0)
    return quatr * theta/sthh;
  else
    return quatr * theta;
}


FaMat33 FaMat33::makeZrotation (double rot)
{
  static FaMat33 r;

  double c = cos(rot);
  double s = sin(rot);

  r(1,1) = c;
  r(2,1) = s;
  r(1,2) = -s;
  r(2,2) = c;

  return r;
}

FaMat33 FaMat33::makeYrotation (double rot)
{
  static FaMat33 r;

  double c = cos(rot);
  double s = sin(rot);

  r(1,1) = c;
  r(3,1) = -s;
  r(1,3) = s;
  r(3,3) = c;

  return r;
}

FaMat33 FaMat33::makeXrotation (double rot)
{
  static FaMat33 r;

  double c = cos(rot);
  double s = sin(rot);

  r(2,2) = c;
  r(3,2) = s;
  r(2,3) = -s;
  r(3,3) = c;

  return r;
}


////////////////////////////////////////////////////////////////////////////////
//
// Global operators
//
////////////////////////////////////////////////////////////////////////////////

FaMat33 operator- (const FaMat33& a)
{
  return FaMat33(-a.v[0], -a.v[1], -a.v[2]);
}


FaMat33 operator+ (const FaMat33& a, const FaMat33& b)
{
  return FaMat33(a.v[0] + b.v[0], a.v[1] + b.v[1], a.v[2] + b.v[2]);
}


FaMat33 operator- (const FaMat33& a, const FaMat33& b)
{
  return FaMat33(a.v[0] - b.v[0], a.v[1] - b.v[1], a.v[2] - b.v[2]);
}


FaMat33 operator* (const FaMat33& a, const FaMat33& b)
{
  return FaMat33(a*b.v[0], a*b.v[1], a*b.v[2]);
}


FaMat33 operator* (const FaMat33& a, double d)
{
  return FaMat33(a.v[0]*d, a.v[1]*d, a.v[2]*d);
}


FaMat33 operator* (double d, const FaMat33& a)
{
  return FaMat33(a.v[0]*d, a.v[1]*d, a.v[2]*d);
}


FaVec3 operator* (const FaMat33& m, const FaVec3& v1)
{
  return m.v[0]*v1[0] + m.v[1]*v1[1] + m.v[2]*v1[2];
}


FaMat33 operator/ (const FaMat33& a, double d)
{
  if (fabs(d) < EPS_ZERO)
  {
#ifdef FFA_DEBUG
    std::cerr <<"FaMat33 operator/(FaMat33&,double): Division by zero "<< d
	      << std::endl;
#endif
    FaVec3 huge_vec(HUGE_VAL, HUGE_VAL, HUGE_VAL);
    return FaMat33(huge_vec, huge_vec, huge_vec);
  }

  return FaMat33(a.v[0] / d, a.v[1] / d, a.v[2] / d);
}


int operator== (const FaMat33& a, const FaMat33& b)
{
  return (a.v[0] == b.v[0]) && (a.v[1] == b.v[1]) && (a.v[2] == b.v[2]);
}

int operator!= (const FaMat33& a, const FaMat33& b)
{
  return !(a == b);
}


std::ostream& operator<< (std::ostream& s, const FaMat33& m)
{
  return s <<'\n'<< m.v[0][0] <<' '<< m.v[1][0] <<' '<< m.v[2][0]
	   <<'\n'<< m.v[0][1] <<' '<< m.v[1][1] <<' '<< m.v[2][1]
	   <<'\n'<< m.v[0][2] <<' '<< m.v[1][2] <<' '<< m.v[2][2];
}

std::istream& operator>> (std::istream& s, FaMat33& m)
{
  FaMat33 m_tmp;
  s >> m_tmp.v[0][0] >> m_tmp.v[1][0] >> m_tmp.v[2][0]
    >> m_tmp.v[0][1] >> m_tmp.v[1][1] >> m_tmp.v[2][1]
    >> m_tmp.v[0][2] >> m_tmp.v[1][2] >> m_tmp.v[2][2];
  if (s) m = m_tmp;
  return s;
}


////////////////////////////////////////////////////////////////////////////////
//
// FORTRAN interface to selected functions
//
////////////////////////////////////////////////////////////////////////////////

SUBROUTINE(ffa_eulerzyx,FFA_EULERZYX) (double* a, double* b, double* angles)
{
  FaVec3 euler(FaMat33(FaMat33(a).transpose() * FaMat33(b)).getEulerZYX());
  angles[0] = euler[0];
  angles[1] = euler[1];
  angles[2] = euler[2];
}


SUBROUTINE(ffa_glbeulerzyx,FFA_GLBEULERZYX) (double* a, double* angles)
{
  FaVec3 euler(FaMat33(a).getEulerZYX());
  angles[0] = euler[0];
  angles[1] = euler[1];
  angles[2] = euler[2];
}
