// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFaLib/FFaAlgebra/FFaMat34.H"


/*!
  Matrix layout :

  [0][0]  [1][0]  [2][0]  [3][0]

  [0][1]  [1][1]  [2][1]  [3][1]

  [0][2]  [1][2]  [2][2]  [3][2]
*/


////////////////////////////////////////////////////////////////////////////////
//
// Local operators
//
////////////////////////////////////////////////////////////////////////////////

FaMat34& FaMat34::operator+= (const FaMat34& m)
{
  r += m.r;
  p += m.p;
  return *this;
}

FaMat34& FaMat34::operator+= (const FaVec3& v1)
{
  p += v1;
  return *this;
}


FaMat34& FaMat34::operator-= (const FaMat34& m)
{
  r -= m.r;
  p -= m.p;
  return *this;
}

FaMat34& FaMat34::operator-= (const FaVec3& v1)
{
  p -= v1;
  return *this;
}


FaMat34& FaMat34::operator*= (double d)
{
  r *= d;
  p *= d;
  return *this;
}


FaMat34& FaMat34::operator/= (double d)
{
  r /= d;
  p /= d;
  return *this;
}


////////////////////////////////////////////////////////////////////////////////
//
// Special functions
//
////////////////////////////////////////////////////////////////////////////////

FaMat34 FaMat34::inverse() const
{
  FaMat33 dirMat = r.transpose();
  return FaMat34(dirMat,-(dirMat*p));
}


FaMat34& FaMat34::setIdentity()
{
  r.setIdentity();
  p.clear();
  return *this;
}


bool FaMat34::isCoincident(const FaMat34& oMat, double tolerance) const
{
  if (p.equals(oMat.p,tolerance))
    return r.isCoincident(oMat.r,tolerance);
  else
    return false;
}


FaMat34& FaMat34::eulerRotateZYX(const FaVec3& anglesXYZ)
{
  r.eulerRotateZYX(anglesXYZ);
  return *this;
}

FaMat34& FaMat34::eulerRotateZYX(const FaVec3& anglesXYZ,
                                 const FaMat34& fromMatrix)
{
  r = fromMatrix.r * r.eulerRotateZYX(anglesXYZ);
  return *this;
}


FaMat34& FaMat34::eulerTransform(const FaVec3& offsetXYZ,
                                 const FaVec3& anglesXYZ,
                                 const FaMat34& fromMatrix)
{
  r = fromMatrix.r * r.eulerRotateZYX(anglesXYZ);
  p = fromMatrix * offsetXYZ;
  return *this;
}


FaMat34& FaMat34::quatrTransform(const FaVec3& offsetXYZ,
                                 const FaVec3& anglesXYZ,
                                 const FaMat34& fromMatrix)
{
  r = fromMatrix.r * r.incRotate(anglesXYZ);
  p = fromMatrix * offsetXYZ;
  return *this;
}


FaMat34& FaMat34::makeGlobalizedCS(const FaVec3& origin,
                                   const FaVec3& p1)
{
  r.makeGlobalizedCS(p1-origin);
  p = origin;
  return *this;
}

FaMat34& FaMat34::makeGlobalizedCS(const FaVec3& origin,
                                   const FaVec3& p1, const FaVec3& p2)
{
  r.makeGlobalizedCS(origin,p1,p2);
  p = origin;
  return *this;
}

FaMat34& FaMat34::makeGlobalizedCS(const FaVec3& origin,
                                   const FaVec3& p1, const FaVec3& p2,
                                   const FaVec3& p3)
{
  r.makeGlobalizedCS(origin,p1,p2,p3);
  p = origin;
  return *this;
}


FaVec3 FaMat34::projectOnXY(const FaVec3& x) const
{
  const FaVec3& eZ = r[VZ];
  return x - (eZ*(x-p))*eZ;
}


FaVec3 FaMat34::getEulerZYX() const
{
  return r.getEulerZYX();
}


double FaMat34::getEulerZYX(int i) const
{
  return r.getEulerZYX()[i];
}

double FaMat34::getEulerZYX(int i, const FaMat34& from) const
{
  return (from.r.transpose()*r).getEulerZYX()[i];
}


double FaMat34::getRotation(int i) const
{
  return r.getRotation()[i];
}

double FaMat34::getRotation(int i, const FaMat34& from) const
{
  return (from.r.transpose()*r).getRotation()[i];
}


FaVec3 FaMat34::getEulerZYX(const FaMat34& from, const FaMat34& to)
{
  return (from.r.transpose()*to.r).getEulerZYX();
}


FaMat34 FaMat34::makeZrotation(double rot)
{
  return FaMat34(FaMat33::makeZrotation(rot),FaVec3());
}

FaMat34 FaMat34::makeYrotation(double rot)
{
  return FaMat34(FaMat33::makeYrotation(rot),FaVec3());
}

FaMat34 FaMat34::makeXrotation(double rot)
{
  return FaMat34(FaMat33::makeXrotation(rot),FaVec3());
}


FaMat34& FaMat34::makeCS_X_XY(const FaVec3& origin,
                              const FaVec3& xpt, const FaVec3& xypl)
{
  r[VX] = xpt - origin;
  r[VZ] = r[VX] ^ (xypl - origin);
  r[VX].normalize();
  r[VZ].normalize();
  r[VY] = r[VZ] ^ r[VX];
  p = origin;
  return *this;
}


FaMat34& FaMat34::makeCS_X_XZ(const FaVec3& origin,
                              const FaVec3& xpt, const FaVec3& xzpl)
{
  r[VX] = xpt - origin;
  r[VY] = (xzpl - origin) ^ r[VX];
  r[VX].normalize();
  r[VY].normalize();
  r[VZ] = r[VX] ^ r[VY];
  p = origin;
  return *this;
}


FaMat34& FaMat34::makeCS_Z_XZ(const FaVec3& origin,
                              const FaVec3& zpt, const FaVec3& xzpl)
{
  r[VZ] = zpt - origin;
  r[VY] = r[VZ] ^ (xzpl - origin);
  r[VZ].normalize();
  r[VY].normalize();
  r[VX] = r[VY] ^ r[VZ];
  p = origin;
  return *this;
}


////////////////////////////////////////////////////////////////////////////////
//
// Reading and writing
//
////////////////////////////////////////////////////////////////////////////////

/*!
  Writes a reduced homogenous transformation matrix. Format:
  R1x  R2x  R3x  Tx
  R1y  R2y  R3y  Ty
  R1z  R2z  R3z  Tz
*/

std::ostream& FaMat34::printStd(std::ostream& os) const
{
  std::ios_base::fmtflags tmpflags = os.flags(std::ios::fixed);
  int tmpprec  = os.precision(8);

  os <<'\n'<< r[VX][VX] <<' '<< r[VY][VX] <<' '<< r[VZ][VX] <<' '<< p[VX]
     <<'\n'<< r[VX][VY] <<' '<< r[VY][VY] <<' '<< r[VZ][VY] <<' '<< p[VY]
     <<'\n'<< r[VX][VZ] <<' '<< r[VY][VZ] <<' '<< r[VZ][VZ] <<' '<< p[VZ];

  os.flags(tmpflags);
  os.precision(tmpprec);
  return os;
}


/*!
  Writes a reduced homogenous transformation matrix. Format:
  R1x  R1y  R1z
  R2x  R2y  R3z
  R3x  R3y  R3z
  Tx   Ty   Tz
*/

std::ostream& FaMat34::printRot(std::ostream& os) const
{
  return os <<'\n'<< r[VX] <<'\n'<< r[VY] <<'\n'<< r[VZ] <<'\n'<< p;
}


/*!
  Reads a reduced homogenous transformation matrix. Format:
  R1x  R2x  R3x  Tx
  R1y  R2y  R3y  Ty
  R1z  R2z  R3z  Tz
*/

bool FaMat34::readStd(std::istream& s)
{
  s >> r[VX][VX] >> r[VY][VX] >> r[VZ][VX] >> p[VX];
  s >> r[VX][VY] >> r[VY][VY] >> r[VZ][VY] >> p[VY];
  s >> r[VX][VZ] >> r[VY][VZ] >> r[VZ][VZ] >> p[VZ];
  return s ? true : false;
}


/*!
  Reads a reduced homogenous transformation matrix. Format:
  R1x  R1y  R1z
  R2x  R2y  R3z
  R3x  R3y  R3z
  Tx   Ty   Tz
*/

bool FaMat34::readRot(std::istream& s)
{
  s >> r[VX] >> r[VY] >> r[VZ] >> p;
  return s ? true : false;
}


////////////////////////////////////////////////////////////////////////////////
//
// Global operators
//
////////////////////////////////////////////////////////////////////////////////

FaMat34 operator- (const FaMat34& a)
{
  return FaMat34(-a.r, -a.p);
}


FaMat34 operator+ (const FaMat34& a, const FaMat34& b)
{
  return FaMat34(a.r + b.r, a.p + b.p);
}

FaMat34 operator- (const FaMat34& a, const FaMat34& b)
{
  return FaMat34(a.r - b.r, a.p - b.p);
}


FaMat34 operator* (const FaMat34& a, const FaMat34& b)
{
  return FaMat34(a.r*b.r, a.r*b.p + a.p);
}


FaMat34 operator* (const FaMat34& a, const FaMat33& b)
{
  return FaMat34(a.r*b, a.p);
}


FaVec3 operator* (const FaMat34& a, const FaVec3& b)
{
  return FaVec3(a.r*b + a.p);
}


bool operator== (const FaMat34& a, const FaMat34& b)
{
  return a.p == b.p && a.r == b.r;
}

bool operator!= (const FaMat34& a, const FaMat34& b)
{
  return !(a == b);
}


std::ostream& operator<< (std::ostream& s, const FaMat34& m)
{
  return m.printStd(s);
}

std::istream& operator>> (std::istream& s, FaMat34& m)
{
  FaMat34 m_tmp;
  if (m_tmp.readStd(s)) m = m_tmp;
  return s;
}
