// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFA_MAT34_H
#define FFA_MAT34_H

#include "FFaLib/FFaAlgebra/FFaMat33.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"


class FaMat34
{
  FaMat33 r;
  FaVec3  p;

public:

  // Constructors

  FaMat34() { this->setIdentity(); }
  FaMat34(const float*  mat) : r(mat), p(mat+9) {}
  FaMat34(const double* mat) : r(mat), p(mat+9) {}
  FaMat34(const FaVec3& pos) : p(pos) {}
  FaMat34(const FaVec3& v0, const FaVec3& v1,
          const FaVec3& v2, const FaVec3& v3) : r(v0,v1,v2), p(v3) {}
  FaMat34(const FaMat33& m, const FaVec3& v) : r(m), p(v) {}
  FaMat34(const FaMat34& m) : r(m.r), p(m.p) {}

  // Local operators

  FaMat34& operator= (const FaMat34& m) { r = m.r; p = m.p; return *this; }
  FaMat34& operator= (const FaVec3&  v) { p = v; return *this; }

  FaMat34& operator+= (const FaMat34& m);
  FaMat34& operator+= (const FaVec3&  v);
  FaMat34& operator-= (const FaMat34& m);
  FaMat34& operator-= (const FaVec3&  v);
  FaMat34& operator*= (double d);
  FaMat34& operator/= (double d);

  // Indexing

  const FaVec3& operator[] (int i) const;
  FaVec3&       operator[] (int i);

  const double& operator() (int i, int j) const;
  double&       operator() (int i, int j);

  const FaMat33& direction() const { return r; }
  const FaVec3& translation() const { return p; }

  // Special functions

  FaMat34  inverse() const;
  FaMat34& setIdentity();

  bool isCoincident(const FaMat34& m, double tolerance = 1.0e-10) const;

  FaMat34& eulerRotateZYX(const FaVec3& anglesXYZ, const FaMat34& fromMatrix);
  FaMat34& eulerRotateZYX(const FaVec3& anglesXYZ);
  FaMat34& eulerTransform(const FaVec3& offsetXYZ, const FaVec3& anglesXYZ,
                          const FaMat34& fromMatrix);
  FaMat34& quatrTransform(const FaVec3& offsetXYZ, const FaVec3& anglesXYZ,
                          const FaMat34& fromMatrix);

  FaMat34& makeGlobalizedCS(const FaVec3& origin,
                            const FaVec3& p1);
  FaMat34& makeGlobalizedCS(const FaVec3& origin,
                            const FaVec3& p1, const FaVec3& p2);
  FaMat34& makeGlobalizedCS(const FaVec3& origin,
                            const FaVec3& p1, const FaVec3& p2,
                            const FaVec3& p3);

  FaMat34& makeCS_X_XY(const FaVec3& origin,
                       const FaVec3& xpt, const FaVec3& xypl);
  FaMat34& makeCS_X_XZ(const FaVec3& origin,
                       const FaVec3& xpt, const FaVec3& xzpl);
  FaMat34& makeCS_Z_XZ(const FaVec3& origin,
                       const FaVec3& zpt, const FaVec3& xzpl);

  FaVec3 projectOnXY(const FaVec3& x) const;
  FaVec3 getEulerZYX() const;
  double getEulerZYX(int i) const;
  double getEulerZYX(int i, const FaMat34& from) const;
  double getRotation(int i) const;
  double getRotation(int i, const FaMat34& from) const;

  static FaVec3  getEulerZYX  (const FaMat34& from, const FaMat34& to);
  static FaMat34 makeZrotation(double angle);
  static FaMat34 makeYrotation(double angle);
  static FaMat34 makeXrotation(double angle);

  // Reading and writing

  std::ostream& printStd(std::ostream& os) const;
  std::ostream& printRot(std::ostream& os) const;

  bool readStd(std::istream& s);
  bool readRot(std::istream& s);

  // Global operators

  friend FaMat34  operator-  (const FaMat34& a);
  friend FaMat34  operator+  (const FaMat34& a, const FaMat34& b);
  friend FaMat34  operator+  (const FaMat34& a, const FaMat33& b);
  friend FaMat34  operator+  (const FaMat34& a, const FaVec3&  b);
  friend FaMat34  operator-  (const FaMat34& a, const FaMat34& b);
  friend FaMat34  operator-  (const FaMat34& a, const FaMat33& b);
  friend FaMat34  operator-  (const FaMat34& a, const FaVec3&  b);

  friend FaMat34  operator*  (const FaMat34& a, const FaMat34& b);
  friend FaMat34  operator*  (const FaMat34& a, const FaMat33& b);
  friend FaVec3   operator*  (const FaMat34& a, const FaVec3&  b);

  friend bool     operator== (const FaMat34& a, const FaMat34& b);
  friend bool     operator!= (const FaMat34& a, const FaMat34& b);

  friend std::ostream& operator<< (std::ostream& s, const FaMat34& m);
  friend std::istream& operator>> (std::istream& s, FaMat34& m);
};


// --- inline functions ---

inline const FaVec3& FaMat34::operator[] (int i) const
{
#ifdef FFA_INDEXCHECK
  if (i < 0 || i > 3)
    std::cerr <<"FaMat34::operator[]: index i="<< i <<" is out of range [0,3]"
              << std::endl;
#endif
  return i == 3 ? p : r[i];
}

inline FaVec3& FaMat34::operator[] (int i)
{
#ifdef FFA_INDEXCHECK
  if (i < 0 || i > 3)
    std::cerr <<"FaMat34::operator[]: index i="<< i <<" is out of range [0,3]"
              << std::endl;
#endif
  return i == 3 ? p : r[i];
}

// Fortran-like 1-based indexing operators.
// The first index is the row number, whereas the second is the column number.

inline const double& FaMat34::operator() (int i, int j) const
{
  return j == 4 ? p[i-1] : r[j-1][i-1];
}

inline double& FaMat34::operator() (int i, int j)
{
  return j == 4 ? p[i-1] : r[j-1][i-1];
}

//TODO: Remove the symbol "makeCS_X_YX", use "makeCS_X_XY" instead
#define makeCS_X_YX makeCS_X_XY
#endif
