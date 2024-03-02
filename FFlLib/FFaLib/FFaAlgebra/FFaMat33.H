// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFA_MAT33_H
#define FFA_MAT33_H

#include "FFaLib/FFaAlgebra/FFaVec3.H"


class FaMat33
{
  FaVec3 v[3];

public:

  // Constructors

  FaMat33() { this->setIdentity(); }
  FaMat33(const float* mat);
  FaMat33(const double* mat);
  FaMat33(const FaVec3& v0, const FaVec3& v1, const FaVec3& v2);
  FaMat33(const FaMat33& m) { for (int i = 0; i < 3; i++) v[i] = m.v[i]; }

  // Local operators

  FaMat33& operator=  (const FaMat33& m);
  FaMat33& operator+= (const FaMat33& m);
  FaMat33& operator-= (const FaMat33& m);
  FaMat33& operator*= (double d);
  FaMat33& operator/= (double d);

  // Indexing

  const FaVec3& operator[] (int i) const;
  FaVec3&       operator[] (int i);

  const double& operator() (int i, int j) const;
  double&       operator() (int i, int j);

  // Special functions

  FaMat33  inverse(double eps = 1.0e-16) const;
  FaMat33& setIdentity();
  FaMat33  transpose() const;
  FaMat33& shift(int delta);

  bool isCoincident(const FaMat33& m, double tolerance = 1.0e-10) const;

  FaMat33& makeGlobalizedCS(const FaVec3& v1);
  FaMat33& makeGlobalizedCS(const FaVec3& v1, const FaVec3& v2);
  FaMat33& makeGlobalizedCS(const FaVec3& v0,
                            const FaVec3& v1, const FaVec3& v2);
  FaMat33& makeGlobalizedCS(const FaVec3& v1, const FaVec3& v2,
                            const FaVec3& v3, const FaVec3& v4);

  FaMat33& eulerRotateZYX(const FaVec3& angles);
  FaMat33& incRotate(const FaVec3& angles);

  FaVec3 getEulerZYX() const;
  FaVec3 getRotation() const;

  static FaMat33 makeZrotation(double angle);
  static FaMat33 makeYrotation(double angle);
  static FaMat33 makeXrotation(double angle);

  // Global operators

  friend FaMat33  operator-  (const FaMat33& a);
  friend FaMat33  operator+  (const FaMat33& a, const FaMat33& b);
  friend FaMat33  operator-  (const FaMat33& a, const FaMat33& b);

  friend FaMat33  operator*  (const FaMat33& a, const FaMat33& b);
  friend FaMat33  operator*  (const FaMat33& a, double d);
  friend FaMat33  operator*  (double d, const FaMat33& a);
  friend FaVec3   operator*  (const FaMat33& a, const FaVec3&  b);
  friend FaMat33  operator/  (const FaMat33& a, double d);

  friend int      operator== (const FaMat33& a, const FaMat33& b);
  friend int      operator!= (const FaMat33& a, const FaMat33& b);

  friend std::ostream& operator<< (std::ostream& s, const FaMat33& m);
  friend std::istream& operator>> (std::istream& s, FaMat33& m);
};


// --- inline functions ---

inline const FaVec3& FaMat33::operator[] (int i) const
{
#ifdef FFA_INDEXCHECK
  if (i < 0 || i > 2)
    std::cerr <<"FaMat33::operator[]: index i="<< i <<" is out of range [0,2]"
              << std::endl;
#endif
  return v[i];
}

inline FaVec3& FaMat33::operator[] (int i)
{
#ifdef FFA_INDEXCHECK
  if (i < 0 || i > 2)
    std::cerr <<"FaMat33::operator[]: index i="<< i <<" is out of range [0,2]"
              << std::endl;
#endif
  return v[i];
}

// Fortran-like 1-based indexing operators.
// The first index is the row number, whereas the second is the column number.

inline const double& FaMat33::operator() (int i, int j) const
{
#ifdef FFA_INDEXCHECK
  if (i < 1 || i > 3)
    std::cerr <<"FaMat33::operator(): index i="<< i <<" is out of range [1,3]"
              << std::endl;
  if (j < 1 || j > 3)
    std::cerr <<"FaMat33::operator(): index j="<< j <<" is out of range [1,3]"
              << std::endl;
#endif
  return v[j-1][i-1];
}

inline double& FaMat33::operator() (int i, int j)
{
#ifdef FFA_INDEXCHECK
   if (i < 1 || i > 3)
     std::cerr <<"FaMat33::operator(): index i="<< i <<" is out of range [1,3]"
               << std::endl;
   if (j < 1 || j > 3)
     std::cerr <<"FaMat33::operator(): index j="<< j <<" is out of range [1,3]"
               << std::endl;
#endif
  return v[j-1][i-1];
}

#endif
