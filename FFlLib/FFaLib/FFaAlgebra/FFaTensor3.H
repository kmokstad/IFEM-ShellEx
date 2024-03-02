// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFA_TENSOR3_H
#define FFA_TENSOR3_H

#include <iostream>

class FaVec3;
class FaMat33;
class FaMat34;
class FFaTensor1;
class FFaTensor2;


/*!
  Class to represent a 3D symmetrical tensor.
  Tensor layout:
  T = {t11, t22, t33, t12, t13, t23}
*/

class FFaTensor3
{
  double myT[6];

public:

  // Constructors

  FFaTensor3()                    { for (int i=0; i<6; i++) myT[i] = 0.0; }
  FFaTensor3(double d)            { for (int i=0; i<6; i++) myT[i] = d; }
  FFaTensor3(const float* t)      { for (int i=0; i<6; i++) myT[i] = t[i]; }
  FFaTensor3(const double* t)     { for (int i=0; i<6; i++) myT[i] = t[i]; }
  FFaTensor3(const FFaTensor3& t) { for (int i=0; i<6; i++) myT[i] = t.myT[i]; }
  FFaTensor3(const FFaTensor2& t);
  FFaTensor3(const FFaTensor1& t);
  FFaTensor3(const FaVec3& v);
  FFaTensor3(double t11      , double t22      , double t33,
             double t12 = 0.0, double t13 = 0.0, double t23 = 0.0)
  {
    myT[0] = t11; myT[1] = t22; myT[2] = t33;
    myT[3] = t12; myT[4] = t13; myT[5] = t23;
  }
  FFaTensor3(const FaVec3& v1, const FaVec3& v2, const FaVec3& v3)
  {
    this->makeInertia(v1,v2,v3);
  }

  // Local operators

  FFaTensor3& operator= (const FFaTensor3& t);
  FFaTensor3& operator= (const FFaTensor1& t);
  FFaTensor3& operator= (const FFaTensor2& t);

  FFaTensor3& operator+= (const FFaTensor3& t)
  {
    for (int i=0; i<6; i++) myT[i] += t.myT[i];
    return *this;
  }

  FFaTensor3& operator-= (const FFaTensor3& t)
  {
    for (int i=0; i<6; i++) myT[i] -= t.myT[i];
    return *this;
  }

  FFaTensor3& operator*= (double d)
  {
    for (int i=0; i<6; i++) myT[i] *= d;
    return *this;
  }

  FFaTensor3& operator/= (double d)
  {
    for (int i=0; i<6; i++) myT[i] /= d;
    return *this;
  }

  // Access to internal array

  const double* getPt() const { return myT; }
  double*       getPt()       { return myT; }

  // Indexing

  const double& operator[] (int i) const;
  double&       operator[] (int i);

  // Special functions

  FFaTensor3& rotate(const FaMat33& rotMx);
  FFaTensor3& rotate(const FaMat34& rotMx);

  FFaTensor3& makeInertia(const FaVec3& v1, const FaVec3& v2, const FaVec3& v3);
  FFaTensor3& translateInertia(const FaVec3& x, double mass);

  double vonMises() const;
  double maxShear() const;
  void   maxShear(FaVec3& v) const;
  double maxPrinsipal(bool absMax = false) const;
  double middlePrinsipal() const;
  double minPrinsipal() const;
  void   prinsipalValues(double& max, double& middle, double& min) const;
  void   prinsipalValues(FaVec3& values, FaMat33& rotation) const;

  // Global operators

  friend FFaTensor3 operator- (const FFaTensor3&);
  friend FFaTensor3 operator+ (const FFaTensor3&, const FFaTensor3&);
  friend FFaTensor3 operator- (const FFaTensor3&, const FFaTensor3&);
  // Scaling by scalar
  friend FFaTensor3 operator* (const FFaTensor3&, double);
  friend FFaTensor3 operator* (double, const FFaTensor3&);
  friend FFaTensor3 operator/ (const FFaTensor3&, double);
  // Equality
  friend bool      operator== (const FFaTensor3&, const FFaTensor3&);
  friend bool      operator!= (const FFaTensor3&, const FFaTensor3&);
  // Rotation
  friend FFaTensor3 operator* (const FFaTensor3&, const FaMat33&);
  friend FFaTensor3 operator* (const FFaTensor3&, const FaMat34&);
  friend FFaTensor3 operator* (const FaMat33&   , const FFaTensor3&);
  friend FFaTensor3 operator* (const FaMat34&   , const FFaTensor3&);
  // Writing and reading
  friend std::ostream& operator<< (std::ostream&, const FFaTensor3&);
  friend std::istream& operator>> (std::istream&, FFaTensor3&);
};


// --- inline functions ---

inline const double& FFaTensor3::operator[] (int i) const
{
#ifdef FFA_INDEXCHECK
  if (i < 0 || i > 5)
    std::cerr <<"FFaTensor3::operator[]: index i="<< i <<" is out of range [0,5]"
	      << std::endl;
#endif
  return myT[i];
}

inline double& FFaTensor3::operator[] (int i)
{
#ifdef FFA_INDEXCHECK
  if (i < 0 || i > 5)
    std::cerr <<"FFaTensor3::operator[]: index i="<< i <<" is out of range [0,5]"
	      << std::endl;
#endif
  return myT[i];
}

#endif
