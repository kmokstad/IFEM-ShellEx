// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFA_TENSOR2_H
#define FFA_TENSOR2_H

#include <iostream>

class FaVec3;
class FaMat33;
class FaMat34;
class FFaTensor1;
class FFaTensor3;


/*!
  Class to represent a 2D symmetrical tensor.
  Tensor layout:
  T = {t11, t22, t12}
*/

class FFaTensor2
{
  double myT[3];

public:

  // Constructors

  FFaTensor2()                    { myT[0] = myT[1] = myT[2] = 0.0; }
  FFaTensor2(double d)            { myT[0] = myT[1] = myT[2] = d; }
  FFaTensor2(const float* t)      { for (int i=0; i<3; i++) myT[i] = t[i]; }
  FFaTensor2(const double* t)     { for (int i=0; i<3; i++) myT[i] = t[i]; }
  FFaTensor2(const FFaTensor2& t) { for (int i=0; i<3; i++) myT[i] = t.myT[i]; }
  FFaTensor2(const FFaTensor1& t);
  FFaTensor2(const FFaTensor3& t);
  FFaTensor2(double t11, double t22, double t12 = 0.0)
  { myT[0] = t11; myT[1] = t22; myT[2] = t12; }

  // Local operators

  FFaTensor2& operator= (const FFaTensor1& t);
  FFaTensor2& operator= (const FFaTensor2& t);
  FFaTensor2& operator= (const FFaTensor3& t);

  FFaTensor2& operator+= (const FFaTensor2& t)
  {
    for (int i=0; i<3; i++) myT[i] += t.myT[i];
    return *this;
  }

  FFaTensor2& operator-= (const FFaTensor2& t)
  {
    for (int i=0; i<3; i++) myT[i] -= t.myT[i];
    return *this;
  }

  FFaTensor2& operator*= (double d)
  {
    for (int i=0; i<3; i++) myT[i] *= d;
    return *this;
  }

  FFaTensor2& operator/= (double d)
  {
    for (int i=0; i<3; i++) myT[i] /= d;
    return *this;
  }

  // Access to internal array

  const double* getPt() const { return myT; }
  double*       getPt()       { return myT; }

  // Indexing

  const double& operator[] (int i) const;
  double&       operator[] (int i);

  // Special functions

  FFaTensor2& rotate(const double ex[2], const double ey[2]);

  double vonMises() const;
  double maxShear() const;
  void   maxShear(FaVec3& v) const;
  double maxPrinsipal(bool absMax = false) const;
  double minPrinsipal() const;
  void   prinsipalValues(double& max, double& min) const;
  void   prinsipalValues(FaVec3& values, FaMat33& rotation) const;

  // Global operators

  friend FFaTensor2 operator- (const FFaTensor2&);
  friend FFaTensor2 operator+ (const FFaTensor2&, const FFaTensor2&);
  friend FFaTensor2 operator- (const FFaTensor2&, const FFaTensor2&);
  // Scaling by scalar
  friend FFaTensor2 operator* (const FFaTensor2&, double);
  friend FFaTensor2 operator* (double, const FFaTensor2&);
  friend FFaTensor2 operator/ (const FFaTensor2&, double);
  // Equality
  friend bool      operator== (const FFaTensor2&, const FFaTensor2&);
  friend bool      operator!= (const FFaTensor2&, const FFaTensor2&);
  // Rotation
  friend FFaTensor3 operator* (const FFaTensor2&, const FaMat33&);
  friend FFaTensor3 operator* (const FFaTensor2&, const FaMat34&);
  friend FFaTensor3 operator* (const FaMat33&   , const FFaTensor2&);
  friend FFaTensor3 operator* (const FaMat34&   , const FFaTensor2&);
  // Writing and reading
  friend std::ostream& operator<< (std::ostream& s, const FFaTensor2& b);
  friend std::istream& operator>> (std::istream& s, FFaTensor2& b);
};


// --- inline functions ---

inline const double& FFaTensor2::operator[] (int i) const
{
#ifdef FFA_INDEXCHECK
  if (i < 0 || i > 2)
    std::cerr <<"FFaTensor2::operator[]: index i="<< i <<" is out of range [0,2]"
	      << std::endl;
#endif
  return myT[i];
}

inline double& FFaTensor2::operator[] (int i)
{
#ifdef FFA_INDEXCHECK
  if (i < 0 || i > 2)
    std::cerr <<"FFaTensor2::operator[]: index i="<< i <<" is out of range [0,2]"
	      << std::endl;
#endif
  return myT[i];
}

#endif
