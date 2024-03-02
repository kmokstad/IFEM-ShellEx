// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFA_TENSOR1_H
#define FFA_TENSOR1_H

#include <iostream>

class FFaTensor2;
class FFaTensor3;


/*!
  Class to represent a 1D simple tensor.
  Tensor layout:
  T = {t11}
*/

class FFaTensor1
{
  double myT;

public:

  // Constructors

  FFaTensor1(double d = 0.0)      { myT = d; }
  FFaTensor1(const FFaTensor1& t) { myT = t.myT; }
  FFaTensor1(const FFaTensor2& t);
  FFaTensor1(const FFaTensor3& t);

  // Local operators

  FFaTensor1& operator= (const FFaTensor1& t);
  FFaTensor1& operator= (const FFaTensor2& t);
  FFaTensor1& operator= (const FFaTensor3& t);

  FFaTensor1& operator+= (const FFaTensor1& t)
  {
    myT += t.myT;
    return *this;
  }

  FFaTensor1& operator-= (const FFaTensor1& t)
  {
    myT -= t.myT;
    return *this;
  }

  FFaTensor1& operator*= (double d)
  {
    myT *= d;
    return *this;
  }

  FFaTensor1& operator/= (double d)
  {
    myT /= d;
    return *this;
  }

  operator double() const { return myT; }

  // Global operators

  friend FFaTensor1 operator- (const FFaTensor1&);
  friend FFaTensor1 operator+ (const FFaTensor1&, const FFaTensor1&);
  friend FFaTensor1 operator- (const FFaTensor1&, const FFaTensor1&);
  // Scaling by scalar
  friend FFaTensor1 operator* (const FFaTensor1&, double);
  friend FFaTensor1 operator* (double, const FFaTensor1&);
  friend FFaTensor1 operator/ (const FFaTensor1&, double);
  // Equality
  friend bool operator==      (const FFaTensor1&, const FFaTensor1&);
  friend bool operator!=      (const FFaTensor1&, const FFaTensor1&);
  // Writing and reading
  friend std::ostream& operator<< (std::ostream&, const FFaTensor1&);
  friend std::istream& operator>> (std::istream&, FFaTensor1&);
};

#endif
