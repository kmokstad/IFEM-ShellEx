// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFaLib/FFaAlgebra/FFaTensor1.H"
#include "FFaLib/FFaAlgebra/FFaTensor2.H"
#include "FFaLib/FFaAlgebra/FFaTensor3.H"
#include <cctype>
#include <cmath>


FFaTensor1::FFaTensor1(const FFaTensor2& t)
{
  myT = t[0];
}


FFaTensor1::FFaTensor1(const FFaTensor3& t)
{
  myT = t[0];
}


FFaTensor1& FFaTensor1::operator= (const FFaTensor1& t)
{
  if (this != &t) myT = t.myT;
  return *this;
}

FFaTensor1& FFaTensor1::operator= (const FFaTensor2& t)
{
  myT = t[0];
  return *this;
}


FFaTensor1& FFaTensor1::operator= (const FFaTensor3& t)
{
  myT = t[0];
  return *this;
}


/*!
  Global operators.
*/

FFaTensor1 operator- (const FFaTensor1& t)
{
  return FFaTensor1(-t.myT);
}


FFaTensor1 operator+ (const FFaTensor1& a, const FFaTensor1& b)
{
  return FFaTensor1(a.myT + b.myT);
}


FFaTensor1 operator- (const FFaTensor1& a, const FFaTensor1& b)
{
  return FFaTensor1(a.myT - b.myT);
}


FFaTensor1 operator* (const FFaTensor1& a, double d)
{
  return FFaTensor1(a.myT * d);
}

FFaTensor1 operator* (double d, const FFaTensor1& a)
{
  return FFaTensor1(a.myT * d);
}


FFaTensor1 operator/ (const FFaTensor1& a, double d)
{
  if (fabs(d) < 1.0e-16)
    return FFaTensor1(HUGE_VAL);

  return FFaTensor1(a.myT / d);
}


bool operator== (const FFaTensor1& a, const FFaTensor1& b)
{
  return (a.myT == b.myT);
}

bool operator!= (const FFaTensor1& a, const FFaTensor1& b)
{
  return (a.myT != b.myT);
}


std::ostream& operator<< (std::ostream& s, const FFaTensor1& t)
{
  s << t.myT;
  return s;
}

std::istream& operator>> (std::istream& s, FFaTensor1& t)
{
  FFaTensor1 tmpT;
  s >> tmpT.myT;
  if (s) t = tmpT;
  return s;
}
