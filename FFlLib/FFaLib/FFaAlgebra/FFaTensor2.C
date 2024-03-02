// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFaLib/FFaAlgebra/FFaTensor2.H"
#include "FFaLib/FFaAlgebra/FFaTensor3.H"
#include "FFaLib/FFaAlgebra/FFaTensor1.H"
#include "FFaLib/FFaAlgebra/FFaMat33.H"
#include "FFaLib/FFaAlgebra/FFaMat34.H"
#include "FFaLib/FFaAlgebra/FFaTensorTransforms.H"
#include <cctype>


FFaTensor2::FFaTensor2(const FFaTensor3& t)
{
  myT[0] = t[0];
  myT[1] = t[1];
  myT[2] = t[3];
}


FFaTensor2::FFaTensor2(const FFaTensor1& t)
{
  myT[0] = t;
  myT[1] = myT[2] = 0.0;
}


FFaTensor2& FFaTensor2::operator= (const FFaTensor3& t)
{
  myT[0] = t[0];
  myT[1] = t[1];
  myT[2] = t[3];
  return *this;
}


FFaTensor2& FFaTensor2::operator= (const FFaTensor2& t)
{
  if (this != &t)
  {
    myT[0] = t.myT[0];
    myT[1] = t.myT[1];
    myT[2] = t.myT[2];
  }
  return *this;
}


FFaTensor2& FFaTensor2::operator= (const FFaTensor1& t)
{
  myT[0] = t;
  myT[1] = myT[2] = 0.0;
  return *this;
}


/*!
  Rotate the tensor to the given CS.
*/

FFaTensor2& FFaTensor2::rotate(const double ex[2], const double ey[2])
{
  FFaTensorTransforms::rotate(myT,ex,ey, myT);
  return *this;
}


/*!
  Get the von Mises value of the tensor.
*/

double FFaTensor2::vonMises() const
{
  return FFaTensorTransforms::vonMises(myT[0],myT[1],myT[2]);
}


/*!
  Get the max shear value of the tensor.
  If it can't be found, HUGE_VAL will be returned.
*/

double FFaTensor2::maxShear() const
{
  double max, min;
  if (FFaTensorTransforms::principalValues(myT[0],myT[1],myT[2],max,min))
    return FFaTensorTransforms::maxShearValue(max,min);
  else
    return HUGE_VAL;
}


/*!
  Get the max shear of the tensor as a directed vector.
  If it can't be found, a 0-vector will be returned.
*/

void FFaTensor2::maxShear(FaVec3& v) const
{
  double values[2], max[2], min[2];
  if (FFaTensorTransforms::principalDirs(myT,values,max,min) == 0)
  {
    v[2] = 0.0;
    FFaTensorTransforms::maxShearDir(2,max,min,v.getPt());
    v *= FFaTensorTransforms::maxShearValue(values[0],values[1]);
  }
  else
    v.clear();
}


/*!
  Get the (absolute) max principal of the tensor.
  If it can't be found, HUGE_VAL will be returned.
*/

double FFaTensor2::maxPrinsipal(bool absMax) const
{
  double max, min;
  if (FFaTensorTransforms::principalValues(myT[0],myT[1],myT[2],max,min))
    return absMax ? (fabs(max) > fabs(min) ? max : min) : max;
  else
    return HUGE_VAL;
}


/*!
  Get the min principal of the tensor.
  If it can't be found, HUGE_VAL will be returned.
*/

double FFaTensor2::minPrinsipal() const
{
  double max, min;
  if (FFaTensorTransforms::principalValues(myT[0],myT[1],myT[2],max,min))
    return min;
  else
    return HUGE_VAL;
}


/*!
  Get the principal values of the tensor.
  If they can't be found, HUGE_VAL will be returned.
*/

void FFaTensor2::prinsipalValues(double& max, double& min) const
{
  if (!FFaTensorTransforms::principalValues(myT[0],myT[1],myT[2],max,min))
    max = min = HUGE_VAL;
}


/*!
  Get a valid rotation matrix corresponding to the principal axes of the tensor.
  The associated principal values are also found in the corresponding order.
  If the matrix can't be found, the identity will be returned,
  along with a vector of HUGE_VAL.
*/

void FFaTensor2::prinsipalValues(FaVec3& values, FaMat33& rotation) const
{
  values[2] = 0.0;
  rotation.setIdentity();
  if (FFaTensorTransforms::principalDirs(myT,
                                         values.getPt(),
                                         rotation[0].getPt(),
                                         rotation[1].getPt()))
  {
    values[0] = values[1] = HUGE_VAL;
    values[2] = 0.0;
  }
  else if ((rotation[0] ^ rotation[1]).isParallell(rotation[2]) != 1)
  {
    // Swap the 1st and 2nd
    std::swap(rotation[0],rotation[1]);
    std::swap(values[0],values[1]);
  }
}


/*!
  Global operators.
*/

FFaTensor2 operator- (const FFaTensor2& t)
{
  return FFaTensor2(-t.myT[0], -t.myT[1], -t.myT[2]);
}


FFaTensor2 operator+ (const FFaTensor2& a, const FFaTensor2& b)
{
  return FFaTensor2(a.myT[0] + b.myT[0],
                    a.myT[1] + b.myT[1],
                    a.myT[2] + b.myT[2]);
}

FFaTensor2 operator- (const FFaTensor2& a, const FFaTensor2& b)
{
  return FFaTensor2(a.myT[0] - b.myT[0],
                    a.myT[1] - b.myT[1],
                    a.myT[2] - b.myT[2]);
}


FFaTensor2 operator* (const FFaTensor2& a, double d)
{
  return FFaTensor2(a.myT[0]*d, a.myT[1]*d, a.myT[2]*d);
}

FFaTensor2 operator* (double d, const FFaTensor2& a)
{
  return FFaTensor2(a.myT[0]*d, a.myT[1]*d, a.myT[2]*d);
}


FFaTensor2 operator/ (const FFaTensor2& a, double d)
{
  if (fabs(d) < 1.0e-16)
    return FFaTensor2(HUGE_VAL);

  return FFaTensor2(a.myT[0]/d, a.myT[1]/d, a.myT[2]/d);
}


bool operator== (const FFaTensor2& a, const FFaTensor2& b)
{
  return (a.myT[0] == b.myT[0] && a.myT[1] == b.myT[1] && a.myT[2] == b.myT[2]);
}


bool operator!= (const FFaTensor2& a, const FFaTensor2& b)
{
  return !(a==b);
}


FFaTensor3 operator* (const FFaTensor2& a, const FaMat33 & m)
{
  FFaTensor3 result;
  FFaTensor3 in(a);
  FFaTensorTransforms::rotate(in.getPt(),
			      m[0].getPt(), m[1].getPt(), m[2].getPt(),
                              result.getPt());
  return result;
}

FFaTensor3 operator* (const FFaTensor2& a, const FaMat34& m)
{
  FFaTensor3 result;
  FFaTensor3 in(a);
  FFaTensorTransforms::rotate(in.getPt(),
			      m[0].getPt(), m[1].getPt(), m[2].getPt(),
                              result.getPt());
  return result;
}


FFaTensor3 operator* (const FaMat33& m, const FFaTensor2& a)
{
  return a*m;
}

FFaTensor3 operator* (const FaMat34& m, const FFaTensor2& a)
{
  return a*m;
}


std::ostream& operator<< (std::ostream& s, const FFaTensor2& t)
{
  s << t.myT[0] << ' ' << t.myT[1] << ' ' << t.myT[2];
  return s;
}

std::istream& operator>> (std::istream& s, FFaTensor2& t)
{
  FFaTensor2 tmpT;
  s >> tmpT[0] >> tmpT[1] >> tmpT[2];
  if (s) t = tmpT;
  return s;
}
