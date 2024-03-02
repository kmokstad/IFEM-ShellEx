// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFaLib/FFaAlgebra/FFaTensor3.H"
#include "FFaLib/FFaAlgebra/FFaTensor2.H"
#include "FFaLib/FFaAlgebra/FFaTensor1.H"
#include "FFaLib/FFaAlgebra/FFaMat33.H"
#include "FFaLib/FFaAlgebra/FFaMat34.H"
#include "FFaLib/FFaAlgebra/FFaTensorTransforms.H"
#include <cctype>


FFaTensor3::FFaTensor3(const FFaTensor2& t)
{
  myT[0] = t[0];
  myT[1] = t[1];
  myT[3] = t[2];
  myT[2] = myT[4] = myT[5] = 0.0;
}


FFaTensor3::FFaTensor3(const FFaTensor1& t)
{
  myT[0] = t;
  myT[1] = myT[2] = myT[3] = myT[4] = myT[5] = 0.0;
}


FFaTensor3::FFaTensor3(const FaVec3& v)
{
  myT[0] = v[0];
  myT[1] = v[1];
  myT[2] = v[2];
  myT[3] = myT[4] = myT[5] = 0.0;
}


FFaTensor3& FFaTensor3::operator= (const FFaTensor3& t)
{
  if (this != &t)
    for (int i = 0; i < 6; i++)
      myT[i] = t.myT[i];
  return *this;
}

FFaTensor3& FFaTensor3::operator= (const FFaTensor2& t)
{
  myT[0] = t[0];
  myT[1] = t[1];
  myT[3] = t[2];
  myT[2] = myT[4] = myT[5] = 0.0;
  return *this;
}


FFaTensor3& FFaTensor3::operator= (const FFaTensor1& t)
{
  myT[0] = t;
  myT[1] = myT[2] = myT[3] = myT[4] = myT[5] = 0.0;
  return *this;
}


/*!
  Rotate the tensor to the given CS.
*/

FFaTensor3& FFaTensor3::rotate(const FaMat33& rotMx)
{
  FFaTensorTransforms::rotate(myT,
                              rotMx[0].getPt(),
                              rotMx[1].getPt(),
                              rotMx[2].getPt(),
                              myT);
  return *this;
}


FFaTensor3& FFaTensor3::rotate(const FaMat34& rotMx)
{
  FFaTensorTransforms::rotate(myT,
                              rotMx[0].getPt(),
                              rotMx[1].getPt(),
                              rotMx[2].getPt(),
                              myT);
  return *this;
}


/*!
  Create the inertia tensor of a tetrahedron spanned by the origin (0,0,0)
  and the given three points \a v1, \a v2 and \a v3.
*/

FFaTensor3& FFaTensor3::makeInertia(const FaVec3& v1, const FaVec3& v2,
                                    const FaVec3& v3)
{
  double x1 = v1[0], y1 = v1[1], z1 = v1[2];
  double x2 = v2[0], y2 = v2[1], z2 = v2[2];
  double x3 = v3[0], y3 = v3[1], z3 = v3[2];
  double Ix = (x1*(x1+x2+x3) + x2*(x2+x3) + x3*x3) * 0.1;
  double Iy = (y1*(y1+y2+y3) + y2*(y2+y3) + y3*y3) * 0.1;
  double Iz = (z1*(z1+z2+z3) + z2*(z2+z3) + z3*z3) * 0.1;

  myT[0] = Iy + Iz;
  myT[1] = Ix + Iz;
  myT[2] = Ix + Iy;
  myT[3] = (x1*y1+x2*y2+x3*y3)*0.1 + (x1*(y2+y3)+x2*(y1+y3)+x3*(y1+y2))*0.05;
  myT[4] = (x1*z1+x2*z2+x3*z3)*0.1 + (x1*(z2+z3)+x2*(z1+z3)+x3*(z1+z2))*0.05;
  myT[5] = (y1*z1+y2*z2+y3*z3)*0.1 + (y1*(z2+z3)+y2*(z1+z3)+y3*(z1+z2))*0.05;

  double vol = (v1*(v2^v3))/6.0;
  for (int i = 0; i < 6; i++)
    myT[i] *= vol;

  return *this;
}


/*!
  Translate an inertia tensor according to the parallel-axis theorem.
*/

FFaTensor3& FFaTensor3::translateInertia(const FaVec3& x, double mass)
{
  double mx = mass*x.sqrLength();
  myT[0] += mx - mass*x[0]*x[0];
  myT[1] += mx - mass*x[1]*x[1];
  myT[2] += mx - mass*x[2]*x[2];
  myT[3] -= mass*x[0]*x[1];
  myT[4] -= mass*x[0]*x[2];
  myT[5] -= mass*x[1]*x[2];
  return *this;
}


/*!
  Get the von Mises value of the tensor.
*/

double FFaTensor3::vonMises() const
{
  return FFaTensorTransforms::vonMises(myT[0],myT[1],myT[2],
                                       myT[3],myT[4],myT[5]);
}


/*!
  Get the max shear value of the tensor.
  If it can't be found, HUGE_VAL will be returned.
*/

double FFaTensor3::maxShear() const
{
  double max, middle, min;
  if (FFaTensorTransforms::principalValues(myT[0],myT[1],myT[2],
                                           myT[3],myT[4],myT[5],
                                           max,middle,min))
    return FFaTensorTransforms::maxShearValue(max,min);
  else
    return HUGE_VAL;
}


/*!
  Get the max shear of the tensor as a directed vector.
  If it can't be found, a 0-vector will be returned.
*/

void FFaTensor3::maxShear(FaVec3& v) const
{
  double values[3], max[3], middle[3], min[3];
  if (FFaTensorTransforms::principalDirs(myT,values,max,middle,min) == 0)
  {
    FFaTensorTransforms::maxShearDir(3,max,min,v.getPt());
    v *= FFaTensorTransforms::maxShearValue(values[0],values[2]);
  }
  else
    v.clear();
}


/*!
  Get the (absolute) max principal of the tensor.
  If it can't be found, HUGE_VAL will be returned.
*/

double FFaTensor3::maxPrinsipal(bool absMax) const
{
  double max, middle, min;
  if (FFaTensorTransforms::principalValues(myT[0],myT[1],myT[2],
                                           myT[3],myT[4],myT[5],
                                           max,middle,min))
    return absMax ? (fabs(max) > fabs(min) ? max : min) : max;
  else
    return HUGE_VAL;
}


/*!
  Get the middle principal of the tensor.
  If it can't be found, HUGE_VAL will be returned.
*/

double FFaTensor3::middlePrinsipal() const
{
  double max, middle, min;
  if (FFaTensorTransforms::principalValues(myT[0],myT[1],myT[2],
                                           myT[3],myT[4],myT[5],
                                           max,middle,min))
    return middle;
  else
    return HUGE_VAL;
}


/*!
  Get the min principal of the tensor.
  If it can't be found, HUGE_VAL will be returned.
*/

double FFaTensor3::minPrinsipal() const
{
  double max, middle, min;
  if (FFaTensorTransforms::principalValues(myT[0],myT[1],myT[2],
                                           myT[3],myT[4],myT[5],
                                           max,middle,min))
    return min;
  else
    return HUGE_VAL;
}


/*!
  Get the principal values of the tensor.
  If they can't be found, HUGE_VAL will be returned.
*/

void FFaTensor3::prinsipalValues(double& max, double& middle, double& min) const
{
  if (!FFaTensorTransforms::principalValues(myT[0],myT[1],myT[2],
                                            myT[3],myT[4],myT[5],
                                            max,middle,min))
    max = middle = min = HUGE_VAL;
}


/*!
  Get a valid rotation matrix corresponding to the principal axes of the tensor.
  The associated principal values are also found in the corresponding order.
  If the matrix can't be found, the identity will be returned,
  along with a vector of HUGE_VAL.
*/

void FFaTensor3::prinsipalValues(FaVec3& values, FaMat33& rotation) const
{
  if (FFaTensorTransforms::principalDirs(myT,
                                         values.getPt(),
                                         rotation[0].getPt(),
                                         rotation[1].getPt(),
                                         rotation[2].getPt()))
  {
    rotation.setIdentity();
    values[0] = values[1] = values[2] = HUGE_VAL;
  }
  else if ((rotation[0] ^ rotation[1]).isParallell(rotation[2]) != 1)
  {
    // Swap the 2nd and 3rd
    std::swap(rotation[1],rotation[2]);
    std::swap(values[1],values[2]);
  }
}


/*!
  Global operators.
*/

FFaTensor3 operator- (const FFaTensor3& t)
{
  return FFaTensor3(-t.myT[0],-t.myT[1],-t.myT[2],
                    -t.myT[3],-t.myT[4],-t.myT[5]);
}


FFaTensor3 operator+ (const FFaTensor3& a, const FFaTensor3& b)
{
  return FFaTensor3(a.myT[0] + b.myT[0],
                    a.myT[1] + b.myT[1],
                    a.myT[2] + b.myT[2],
                    a.myT[3] + b.myT[3],
                    a.myT[4] + b.myT[4],
                    a.myT[5] + b.myT[5]);
}


FFaTensor3 operator- (const FFaTensor3& a, const FFaTensor3& b)
{
  return FFaTensor3(a.myT[0] - b.myT[0],
                    a.myT[1] - b.myT[1],
                    a.myT[2] - b.myT[2],
                    a.myT[3] - b.myT[3],
                    a.myT[4] - b.myT[4],
                    a.myT[5] - b.myT[5]);
}


FFaTensor3 operator* (const FFaTensor3& a, double d)
{
  return FFaTensor3(a.myT[0] * d,
                    a.myT[1] * d,
                    a.myT[2] * d,
                    a.myT[3] * d,
                    a.myT[4] * d,
                    a.myT[5] * d);
}

FFaTensor3 operator* (double d, const FFaTensor3& a)
{
  return a*d;
}


FFaTensor3 operator/ (const FFaTensor3& a, double d)
{
  if (fabs(d) < 1.0e-16)
    return FFaTensor3(HUGE_VAL);

  return FFaTensor3(a.myT[0] / d, a.myT[1] / d, a.myT[2] / d,
                    a.myT[3] / d, a.myT[4] / d, a.myT[5] / d);
}


bool operator== (const FFaTensor3& a, const FFaTensor3& b)
{
  return ((a.myT[0] == b.myT[0]) &&
          (a.myT[1] == b.myT[1]) &&
          (a.myT[2] == b.myT[2]) &&
          (a.myT[3] == b.myT[3]) &&
          (a.myT[4] == b.myT[4]) &&
          (a.myT[5] == b.myT[5]));
}


bool operator!= (const FFaTensor3& a, const FFaTensor3& b)
{
  return !(a==b);
}


FFaTensor3 operator* (const FFaTensor3& a, const FaMat33& m)
{
  FFaTensor3 result;
  FFaTensorTransforms::rotate(a.myT, m[0].getPt(), m[1].getPt(), m[2].getPt(),
                              result.myT);
  return result;
}

FFaTensor3 operator* (const FFaTensor3& a, const FaMat34& m)
{
  FFaTensor3 result;
  FFaTensorTransforms::rotate(a.myT, m[0].getPt(), m[1].getPt(), m[2].getPt(),
                              result.myT);
  return result;
}


FFaTensor3 operator* (const FaMat33& m, const FFaTensor3& a)
{
  return a*m;
}

FFaTensor3 operator* (const FaMat34& m, const FFaTensor3& a)
{
  return a*m;
}


std::ostream& operator<< (std::ostream& s, const FFaTensor3& t)
{
  s << t.myT[0] << ' ' << t.myT[1] << ' ' << t.myT[2]<< ' '
    << t.myT[3] << ' ' << t.myT[4] << ' ' << t.myT[5];
  return s;
}

std::istream& operator>> (std::istream& s, FFaTensor3& t)
{
  FFaTensor3 tmpT;
  s >> tmpT[0] >> tmpT[1] >> tmpT[2] >> tmpT[3] >> tmpT[4] >> tmpT[5];
  if (s) t = tmpT;
  return s;
}
