// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFaVec3.H
  \brief Point vectors in 3D space.
*/

#ifndef FFA_VEC3_H
#define FFA_VEC3_H

#include <iostream>
#include <cmath>


/*!
  \brief Enums representing vector and matrix indices.
*/
enum FFaVec3IdxEnum
{
  VX, //!< The X-component
  VY, //!< The Y-component
  VZ, //!< The Z-component
  VW  //!< Translation part of a FaMat34 object
};


/*!
  \brief Class for point vectors in 3D space.
  \details This class is used to represent 3D point vectors in all type of
  geometry calculations. The class is equipped with a lot of convenience methods
  for the most common operations on such objects.
*/

class FaVec3
{
public:
  //! \brief Default constructor.
  FaVec3() { this->clear(); }
  //! \brief Constructor creating a vector from three double values.
  FaVec3(double x, double y, double z) { n[VX] = x; n[VY] = y; n[VZ] = z; }
  //! \brief Constructor creating a vector from a double array.
  FaVec3(const double* v) { n[VX] = v[VX];   n[VY] = v[VY];   n[VZ] = v[VZ]; }
  //! \brief Constructor creating a vector from a float array.
  FaVec3(const float* v)  { n[VX] = v[VX];   n[VY] = v[VY];   n[VZ] = v[VZ]; }
  //! \brief Copy constructor.
  FaVec3(const FaVec3& v) { n[VX] = v.n[VX]; n[VY] = v.n[VY]; n[VZ] = v.n[VZ]; }
  //! \brief Empty destructor.
  virtual ~FaVec3() {}

  //! \brief Assignment operator.
  FaVec3& operator=  (const FaVec3& v);
  //! \brief Addition operator.
  FaVec3& operator+= (const FaVec3& v);
  //! \brief Subtraction operator.
  FaVec3& operator-= (const FaVec3& v);
  //! \brief Scaling operator (multiplication by a scalar).
  FaVec3& operator*= (double d);
  //! \brief Scaling operator (division by a scalar).
  FaVec3& operator/= (double d);

  //! \brief Access to internal array (read only).
  const double* getPt() const { return n; }
  //! \brief Access to internal array (for manipulation).
  double*       getPt()       { return n; }

  //! \brief Zero-based indexing operator (read only).
  const double& operator[] (int i) const;
  //! \brief Zero-based indexing operator (for manipulation).
  double&       operator[] (int i);

  //! \brief One-based indexing operator (read only).
  const double& operator() (int i) const;
  //! \brief One-based indexing operator (for manipulation).
  double&       operator() (int i);

  //! \brief Returns the X-component of the vector.
  double x() const { return n[VX]; }
  //! \brief Returns the Y-component of the vector.
  double y() const { return n[VY]; }
  //! \brief Returns the Z-component of the vector.
  double z() const { return n[VZ]; }

  //! \brief Assigns the X-component of the vector.
  void x(double v) { n[VX] = v; }
  //! \brief Assigns the Y-component of the vector.
  void y(double v) { n[VY] = v; }
  //! \brief Assigns the Z-component of the vector.
  void z(double v) { n[VZ] = v; }
  //! \brief Sets the vector to zero.
  void clear() { n[VX] = n[VY] = n[VZ] = 0.0; }

  //! \brief Checks if this vector is parallel to \a otherVec.
  int isParallell (const FaVec3& otherVec, double tolerance = 1.0e-10) const;
  //! \brief Checks if this vector is equal to \a otherVec.
  bool equals     (const FaVec3& otherVec, double tolerance = 0.0) const;
  //! \brief Calculates angle between this vector and \a otherVec.
  double angle    (const FaVec3& otherVec) const;
  //! \brief Checks if this vector is zero.
  bool isZero     (double tolerance = 1.0e-10) const;

  //! \brief Returns the square of the length of this vector.
  double sqrLength() const { return n[VX]*n[VX] + n[VY]*n[VY] + n[VZ]*n[VZ]; }
  //! \brief Returns the length of this vector.
  double length() const { return sqrt(this->sqrLength()); }

  //! \brief Truncates small components of this vector to zero.
  FaVec3& truncate(double tolerance = 1.0e-10);
  //! \brief Normalizes this vector to unit length.
  FaVec3& normalize(double truncTol = 0.0);
  //! \brief Rounds off the components down to \a precision significant digits.
  FaVec3& round(int precision);

  //! \brief Sets this vector to the Cartesian equivalent
  //! of the cylindrical coordinates \a cylCoords.
  FaVec3& setByCylCoords(const FaVec3& cylCoords, FFaVec3IdxEnum axis = VZ);
  //! \brief Returns the cylindrical coordinates corresponding to
  //! the Cartesian position of this vector.
  FaVec3  getAsCylCoords(FFaVec3IdxEnum axis = VZ) const;

  //! \brief Sets this vector to the Cartesian equivalent
  //! of the spherical coordinates \a sphCoords.
  FaVec3& setBySphCoords(const FaVec3& sphCoords, FFaVec3IdxEnum axis = VZ);
  //! \brief Returns the spherical coordinates corresponding to
  //! the Cartesian position of this vector.
  FaVec3  getAsSphCoords(FFaVec3IdxEnum axis = VZ) const;

  //! \brief Unary negation operator.
  friend FaVec3   operator-  (const FaVec3& v);
  //! \brief Global addition operator.
  friend FaVec3   operator+  (const FaVec3& a, const FaVec3& b);
  //! \brief Global subtraction operator.
  friend FaVec3   operator-  (const FaVec3& a, const FaVec3& b);
  //! \brief Scaling by scalar.
  friend FaVec3   operator*  (const FaVec3& a,       double  d);
  //! \brief Scaling by scalar.
  friend FaVec3   operator*  (      double  d, const FaVec3& a);
  //! \brief Division by scalar.
  friend FaVec3   operator/  (const FaVec3& a,       double  d);
  //! \brief Dot product.
  friend double   operator*  (const FaVec3& a, const FaVec3& b);
  //! \brief Cross product.
  friend FaVec3   operator^  (const FaVec3& a, const FaVec3& b);
  //! \brief Equality operator.
  friend bool     operator== (const FaVec3& a, const FaVec3& b);
  //! \brief Unequality operator.
  friend bool     operator!= (const FaVec3& a, const FaVec3& b);

  //! \brief Printing operator.
  friend std::ostream& operator<< (std::ostream& s, const FaVec3& v);
  //! \brief Reading operator.
  friend std::istream& operator>> (std::istream& s, FaVec3& v);

private:
  double n[3]; //!< The actual vector components
};


// --- inline functions ---

inline FaVec3& FaVec3::operator= (const FaVec3& v)
{
  if (this != &v)
    for (int i = 0; i < 3; i++)
      n[i] = v.n[i];

  return *this;
}

inline const double& FaVec3::operator[] (int i) const
{
#ifdef FFA_INDEXCHECK
  if (i < 0 || i > 2)
    std::cerr <<"FaVec3::operator[]: index i="<< i <<" is out of range [0,2]"
              << std::endl;
#endif
  return n[i];
}

inline double& FaVec3::operator[] (int i)
{
#ifdef FFA_INDEXCHECK
  if (i < 0 || i > 2)
    std::cerr <<"FaVec3::operator[]: index i="<< i <<" is out of range [0,2]"
              << std::endl;
#endif
  return n[i];
}

inline const double& FaVec3::operator() (int i) const
{
#ifdef FFA_INDEXCHECK
  if (i < 1 || i > 3)
    std::cerr <<"FaVec3::operator(): index i="<< i <<" is out of range [1,3]"
              << std::endl;
#endif
  return n[i-1];
}

inline double& FaVec3::operator() (int i)
{
#ifdef FFA_INDEXCHECK
  if (i < 1 || i > 3)
    std::cerr <<"FaVec3::operator(): index i="<< i <<" is out of range [1,3]"
              << std::endl;
#endif
  return n[i-1];
}

#endif
