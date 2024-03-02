// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFA_ALGEBRA_H
#define FFA_ALGEBRA_H

class FaVec3;
class FaMat33;


/*!
  \brief Class with static matrix-vector methods used by link file parsers, etc.
*/

class FFaAlgebra
{
public:
  static void eccTransform6 (double mat[6][6], const FaVec3& X);
  static bool congruenceTransform (double** mat, const FaMat33& T,
				   int N, int node = 0);
};

#endif
