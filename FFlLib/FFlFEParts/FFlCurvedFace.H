// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_CURVED_FACE_H
#define FFL_CURVED_FACE_H

#include <vector>

class FaVec3;

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif

class FFlNode;


class FFlCurvedFace
{
public:
  static bool faceNormals(const std::vector<FFlNode*>& nodes,
                          std::vector<FaVec3>& normals);

private:
  static void shapeDerivs6(double* N1, double* N2, double X1, double X2);
  static void shapeDerivs8(double* N1, double* N2, double X1, double X2);
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
