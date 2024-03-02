// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef _FFAVOLUME_H_
#define _FFAVOLUME_H_

class FaVec3;
class FFaTensor3;


/*!
  Class with static functions for calculation of volume, volume center
  and volume moment (inertia) for tetrahedron, wedge and hexahedron object.
*/

class FFaVolume
{
public:
  // Volume only calculation.
  static void tetVolume (const FaVec3& v1, const FaVec3& v2,
                         const FaVec3& v3, const FaVec3& v4, double& vol);
  static void wedVolume (const FaVec3& v1, const FaVec3& v2,
                         const FaVec3& v3, const FaVec3& v4,
                         const FaVec3& v5, const FaVec3& v6, double& vol);
  static void hexVolume (const FaVec3& v1, const FaVec3& v2,
                         const FaVec3& v3, const FaVec3& v4,
                         const FaVec3& v5, const FaVec3& v6,
                         const FaVec3& v7, const FaVec3& v8, double& vol);

  // Calculation of volume center. The input vertices are then shifted
  // such that they become relative to the computed volume center.
  // The return value is the volume of the object.
  static double tetCenter (FaVec3& v1, FaVec3& v2,
                           FaVec3& v3, FaVec3& v4, FaVec3& vc);
  static double wedCenter (FaVec3& v1, FaVec3& v2,
                           FaVec3& v3, FaVec3& v4,
                           FaVec3& v5, FaVec3& v6, FaVec3& vc);
  static double hexCenter (FaVec3& v1, FaVec3& v2,
                           FaVec3& v3, FaVec3& v4,
                           FaVec3& v5, FaVec3& v6,
                           FaVec3& v7, FaVec3& v8, FaVec3& vc);

  // Calculation of volume moments (inertias). The input vertices
  // are here assumed to be relative to the objects volume center.
  static void tetMoment (const FaVec3& v1, const FaVec3& v2,
                         const FaVec3& v3, const FaVec3& v4, FFaTensor3& vm);
  static void wedMoment (const FaVec3& v1, const FaVec3& v2,
                         const FaVec3& v3, const FaVec3& v4,
                         const FaVec3& v5, const FaVec3& v6, FFaTensor3& vm);
  static void hexMoment (const FaVec3& v1, const FaVec3& v2,
                         const FaVec3& v3, const FaVec3& v4,
                         const FaVec3& v5, const FaVec3& v6,
                         const FaVec3& v7, const FaVec3& v8, FFaTensor3& vm);

private:
  static double tetVolume (const FaVec3& v1, const FaVec3& v2,
                           const FaVec3& v3, const FaVec3& v4);

  static double pyrVolume (const FaVec3& v1, const FaVec3& v2,
                           const FaVec3& v3, const FaVec3& v4,
                           const FaVec3& v5);

  static double pyrCenter (const FaVec3& v1, const FaVec3& v2,
                           const FaVec3& v3, const FaVec3& v4,
                           const FaVec3& v5, FaVec3& vc);

  static FFaTensor3 pyrMoment (const FaVec3& v1, const FaVec3& v2,
                               const FaVec3& v3, const FaVec3& v4);
};

#endif
