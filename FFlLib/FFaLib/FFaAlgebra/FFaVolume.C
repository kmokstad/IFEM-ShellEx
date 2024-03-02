// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFaLib/FFaAlgebra/FFaVolume.H"
#include "FFaLib/FFaAlgebra/FFaTensor3.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"


/*!
  Volume calculation. The volume of an object is computed as a sum of
  tetrahedron and pyramid contributions to account for possibly warped faces.
*/

void FFaVolume::tetVolume (const FaVec3& v1, const FaVec3& v2,
                           const FaVec3& v3, const FaVec3& v4, double& vol)
{
  // Single tetrahedron, cannot be warped
  vol = tetVolume(v1,v2,v3,v4);
}


void FFaVolume::wedVolume (const FaVec3& v1, const FaVec3& v2,
                           const FaVec3& v3, const FaVec3& v4,
                           const FaVec3& v5, const FaVec3& v6, double& vol)
{
  // Estimated volume center
  FaVec3 v0 = (v1 + v2 + v3 + v4 + v5 + v6) / 6.0;

  // Three pyramids and two tetrahedrons
  vol = tetVolume(v1,v2,v3,v0)
    +   tetVolume(v6,v5,v4,v0)
    +   pyrVolume(v1,v4,v5,v2,v0)
    +   pyrVolume(v2,v5,v6,v3,v0)
    +   pyrVolume(v3,v6,v4,v1,v0);
}


void FFaVolume::hexVolume (const FaVec3& v1, const FaVec3& v2,
                           const FaVec3& v3, const FaVec3& v4,
                           const FaVec3& v5, const FaVec3& v6,
                           const FaVec3& v7, const FaVec3& v8, double& vol)
{
  // Estimated volume center
  FaVec3 v0 = (v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8) * 0.125;

  // Six pyramids
  vol = pyrVolume(v1,v2,v3,v4,v0)
    +   pyrVolume(v8,v7,v6,v5,v0)
    +   pyrVolume(v1,v5,v6,v2,v0)
    +   pyrVolume(v2,v6,v7,v3,v0)
    +   pyrVolume(v3,v7,v8,v4,v0)
    +   pyrVolume(v4,v8,v5,v1,v0);
}


/*!
  Calculation of volume center. The input vertices are then shifted
  such that they become relative to the computed volume center.
  The volume of the object is retured.
*/

double FFaVolume::tetCenter (FaVec3& v1, FaVec3& v2,
                             FaVec3& v3, FaVec3& v4, FaVec3& vc)
{
  double volT = tetVolume(v1,v2,v3,v4);

  vc  = (v1 + v2 + v3 + v4) * 0.25;
  v1 -= vc;
  v2 -= vc;
  v3 -= vc;
  v4 -= vc;

  return volT;
}


double FFaVolume::wedCenter (FaVec3& v1, FaVec3& v2,
                             FaVec3& v3, FaVec3& v4,
                             FaVec3& v5, FaVec3& v6, FaVec3& vc)
{
  FaVec3 v0 = (v1 + v2 + v3 + v4 + v5 + v6) / 6.0;
  FaVec3 x1 = (v1 + v2 + v3 + v0) * 0.25;
  FaVec3 x2 = (v6 + v5 + v4 + v0) * 0.25;
  FaVec3 x3, x4, x5;

  double vol1 = tetVolume(v1,v2,v3,v0);
  double vol2 = tetVolume(v6,v5,v4,v0);
  double vol3 = pyrCenter(v1,v4,v5,v2,v0,x3);
  double vol4 = pyrCenter(v2,v5,v6,v3,v0,x4);
  double vol5 = pyrCenter(v3,v6,v4,v1,v0,x5);
  double volW = vol1 + vol2 + vol3 + vol4 + vol5;

  vc = (x1*vol1 + x2*vol2 + x3*vol3 + x4*vol4 + x5*vol5) / volW;

  v1 -= vc;
  v2 -= vc;
  v3 -= vc;
  v4 -= vc;
  v5 -= vc;
  v6 -= vc;

  return volW;
}


double FFaVolume::hexCenter (FaVec3& v1, FaVec3& v2,
                             FaVec3& v3, FaVec3& v4,
                             FaVec3& v5, FaVec3& v6,
                             FaVec3& v7, FaVec3& v8, FaVec3& vc)
{
  FaVec3 v0 = (v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8) * 0.125;
  FaVec3 x1, x2, x3, x4, x5, x6;

  double vol1 = pyrCenter(v1,v2,v3,v4,v0,x1);
  double vol2 = pyrCenter(v8,v7,v6,v5,v0,x2);
  double vol3 = pyrCenter(v1,v5,v6,v2,v0,x3);
  double vol4 = pyrCenter(v2,v6,v7,v3,v0,x4);
  double vol5 = pyrCenter(v3,v7,v8,v4,v0,x5);
  double vol6 = pyrCenter(v4,v8,v5,v1,v0,x6);
  double volH = vol1 + vol2 + vol3 + vol4 + vol5 + vol6;

  vc = (x1*vol1 + x2*vol2 + x3*vol3 + x4*vol4 + x5*vol5 + x6*vol6) / volH;

  v1 -= vc;
  v2 -= vc;
  v3 -= vc;
  v4 -= vc;
  v5 -= vc;
  v6 -= vc;
  v7 -= vc;
  v8 -= vc;

  return volH;
}


/*!
  Calculation of volume moments (inertias). The input vertices
  are here assumed to be relative to the objects volume center.
  The moments are computed as a sum of pyramid and tetrahedron contributions,
  where the first vertex is at the volume center and the four/three other
  vertices are on a face.
*/

void FFaVolume::tetMoment (const FaVec3& v1, const FaVec3& v2,
                           const FaVec3& v3, const FaVec3& v4, FFaTensor3& vm)
{
  // Four tetrahedrons
  vm = FFaTensor3(v1,v3,v2) + FFaTensor3(v1,v2,v4)
     + FFaTensor3(v2,v3,v4) + FFaTensor3(v1,v4,v3);
}


void FFaVolume::wedMoment (const FaVec3& v1, const FaVec3& v2,
                           const FaVec3& v3, const FaVec3& v4,
                           const FaVec3& v5, const FaVec3& v6, FFaTensor3& vm)
{
  // Three pyramids and two tetrahedrons
  vm = FFaTensor3(v1,v3,v2)
    +  FFaTensor3(v4,v5,v6)
    +  pyrMoment(v1,v2,v5,v4)
    +  pyrMoment(v2,v3,v6,v5)
    +  pyrMoment(v3,v1,v4,v6);
}


void FFaVolume::hexMoment (const FaVec3& v1, const FaVec3& v2,
                           const FaVec3& v3, const FaVec3& v4,
                           const FaVec3& v5, const FaVec3& v6,
                           const FaVec3& v7, const FaVec3& v8, FFaTensor3& vm)
{
  // Six pyramids
  vm = pyrMoment(v4,v3,v2,v1)
    +  pyrMoment(v5,v6,v7,v8)
    +  pyrMoment(v1,v2,v6,v5)
    +  pyrMoment(v2,v3,v7,v6)
    +  pyrMoment(v3,v4,v8,v7)
    +  pyrMoment(v4,v1,v5,v8);
}


/*!
  Private methods.
*/

double FFaVolume::tetVolume (const FaVec3& v1, const FaVec3& v2,
                             const FaVec3& v3, const FaVec3& v4)
{
  // Volume of a tetrahedron as a triple vector product
  return (((v2-v1) ^ (v3-v1)) * (v4-v1)) / 6.0;
}


double FFaVolume::pyrVolume (const FaVec3& v1, const FaVec3& v2,
                             const FaVec3& v3, const FaVec3& v4,
                             const FaVec3& v5)
{
  FaVec3 v0 = (v1 + v2 + v3 + v4) * 0.25;

  // Volume of a pyramid as a sum of four tetrahedron contributions
  return tetVolume(v0,v1,v2,v5) + tetVolume(v0,v2,v3,v5)
    +    tetVolume(v0,v3,v4,v5) + tetVolume(v0,v4,v1,v5);
}


double FFaVolume::pyrCenter (const FaVec3& v1, const FaVec3& v2,
                             const FaVec3& v3, const FaVec3& v4,
                             const FaVec3& v5, FaVec3& vc)
{
  FaVec3 v0 = (v1 + v2 + v3 + v4) * 0.25;

  // Volume center of a pyramid as a sum of four tetrahedron contributions
  double vol1 = tetVolume(v0,v1,v2,v5);
  double vol2 = tetVolume(v0,v2,v3,v5);
  double vol3 = tetVolume(v0,v3,v4,v5);
  double vol4 = tetVolume(v0,v4,v1,v5);
  double volP = vol1 + vol2 + vol3 + vol4;

  if (volP > 0.0)
    vc = (v1*(vol4+vol1) + v2*(vol1+vol2) + v3*(vol2+vol3) + v4*(vol3+vol4) +
          (v0+v5)*volP) / (4.0*volP);
  else
    vc = v0*0.8 + v5*0.2;

  return volP;
}


FFaTensor3 FFaVolume::pyrMoment (const FaVec3& v1, const FaVec3& v2,
                                 const FaVec3& v3, const FaVec3& v4)
{
  FaVec3 v0 = (v1 + v2 + v3 + v4) * 0.25;

  // Volume moment of a pyramid as a sum of four tetrahedron contributions
  return FFaTensor3(v0,v1,v2) + FFaTensor3(v0,v2,v3)
    +    FFaTensor3(v0,v3,v4) + FFaTensor3(v0,v4,v1);
}
