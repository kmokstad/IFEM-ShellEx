// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFaTensorTransforms.H
  \brief Tensor transformation utilities.
*/

#ifndef FFA_TENSOR_TRANSFORMS_H
#define FFA_TENSOR_TRANSFORMS_H


namespace FFaTensorTransforms //! Various tensor transformation utilities
{
  //! \brief Returns the von Mises value of a symmetric stress/strain tensor.
  double vonMises(int N, const double* S);
  //! \brief Returns the von Mises value of a 2D stress/strain tensor.
  double vonMises(double s11, double s22, double s12);
  //! \brief Returns the von Mises value of a 3D stress/strain tensor.
  double vonMises(double s11, double s22, double s33,
                  double s12, double s13, double s23);

  //! \brief Calculates the principal directions and values of a tensor.
  int principalDirs(int N, const double* S, double* pVals, double* pDirs);
  //! \brief Calculates the principal directions and values of a 2D tensor.
  int principalDirs(const double* S,
                    double* pVal, double* p1Dir, double* p2Dir);
  //! \brief Calculates the principal directions and values of a 3D tensor.
  int principalDirs(const double* S,
                    double* pVal, double* p1Dir, double* p2Dir, double* p3Dir);

  //! \brief Calculates the principal values of a tensor.
  bool principalValues(int N, const double* S, double* P);
  //! \brief Calculates the principal values of a 2D tensor.
  bool principalVals2D(double s11, double s22, double s12, double* pVal);
  //! \brief Calculates the principal values of a 3D tensor.
  bool principalVals3D(double s11, double s22, double s33,
                       double s12, double s13, double s23, double* pVal);

  // For backward compatibility, purge these two later..
  inline bool principalValues(double s11, double s22, double s12,
                              double& p1, double& p2)
  {
    double p[3] = { 0.0, 0.0, 0.0 };
    bool ok = principalVals2D(s11,s22,s12,p);
    p1 = p[0];
    p2 = p[1];
    return ok;
  }
  inline bool principalValues(double s11, double s22, double s33,
                              double s12, double s13, double s23,
                              double& p1, double& p2, double& p3)
  {
    double p[3] = { 0.0, 0.0, 0.0 };
    bool ok = principalVals3D(s11,s22,s33,s12,s13,s23,p);
    p1 = p[0];
    p2 = p[1];
    p3 = p[2];
    return ok;
  }

  //! \brief Returns the maximum shear value.
  double maxShearValue(double pMax, double pMin);
  //! \brief Returns the direction of maximum shear.
  void   maxShearDir  (int N, const double* pMaxDir, const double* pMinDir,
                       double* maxShearDir);

  //! \brief Coordinate transformation for a 2D tensor.
  void rotate2D(const double* inTensor,
                const double* rotMx,
                double* outTensor);
  //! \brief Coordinate transformation for a 2D tensor.
  void rotate  (const double* inTensor,
                const double* eX, const double* eY,
                double* outTensor);
  //! \brief Coordinate transformation for a 3D tensor.
  void rotate3D(const double* inTensor,
                const double* rotMx,
                double* outTensor);
  //! \brief Coordinate transformation for a 3D tensor.
  void rotate  (const double* inTensor,
                const double* eX, const double* eY, const double* eZ,
                double* outTensor);

  //! \brief Plain 2D-to-3D tensor conversion.
  void from2Dto3D(const double* tensor2D, double* tensor3D);
  //! \brief Plain 3D-to-2D tensor conversion.
  void from3Dto2D(const double* tensor3D, double* tensor2D);
}

#endif
