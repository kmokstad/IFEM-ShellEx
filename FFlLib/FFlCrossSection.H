// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_CROSS_SECTION_H
#define FFL_CROSS_SECTION_H

#include <vector>
#include <string>
#include <iostream>


/*!
  \brief Struct with beam cross section parameters.
*/

struct FFlCrossSection
{
  double A   = 0.0; //!< Cross section area
  double Iyy = 0.0; //!< 2nd moment of inertia about local Y-axis
  double Izz = 0.0; //!< 2nd moment of inertia about local Z-axis
  double Izy = 0.0; //!< Coupling moment of inertia for non-symmetric sections
  double J   = 0.0; //!< Torsional stiffness parameter
  double NSM = 0.0; //!< Non-structural mass per length unit
  double K1  = 0.0; //!< Shear reduction factor in local Y-direction
  double K2  = 0.0; //!< Shear reduction factor in local Z-direction
  double S1  = 0.0; //!< Shear center position in local Y-direction
  double S2  = 0.0; //!< Shear center position in local Z-direction

  //! \brief Default constructor.
  FFlCrossSection(double K = 0.0) { K1 = K2 = K; }
  //! \brief Constructor calculating cross section parameters from given data.
  //! \param[in] Type Geometric cross section type
  //! \param[in] Dim Parameter values defining the cross section geometry
  FFlCrossSection(const std::string& Type, const std::vector<double>& Dim);
};


//! \brief Output stream operator.
std::ostream& operator<<(std::ostream& os, const FFlCrossSection& cs);

#endif
