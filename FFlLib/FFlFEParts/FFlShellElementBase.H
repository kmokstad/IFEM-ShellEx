// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_SHELL_ELEMENT_BASE_H
#define FFL_SHELL_ELEMENT_BASE_H

#include "FFlLib/FFlElementBase.H"


class FFlShellElementBase : public FFlElementBase
{
public:
  FFlShellElementBase(int ID) : FFlElementBase(ID) {}
  virtual ~FFlShellElementBase() {}

  double getThickness() const;

  virtual double getMassDensity() const;

  static double offPlaneTol; //!< Relative out-of-plane point search tolerance
};

#endif
