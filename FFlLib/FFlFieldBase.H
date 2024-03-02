// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_FIELD_BASE_H
#define FFL_FIELD_BASE_H

#include <vector>
#include <string>
#include <iostream>

class FFaCheckSum;


class FFlFieldBase
{
public:
  FFlFieldBase(){}
  virtual ~FFlFieldBase(){}

  virtual void calculateChecksum(FFaCheckSum* cs) const = 0;
  virtual bool parse(std::vector<std::string>::const_iterator& begin,
		     const std::vector<std::string>::const_iterator& end) = 0;

  static bool parseNumericField(int&, const std::string& field, bool silent = false);
  static bool parseNumericField(double&, const std::string& field);

protected:
  virtual void write(std::ostream& os) const = 0;

  friend std::ostream& operator<<(std::ostream& os, const FFlFieldBase& field);
};

#endif
