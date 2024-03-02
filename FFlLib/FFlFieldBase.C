// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlFieldBase.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"
#include <cstdlib>
#include <cstring>
#include <climits>


bool FFlFieldBase::parseNumericField(int& val, const std::string& field, bool silent)
{
  if (field.empty())
    return true;

  // Convert the string to an int
  const char* startPtr = field.c_str();
  char* endPtr = NULL;
  long int lval = strtol(startPtr,&endPtr,10);
  if (endPtr == startPtr+strlen(startPtr))
  {
    if (lval <= INT_MAX && -lval <= INT_MAX)
    {
      val = static_cast<int>(lval);
      return true;
    }
    else if (!silent)
    {
      ListUI <<"\n *** Error: Integer field overflow \""<< startPtr <<"\"\n";
      return false;
    }
  }

  if (silent)
    return false; // skip the error message

  // A decimal point is allowed only if it is at the end or followed by 0's only
  bool integer = true;
  bool haveDot = false;
  for (char* p = endPtr; p < startPtr+strlen(startPtr); p++)
    if (!haveDot && p[0] == '.')
      haveDot = true;
    else if (!(haveDot && p[0] == '0'))
      integer = false;
  if (integer) return true;

  ListUI <<"\n *** Error: Can not convert field to int: \""<< startPtr <<"\"\n"
         <<"                                           ";
  for (const char* q = startPtr; q < endPtr; q++) ListUI <<" ";
  ListUI <<"^\n";
  return false;
}


bool FFlFieldBase::parseNumericField(double& val, const std::string& field)
{
  if (field.empty())
    return true;

  // Preprocess the field adding the exponential 'E' if missing
  char c1, c2;
  for (int i = field.length()-1; i > 0; i--)
  {
    c1 = field[i];
    if (c1 == '-' || c1 == '+')
    {
      c2 = field[i-1];
      if (c2 != 'e' && c2 != 'E')
        const_cast<std::string&>(field).insert(i,1,'E');
      break;
    }
  }

  // Convert the string to a double
  const char* startPtr = field.c_str();
  char* endPtr = NULL;
  val = strtod(startPtr,&endPtr);
  if (endPtr == startPtr+strlen(startPtr))
    return true;

  ListUI <<"\n *** Error: Can not convert field to double: \""<< field <<"\"\n"
         <<"                                              ";
  for (const char* q = startPtr; q < endPtr; q++) ListUI <<" ";
  ListUI <<"^\n";
  return false;
}


std::ostream& operator<<(std::ostream& os, const FFlFieldBase& field)
{
  field.write(os);
  return os;
}
