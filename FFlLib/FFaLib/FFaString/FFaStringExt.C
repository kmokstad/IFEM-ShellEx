// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFaStringExt.C
  \brief Extensions of the STL string class.
*/

#include <cstdlib>
#include <cstring>

#include "FFaLib/FFaString/FFaStringExt.H"


/*!
  \brief Constructs a string from a double value with a given precision.
  \details This is the default double-to-string method.
  \param[in] val The value to construct the string from
  \param[in] f Format specifier (typically, 'f' or 'g')
  \param[in] precision Max number of digits after decimal point
*/

FFaNumStr::FFaNumStr(double val, char f, int precision)
{
  char buf[64];
  char format[8];
  strncpy(format,"% .*lf",8); format[5] = f;
  snprintf(buf, 64, format, precision, val);
  this->append(buf);
}


/*!
  \brief Constructs a string from a double value with a given precision.
  \param[in] val The value to construct the string from
  \param[in] integerDigits Number of digits after the decimal point
  when \a val is an integer value (0, 1, 2...)
  \param[in] precision Max number of digits after decimal point
  \param[in] ceiling If the absolute value of \a val is greater than or equal to
  this value, scientific notation is used (e.g., 1.98e+09)
  \param[in] floor If the absolute value of \a val is less than or equal to
  this value, scientific notation (e.g., 1.98e-09).
  If zero, decimal notation is enforced.
  \param[in] useDigitGrouping If \e true, any digits before the decimal point
  are grouped in space-separated triplets
*/

FFaNumStr::FFaNumStr(double val, int integerDigits, int precision,
		     double ceiling, double floor, bool useDigitGrouping)
{
  char buf[128];

  // Maximum precision limit
  if (val > -1.0e-15 && val < 1.0e-15) // truncate to 0.0
    val = 0.0;

  // Print the string
  if (val >= ceiling || -val >= ceiling)
    snprintf(buf, 128, "% .*g", precision, val);
  else if ((double)((int)val) == val)
  {
    if (integerDigits < 0)
      snprintf(buf, 128, "% d", (int)val);
    else
      snprintf(buf, 128, "% #.*f", integerDigits, val);
  }
  else if (val <= floor && -val <= floor)
    snprintf(buf, 128, "% .*g", precision, val);
  else if (floor <= 0.0)
    snprintf(buf, 128, "% #.*f", precision, val);
  else
    snprintf(buf, 128, "% .*g", precision, val);

  this->append(buf);

  // Always decimal dot on real numbers
  size_t posDot = this->find('.');
  if (posDot == npos && integerDigits >= 0)
  {
    size_t posExp = this->find('e');
    if (posExp == npos) {
      posDot = this->size();
      this->append(".0"); // append
    }
    else {
      posDot = posExp;
      this->insert(posExp,".0"); // insert before e
    }
  }

  // Digit grouping?
  if (useDigitGrouping && posDot != npos)
    for (int i = posDot-3; i > 0; i -= 3)
      this->insert(i,1,' ');
}


/*
  Special functions used by the solver parser.
  Mainly for parsing description-field commands (beta features).
*/

size_t FFaString::getPosAfterString(const char* s) const
{
  if (!s) return npos;

  size_t ipos = this->find(s);
  if (ipos == npos)
  {
    // Check for case-insensitive match
    FFaUpperCaseString THIS(this->c_str());
    ipos = THIS.find(FFaUpperCaseString(s));
    if (ipos == npos) return npos;
  }

  return ipos + strlen(s);
}


bool FFaString::hasSubString(const char* s) const
{
  return this->getPosAfterString(s) == npos ? false : true;
}


int FFaString::getIntAfter(const char* s) const
{
  size_t ipos = this->getPosAfterString(s);
  if (ipos >= this->size()-1) return 0;

  return atoi(this->c_str()+ipos);
}


double FFaString::getDoubleAfter(const char* s) const
{
  size_t ipos = this->getPosAfterString(s);
  if (ipos >= this->size()-1) return 0.0;

  return atof(this->c_str()+ipos);
}


int FFaString::getIntsAfter(const char* s, const int n, int* v) const
{
  int i;
  for (i = 0; i < n; i++) v[i] = 0;
  if (!s || n < 1) return 0;

  size_t ipos = this->getPosAfterString(s);
  if (ipos >= this->size()-1) return 0;

  for (i = 0; i < n; i++)
  {
    ipos = this->find_first_not_of(' ',ipos);
    if (ipos == npos) return i;
    char c = this->operator[](ipos);
    if (c != '-' && !isdigit(c)) return i;
    v[i] = atoi(this->c_str()+ipos);
    ipos = this->find_first_of(' ',ipos);
    if (ipos == npos) return i+1;
  }

  return n;
}


int FFaString::getDoublesAfter(const char* s, const int n, double* v) const
{
  int i;
  for (i = 0; i < n; i++) v[i] = 0.0;
  if (!s || n < 1) return 0;

  size_t ipos = this->getPosAfterString(s);
  if (ipos >= this->size()-1) return 0;

  for (i = 0; i < n; i++)
  {
    ipos = this->find_first_not_of(' ',ipos);
    if (ipos == npos) return i;
    char c = this->operator[](ipos);
    if (c != '-' && c != '.' && !isdigit(c)) return i;
    v[i] = atof(this->c_str()+ipos);
    ipos = this->find_first_of(' ',ipos);
    if (ipos == npos) return i+1;
  }

  return n;
}


std::string FFaString::getTextAfter(const char* s, const char* end) const
{
  size_t ipos = this->getPosAfterString(s);
  while (ipos < this->size() && isspace(this->operator[](ipos)))
    ipos++;

  if (ipos >= this->size()-1)
    return "";
  else if (end)
    return this->substr(ipos,this->find(end,ipos)-ipos);
  else
    return this->substr(ipos);
}
