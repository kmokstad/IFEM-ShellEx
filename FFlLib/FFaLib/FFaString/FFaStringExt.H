// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFaStringExt.H
  \brief Extensions of the STL string class.
*/

#ifndef FFA_STRING_EXT_H
#define FFA_STRING_EXT_H

#include <string>
#include <cstdio>
#include <cctype>


/*!
  \brief Extension of the STL string class converting numbers to a text string.
*/

class FFaNumStr : public std::string
{
public:
  //! \brief Constructor converting an integer value with a format string.
  FFaNumStr(const char* format, int val)
  {
    char buf[128];
    snprintf(buf,128,format,val);
    append(buf);
  }

  //! \brief Constructor converting an integer value with a format string.
  FFaNumStr(const char* format, size_t val)
  {
    char buf[128];
    snprintf(buf,128,format,val);
    append(buf);
  }

  //! \brief Constructor converting a float value with a format string.
  FFaNumStr(const char* format, float val)
  {
    char buf[128];
    snprintf(buf,128,format,val);
    append(buf);
  }

  //! \brief Constructor converting a double value with a format string.
  FFaNumStr(const char* format, double val)
  {
    char buf[128];
    snprintf(buf,128,format,val);
    append(buf);
  }

  //! \brief Constructor converting an integer value.
  FFaNumStr(int val)          { this->append(std::to_string(val)); }
  //! \brief Constructor converting an unsigned integer value.
  FFaNumStr(unsigned int val) { this->append(std::to_string(val)); }
  //! \brief Constructor converting a long integer value.
  FFaNumStr(long int val)     { this->append(std::to_string(val)); }

  // These two constructors are documented in the C-file.
  FFaNumStr(double val, char f = 'g', int precision = 6);
  FFaNumStr(double val, int integerDigits, int precision = 12,
	    double ceiling = 1.0e+7, double floor = 1.0e-5,
	    bool useDigitGrouping = false);
};


/*!
  \brief Extension of the STL string class for lower-case strings only.
*/

class FFaLowerCaseString : public std::string
{
public:
  //! \brief The copy constructor converts the string to lower case.
  FFaLowerCaseString(const std::string& s)
  {
    this->reserve(s.length());
    for (size_t i = 0; i < s.length(); i++)
      this->operator+=(tolower(s[i]));
  }
};


/*!
  \brief Extension of the STL string class for upper-case strings only.
*/

class FFaUpperCaseString : public std::string
{
public:
  //! \brief The copy constructor converts the string to upper case.
  FFaUpperCaseString(const std::string& s)
  {
    this->reserve(s.length());
    for (size_t i = 0; i < s.length(); i++)
      this->operator+=(toupper(s[i]));
  }
};


/*!
  \brief Extension of the STL string class with some parsing methods.
  \details This is mainly used for description-field command parsing.
*/

class FFaString : public std::string
{
public:
  //! \brief Copy constructor.
  FFaString(const std::string& s) : std::string(s) {}
  //! \brief Checks if the string contains the specified substring.
  bool hasSubString(const char*) const;
  //! \brief Parses an integer after the specified substring.
  int getIntAfter(const char*) const;
  //! \brief Parses a double after the specified substring.
  double getDoubleAfter(const char*) const;
  //! \brief Parses a list of integers after the specified substring.
  int  getIntsAfter(const char*, const int, int*) const;
  //! \brief Parses a list of doubles after the specified substring.
  int getDoublesAfter(const char*, const int, double*) const;
  //! \brief Parses a text string after the specified substring.
  std::string getTextAfter(const char*, const char* end = NULL) const;

private:
  //! \brief Returns the position of the specified substring.
  size_t getPosAfterString(const char*) const;
};

#endif
