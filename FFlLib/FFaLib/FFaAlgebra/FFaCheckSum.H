// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFA_CHECK_SUM_H
#define FFA_CHECK_SUM_H

#include <vector>
#include <string>
#include <iostream>

class FaVec3;


class FFaCheckSum
{
public:
  FFaCheckSum();
  FFaCheckSum(const FFaCheckSum&) = delete;
  ~FFaCheckSum() { delete[] crc32_table; }

  FFaCheckSum& operator=(const FFaCheckSum& cs);
  bool operator==(const FFaCheckSum& cs) const;

  void reset() { checksum = addval = 0; }
  unsigned int getCurrent() const { return checksum; }

  void add(const std::string& e) { for (char c : e) checksum += doCRC(c); }
  void add(char e)               { checksum += doCRC(e); }
  void add(unsigned char e)      { checksum += doCRC(e); }
  void add(short int e)          { checksum += doCRC(e); }
  void add(unsigned short int e) { checksum += doCRC(e); }
  void add(int e)                { checksum += doCRC(e); }
  void add(unsigned int e)       { checksum += doCRC(e); }
  void add(float e);
  void add(double e, int precision = 0);
  void add(const FaVec3& e, int precision = 0);

  template<class T> void add(const std::vector<T>& v)
  {
    for (const T& e : v) this->add(e);
  }

private:
  unsigned int doCRC(unsigned int data);

  friend std::ostream& operator<<(std::ostream& s, const FFaCheckSum& cs);
  friend std::istream& operator>>(std::istream& s, FFaCheckSum& cs);

  unsigned int* crc32_table;
  unsigned int  checksum, addval;
  bool          isBigEndian;
};

#if FFA_DEBUG > 2
template<> void FFaCheckSum::add(const std::vector<double>& v);
#endif
#endif
