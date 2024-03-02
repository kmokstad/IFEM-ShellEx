// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFaLib/FFaAlgebra/FFaCheckSum.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"

#define CRC32_POLY 0x04c11db7


FFaCheckSum::FFaCheckSum()
{
  this->reset();

  crc32_table = new unsigned int[256];
  for (int i = 0; i < 256; i++)
  {
    unsigned int c = i<<24;
    for (int j = 0; j < 8; j++)
      c = c & 0x80000000 ? (c<<1)^CRC32_POLY : (c<<1);
    crc32_table[i] = c;
  }

  // determine endian format
  int nl = 0x12345678;
  unsigned char* p = (unsigned char*)&nl; // 32-bit integer
  if (p[0] == 0x12 && p[1] == 0x34 && p[2] == 0x56 && p[3] == 0x78)
    isBigEndian = true;
  else //if (p[0] == 0x78 && p[1] == 0x56 && p[2] == 0x34 && p[3] == 0x12)
    isBigEndian = false;

#if FFA_DEBUG > 1
  std::cout <<"FFaCheckSum: system is "
            << (isBigEndian? "Big" : "Little")
            <<"Endian"<< std::endl;
#endif
}


FFaCheckSum& FFaCheckSum::operator=(const FFaCheckSum& cs)
{
  if (this == &cs)
    return *this;

  checksum = cs.checksum;
  addval   = cs.addval;

  return *this;
}


bool FFaCheckSum::operator==(const FFaCheckSum& cs) const
{
  if (this == &cs)
    return true;

  return checksum == cs.checksum;
}


unsigned int FFaCheckSum::doCRC(unsigned int data)
{
  data += addval++;

  unsigned int crc = 0xffffffff;
  crc=(crc<<8)^crc32_table[(crc>>24)^((unsigned char)(data>>24))];
  crc=(crc<<8)^crc32_table[(crc>>24)^((unsigned char)(data>>16))];
  crc=(crc<<8)^crc32_table[(crc>>24)^((unsigned char)(data>>8))];
  crc=(crc<<8)^crc32_table[(crc>>24)^((unsigned char)(data))];

  return crc ? crc : 0xdeadbeef;
}


void FFaCheckSum::add(const FaVec3& e, int precision)
{
#if FFA_DEBUG > 2
  int savePrecision = std::cout.precision(16);
#endif

  // cast to float to avoid round-of problems since FaVec3 objects may contain
  // results of computations that are written to ftl-files with lower precision
  for (int i = 0; i < 3; i++)
  {
#if FFA_DEBUG > 2
    std::cout <<"Checksum after FaVec3 element #"<< i <<", "<< e[i];
#endif
    if (precision > 0) // round to <precision> significant digits
      this->add(e[i],precision);
    else // cast to float
      this->add((float)e[i]);
#if FFA_DEBUG > 2
    std::cout <<" : "<< checksum << std::endl;
#endif
  }

#if FFA_DEBUG > 2
  std::cout.precision(savePrecision);
#endif
}


#if FFA_DEBUG > 2
template<> void FFaCheckSum::add(const std::vector<double>& v)
{
  int idx = 0, savePrecision = std::cout.precision(16);
  for (double e : v)
  {
    std::cout <<"Checksum after vector<double> element #"<< idx++ <<", "<< e;
    this->add(e);
    std::cout <<" : "<< checksum << std::endl;
  }
  std::cout.precision(savePrecision);
}
#endif


void FFaCheckSum::add(float e)
{
  unsigned char* p = (unsigned char*)&e;

  // ommit the 4th byte containing the least significant digits

  int i;
  if (isBigEndian)
  {
    // check if only the sign bit is set, if so reset it to zero
    bool onlySignBit = (p[0] == 128);
    for (i = 1; onlySignBit && i < 3; i++)
      if (p[i] != 0) onlySignBit = false;
    if (onlySignBit) p[0] = 0;

    for (i = 0; i < 3; i++)
      checksum += doCRC(p[i]);
  }
  else // Little Endian, swap bytes
  {
    // check if only the sign bit is set, if so reset it to zero
    bool onlySignBit = (p[3] == 128);
    for (i = 2; onlySignBit && i > 0; i--)
      if (p[i] != 0) onlySignBit = false;
    if (onlySignBit) p[3] = 0;

    for (i = 3; i > 0; i--)
      checksum += doCRC(p[i]);
  }
}


void FFaCheckSum::add(double e, int precision)
{
  if (precision > 0 && e != 0.0)
  {
    // Round the value <e> to <precision> significant digits
    int expon = 1-precision;
    double ae = e > 0.0 ? e : -e;
    for (; ae >= 10.0; expon++) ae /= 10.0;
    for (; ae <   1.0; expon--) ae *= 10.0;
    long long int le = round(e/pow(10.0,expon));
#if FFA_DEBUG > 3
    std::cout <<" = "<< le <<"E"<< expon;
#endif
    // Calculate the checksum from the mantissa and exponent
    checksum += doCRC(le) + doCRC(expon);
    return;
  }

  unsigned char* p = (unsigned char*)&e;

  // ommit the 7th and 8th byte containing the least significant digits

  int i;
  if (isBigEndian)
  {
    // check if only the sign bit is set, if so reset it to zero
    bool onlySignBit = (p[0] == 128);
    for (i = 1; onlySignBit && i < 6; i++)
      if (p[i] != 0) onlySignBit = false;
    if (onlySignBit) p[0] = 0;

    for (i = 0; i < 6; i++)
      checksum += doCRC(p[i]);
  }
  else // LittleEndian, swap bytes
  {
    // check if only the sign bit is set, if so reset it to zero
    bool onlySignBit = (p[7] == 128);
    for (i = 6; onlySignBit && i > 1; i--)
      if (p[i] > 0) onlySignBit = false;
    if (onlySignBit) p[7] = 0;

    for (i = 7; i > 1; i--)
      checksum += doCRC(p[i]);
  }
}


std::ostream& operator<<(std::ostream& s, const FFaCheckSum& cs)
{
  return s << cs.checksum <<" "<< cs.addval;
}


std::istream& operator>>(std::istream& s, FFaCheckSum& cs)
{
  s >> cs.checksum >> cs.addval;
  return s;
}
