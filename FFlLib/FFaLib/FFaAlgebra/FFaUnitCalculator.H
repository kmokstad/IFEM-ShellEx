#ifndef FFA_UNIT_CALCULATOR_H
#define FFA_UNIT_CALCULATOR_H

#include <string>


//! \brief Dummy class doing nothing.
class FFaUnitCalculator
{
public:
  template<class T> bool convert(T&, const std::string&) const { return true; }
  bool convert(double&, const std::string&, int) const { return true; }
};

#endif
