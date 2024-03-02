// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFA_UNIT_CALCULATOR_H
#define FFA_UNIT_CALCULATOR_H

#include <string>
#include <vector>
#include <map>

#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFaUnitCalculator
{
  struct SingleUnit
  {
    double factor;
    std::string origUnit;
    std::string convUnit;
  };

  typedef std::map<std::string,SingleUnit> SingleUnitMap;
  typedef SingleUnitMap::const_iterator    SingleUnitMapIter;

public:
  FFaUnitCalculator(const char* n = NULL) : origGroup("SI"), convGroup("SI")
  {
    if (n) name = n;
  }

  bool operator==(const FFaUnitCalculator& cal) const;

  bool isValid() const { return !name.empty(); }

  const std::string& getName()      const { return name; }
  const std::string& getOrigUnits() const { return origGroup; }
  const std::string& getConvUnits() const { return convGroup; }

  template <class T> bool convert(T& value, const std::string& aName) const
  {
    SingleUnitMapIter it = myConvFactors.find(aName);
    if (it == myConvFactors.end()) return false;

    value *= it->second.factor;
    return true;
  }

  //! \brief Conversion with rounding to \a prec significant digits.
  bool convert(double& value, const std::string& aName, int prec) const;

  //! \brief Reversed conversion.
  template <class T> bool inverse(T& value, const std::string& aName) const
  {
    SingleUnitMapIter it = myConvFactors.find(aName);
    if (it == myConvFactors.end()) return false;

    value /= it->second.factor;
    return true;
  }

  double convFactor(const std::string& aName) const
  {
    SingleUnitMapIter it = myConvFactors.find(aName);
    return it == myConvFactors.end() ?  0.0 : it->second.factor;
  }

  std::string origUnitName(const std::string& aName) const
  {
    SingleUnitMapIter it = myConvFactors.find(aName);
    return it == myConvFactors.end() ? std::string("") : it->second.origUnit;
  }

  std::string convUnitName(const std::string& aName) const
  {
    SingleUnitMapIter it = myConvFactors.find(aName);
    return it == myConvFactors.end() ? std::string("") : it->second.convUnit;
  }

  void addConversion(const std::string& propName, double sf,
		     const std::string& origUnit = "",
		     const std::string& convUnit = "");

private:
  std::string name;
  std::string origGroup;
  std::string convGroup;
  SingleUnitMap myConvFactors;

  friend std::ostream& operator<< (std::ostream& os, const FFaUnitCalculator& field);
  friend std::istream& operator>> (std::istream& is, FFaUnitCalculator& field);
};


class FFaUnitCalculatorProvider : public FFaSingelton<FFaUnitCalculatorProvider>
{
public:
  const FFaUnitCalculator* getCalculator(const std::string& calcName) const;

  void getCalculatorNames(std::vector<std::string>& definedCalcs) const;
  std::vector<const FFaUnitCalculator*> getCalculators() const;

  void addCalculator(const FFaUnitCalculator& calc);

  //! \brief Reads/writes calculator data from file.
  bool readCalculatorDefs(const std::string& filename);
  //! \brief Writes calculator data to file.
  bool printCalculatorDefs(const std::string& filename) const;

protected:
  FFaUnitCalculatorProvider() {}
  virtual ~FFaUnitCalculatorProvider() {}

private:
  std::map<std::string,FFaUnitCalculator> myCalcs;

  friend FFaUnitCalculatorProvider* FFaSingelton<FFaUnitCalculatorProvider>::instance();
};

#endif
