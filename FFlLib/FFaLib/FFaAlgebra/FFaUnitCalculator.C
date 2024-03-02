// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <cstdlib>

#include "FFaLib/FFaAlgebra/FFaUnitCalculator.H"
#include "FFaLib/FFaAlgebra/FFaMath.H"
#include "FFaLib/FFaString/FFaTokenizer.H"


bool FFaUnitCalculator::operator==(const FFaUnitCalculator& cal) const
{
  if (this == &cal) return true;

  if (name != cal.name) return false;
  if (origGroup != cal.origGroup) return false;
  if (convGroup != cal.convGroup) return false;
  if (myConvFactors.size() != cal.myConvFactors.size()) return false;

  SingleUnitMapIter it = myConvFactors.begin();
  SingleUnitMapIter jt = cal.myConvFactors.begin();
  for (; it != myConvFactors.end(); it++, jt++)
    if (it->first != jt->first)
      return false;
    else if (it->second.factor != jt->second.factor)
      return false;
    else if (it->second.origUnit != jt->second.origUnit)
      return false;
    else if (it->second.convUnit != jt->second.convUnit)
      return false;

  return true;
}


bool FFaUnitCalculator::convert(double& v, const std::string& aName,
                                int prec) const
{
  SingleUnitMapIter it = myConvFactors.find(aName);
  if (it == myConvFactors.end())
    return false;

  if (prec > 0) // round to <prec> significant digits
    v = round(v*it->second.factor,prec);
  else
    v *= it->second.factor;

  return true;
}


void FFaUnitCalculator::addConversion(const std::string& propName, double sf,
                                      const std::string& origUnit,
                                      const std::string& convUnit)
{
  myConvFactors[propName].factor   = sf;
  myConvFactors[propName].origUnit = origUnit;
  myConvFactors[propName].convUnit = convUnit;
}


std::ostream& operator<<(std::ostream& os, const FFaUnitCalculator& ucal)
{
  os <<  "<\""<< ucal.name
     <<"\",\""<< ucal.origGroup
     <<"\",\""<< ucal.convGroup
     <<"\"";

  FFaUnitCalculator::SingleUnitMapIter cit;
  for (cit = ucal.myConvFactors.begin(); cit != ucal.myConvFactors.end(); cit++)
    os <<",\n    <"
       << cit->first <<","
       << cit->second.factor <<",\""
       << cit->second.origUnit <<"\",\""
       << cit->second.convUnit <<"\">";
  os << ">";

  return os;
}


std::istream& operator>>(std::istream& is, FFaUnitCalculator& ucal)
{
  char c; // read until the first non-space char
  while (is.get(c) && isspace(c));

  // check if the first char is a <
  if (c == '<')
  {
    FFaTokenizer tokens(is,'<','>');
    if (tokens.size() < 3)
    {
       std::cerr <<"Error in unit conversion tokens - check token definition"
                 << std::endl;
       return is;
    }

    ucal.name = tokens[0];
    ucal.origGroup = tokens[1];
    ucal.convGroup = tokens[2];

    for (size_t i = 3; i < tokens.size(); i++)
    {
      // sub-tokens for units
      FFaTokenizer unitToken(tokens[i],'<','>');
      if (unitToken.size() != 4)
        std::cerr <<"Error in unit conversion tokens - check token definition:"
                  << "\n\t"<< tokens[i] << std::endl;
      else
      {
        FFaUnitCalculator::SingleUnit readUnit;
        readUnit.factor = atof(unitToken[1].c_str());
        readUnit.origUnit = unitToken[2];
        readUnit.convUnit = unitToken[3];
        ucal.myConvFactors[unitToken[0]] = readUnit;
      }
    }
  }

  return is;
}


//////////////////////////////////////////////////////////////////////
// Calculator provider

const FFaUnitCalculator*
FFaUnitCalculatorProvider::getCalculator(const std::string& calcName) const
{
  if (calcName.empty() || calcName == "none")
    return NULL;

  std::map<std::string,FFaUnitCalculator>::const_iterator it = myCalcs.find(calcName);
  if (it != myCalcs.end())
    return &(it->second);

  std::cerr <<"Non-existing unit conversion \""<< calcName <<"\""<< std::endl;
  return NULL;
}


std::vector<const FFaUnitCalculator*>
FFaUnitCalculatorProvider::getCalculators() const
{
  std::vector<const FFaUnitCalculator*> tmpVec;
  tmpVec.reserve(myCalcs.size());
  for (const std::pair<const std::string,FFaUnitCalculator>& cit : myCalcs)
    tmpVec.push_back(&cit.second);

  return tmpVec;
}


void FFaUnitCalculatorProvider::getCalculatorNames(std::vector<std::string>& definedCalcs) const
{
  definedCalcs.reserve(myCalcs.size());
  for (const std::pair<const std::string,FFaUnitCalculator>& cit : myCalcs)
    definedCalcs.push_back(cit.first);
}


void FFaUnitCalculatorProvider::addCalculator(const FFaUnitCalculator& calc)
{
  if (calc.isValid())
    myCalcs[calc.getName()] = calc;
  else
    std::cerr <<"  ** FFaUnitCalculatorProvider::addCalculator:"
              <<" Invalid unit calculator (no name)."<< std::endl;
}


bool FFaUnitCalculatorProvider::readCalculatorDefs(const std::string& filename)
{
  std::ifstream is(filename.c_str());
  if (!is)
  {
    std::cerr <<"Can't open unit-conversion file "<< filename << std::endl;
    return false;
  }

  // search for begin-char:
  char c;
  while (is.get(c) && !is.eof())
  {
    while (!is.eof() && isspace(c))
      is.get(c);

    if (c == '#')
      is.ignore(BUFSIZ,'\n');
    else if (c == '<')
    {
      // found valid start of entry
      is.putback(c);
      FFaUnitCalculator aCal;
      is >> aCal;
      this->addCalculator(aCal);
    }
    else if (!is.eof())
    {
      std::cerr <<"Error in calculator definition file"<< std::endl;
      is.ignore(BUFSIZ,'\n');
    }
  }

  return true;
}


bool FFaUnitCalculatorProvider::printCalculatorDefs(const std::string& filename) const
{
  std::ofstream os(filename.c_str());
  if (!os)
  {
    std::cerr <<"Can't open unit conversion output file "<< filename
              << std::endl;
    return false;
  }

  os <<"#FEDEM converter defs";
  for (const std::pair<const std::string,FFaUnitCalculator>& cit : myCalcs)
    os <<'\n'<< cit.second;
  os << std::endl;

  return true;
}
