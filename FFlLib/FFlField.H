// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_FIELD_H
#define FFL_FIELD_H

#include "FFlFieldBase.H"
#include "FFaLib/FFaAlgebra/FFaCheckSum.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"


template <class value>
class FFlField : public FFlFieldBase
{
public:
  FFlField() {}
  FFlField(const FFlField&) = delete;
  virtual ~FFlField() {}

  value& operator=(const value& data)
  {
    myData = data;
    return myData;
  }

  FFlField& operator=(const FFlField&) = default;

  const value& getValue() const { return myData; }

  value& data() { return myData; }

  void setValue(const value& newData) { myData = newData; }

  virtual void calculateChecksum(FFaCheckSum* cs) const { cs->add(myData); }

  virtual bool parse(std::vector<std::string>::const_iterator& begin,
                     const std::vector<std::string>::const_iterator& end);

protected:
  virtual void write(std::ostream& os) const { os << myData; }

  bool isEqual(const FFlField& ref) const { return myData == ref.myData; }

private:
  value myData;

  friend bool operator==(const FFlField& f1, const FFlField& f2) { return f1.isEqual(f2); }
  friend bool operator!=(const FFlField& f1, const FFlField& f2) { return !f1.isEqual(f2); }
};


template<> inline FFlField<int>::FFlField()            { myData = 0; }
template<> inline FFlField<unsigned short>::FFlField() { myData = 0; }
template<> inline FFlField<double>::FFlField()         { myData = 0.0; }


template<> inline void
FFlField<std::string>::write(std::ostream& os) const
{
  os << "\"" << myData << "\"";
}


template<> inline bool
FFlField<double>::parse(std::vector<std::string>::const_iterator& begin,
			const std::vector<std::string>::const_iterator& end)
{
  if (begin == end)
    myData = 0.0;
  else if (!parseNumericField(myData,*begin++))
    return false;

  return true;
}


template<> inline bool
FFlField<FaVec3>::parse(std::vector<std::string>::const_iterator& begin,
			const std::vector<std::string>::const_iterator& end)
{
  for (int i = 0; i < 3; i++)
    if (begin == end)
      myData[i] = 0.0;
    else if (!parseNumericField(myData[i],*begin++))
      return false;

  return true;
}


template<> inline bool
FFlField<int>::parse(std::vector<std::string>::const_iterator& begin,
		     const std::vector<std::string>::const_iterator& end)
{
  if (begin == end)
    myData = 0;
  else if (!parseNumericField(myData,*begin++))
    return false;

  return true;
}


template<> inline bool
FFlField<unsigned short>::parse(std::vector<std::string>::const_iterator& begin,
				const std::vector<std::string>::const_iterator& end)
{
  int intVal = 0;
  if (begin != end)
    if (!parseNumericField(intVal,*begin++))
      return false;

  myData = (unsigned short)intVal;
  return true;
}


template<> inline bool
FFlField<std::string>::parse(std::vector<std::string>::const_iterator& begin,
			     const std::vector<std::string>::const_iterator& end)
{
  if (begin != end)
    myData = *begin++;
  else
    myData = "";

  return true;
}


template<> inline bool
FFlField<bool>::parse(std::vector<std::string>::const_iterator& begin,
		      const std::vector<std::string>::const_iterator& end)
{
  if (begin == end)
    myData = false;
  else
  {
    myData = (*begin == "1" || *begin == "true" || *begin == "TRUE");
    begin++;
  }

  return true;
}


template<> inline bool
FFlField< std::vector<double> >::parse(std::vector<std::string>::const_iterator& begin,
				       const std::vector<std::string>::const_iterator& end)
{
  double dVal = 0.0;
  while (begin != end)
    if (!parseNumericField(dVal,*begin++))
      return false;
    else
      myData.push_back(dVal);

  return true;
}


template<> inline void
FFlField< std::vector<double> >::write(std::ostream& os) const
{
  if (myData.size() > 0) os << myData[0];
  for (size_t i = 1; i < myData.size(); i++)
    os << " " << myData[i];
}


inline bool dEqual(const double a, const double b, const double tol = 1.0e-10)
{
  return fabs(a-b) <= tol*(fabs(a) > fabs(b) ? fabs(a) : fabs(b));
}


template<> inline bool
FFlField<double>::isEqual(const FFlField<double>& ref) const
{
  return dEqual(myData,ref.myData);
}


template<> inline bool
FFlField<FaVec3>::isEqual(const FFlField<FaVec3>& ref) const
{
  if (dEqual(myData[0],ref.myData[0]))
    if (dEqual(myData[1],ref.myData[1]))
      return dEqual(myData[2],ref.myData[2]);

  return false;
}


template<> inline bool
FFlField< std::vector<double> >::isEqual(const FFlField< std::vector<double> >& ref) const
{
  size_t n = myData.size();
  if (n != ref.myData.size()) return false;

  for (size_t i = 0; i < n; i++)
    if (!dEqual(myData[i],ref.myData[i])) return false;

  return true;
}

#endif
