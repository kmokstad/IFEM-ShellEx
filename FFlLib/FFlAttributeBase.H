// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_ATTRIBUTE_BASE_H
#define FFL_ATTRIBUTE_BASE_H

#include "FFlLib/FFlNamedPartBase.H"
#include "FFlLib/FFlFEAttributeRefs.H"
#include "FFaLib/FFaPatterns/FFaGenericFactory.H"

class FFlFieldBase;
class FFaCheckSum;
class FFaUnitCalculator;


class FFlAttributeBase : public FFlNamedPartBase,
			 public FFlFEAttributeRefs
{
protected:
  FFlAttributeBase(int ID = 0) : FFlNamedPartBase(ID) {}
  FFlAttributeBase(const FFlAttributeBase& obj) : FFlNamedPartBase(obj),
						  FFlFEAttributeRefs(obj) {}
  FFlAttributeBase& operator=(const FFlAttributeBase&) = delete;

public:
  virtual ~FFlAttributeBase() {}

  const std::string& getTypeName() const;
  const std::string& getDescription() const;

  size_t size() const { return myFields.size(); }

  std::vector<FFlFieldBase*>::const_iterator begin() const { return myFields.begin(); }
  std::vector<FFlFieldBase*>::const_iterator end()   const { return myFields.end(); }

  void print(const char* prefix = "Attribute ") const;

  virtual void calculateChecksum(FFaCheckSum* cs, int csMask) const;

  virtual FFlAttributeBase* clone() const = 0;

  virtual void convertUnits(const FFaUnitCalculator*) {}

  virtual bool isIdentic(const FFlAttributeBase* other) const;

  virtual void resize(size_t n) { if (n < myFields.size()) myFields.resize(n); }

protected:
  int getAttributeTypeID() const;

  void addField(FFlFieldBase& field) { myFields.push_back(&field); }

private:
  std::vector<FFlFieldBase*> myFields;
};


typedef FFaGenericFactory<FFlAttributeBase> AttributeFactory;


#define FFL_FE_ATTRIBUTE_FACTORY_INIT(clName)			\
public:								\
  static void create(int newID, FFlAttributeBase*& created)	\
  { created = new clName(newID); }				\
  clName& operator=(const clName&) = delete;			\
  virtual FFlAttributeBase* clone() const 	     	 	\
  { return new clName(*this); }

#endif
