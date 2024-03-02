// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_LOAD_BASE_H
#define FFL_LOAD_BASE_H

#include "FFlLib/FFlPartBase.H"
#include "FFlLib/FFlFEAttributeRefs.H"
#include "FFaLib/FFaPatterns/FFaGenericFactory.H"

class FFlFieldBase;
class FFlElementBase;
class FFlNode;
class FFaCheckSum;
class FFaUnitCalculator;
class FaVec3;


class FFlLoadBase : public FFlPartBase,
		    public FFlFEAttributeRefs
{
protected:
  FFlLoadBase(int ID) : FFlPartBase(ID) {}
  FFlLoadBase(const FFlLoadBase& obj) : FFlPartBase(obj),
					FFlFEAttributeRefs(obj) {}
  FFlLoadBase& operator=(const FFlLoadBase&) = delete;

public:
  virtual ~FFlLoadBase() {}

  const std::string& getTypeName() const;

  virtual FFlLoadBase* clone() const = 0;

  virtual bool resolveElmRef(const std::vector<FFlElementBase*>&,
			     bool = false) { return true; }
  virtual bool resolveNodeRef(const std::vector<FFlNode*>&,
			      bool = false) { return true; }

  std::vector<FFlFieldBase*>::const_iterator begin() const { return myFields.begin(); }
  std::vector<FFlFieldBase*>::const_iterator end()   const { return myFields.end(); }

  virtual void setValue(const std::vector<double>& p) = 0;
  virtual void setValue(const FaVec3& p) = 0;

  virtual void setTarget(const std::vector<int>& t) = 0;
  virtual void setTarget(int, int = 0) {}
  virtual bool getTarget(int& eid, int& face) const = 0;
  virtual size_t getTargetCount() const = 0;

  virtual int getLoad(std::vector<FaVec3>& p, int& face) const = 0;

  virtual void calculateChecksum(FFaCheckSum* cs, int csMask = 0) const;

  virtual void convertUnits(const FFaUnitCalculator* convCal) = 0;

protected:
  void addField(FFlFieldBase& field) { myFields.push_back(&field); }

private:
  std::vector<FFlFieldBase*> myFields;
};


typedef FFaGenericFactory<FFlLoadBase> LoadFactory;


#define FFL_LOAD_FACTORY_INIT(clName)			\
public:							\
  static void create(int newID, FFlLoadBase*& created)	\
  { created = new clName(newID); }			\
  clName& operator=(const clName&) = delete;		\
  virtual FFlLoadBase* clone() const			\
  { return new clName(*this); }

#endif
