// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_CLOAD_H
#define FFL_CLOAD_H

#include "FFlLib/FFlLoadBase.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlCLOAD : public FFlLoadBase
{
protected:
  FFlCLOAD(int ID, const char type);
  FFlCLOAD(const FFlCLOAD& obj);
  FFlCLOAD& operator=(const FFlCLOAD&) = delete;

public:
  virtual ~FFlCLOAD() {}

  virtual bool resolveNodeRef(const std::vector<FFlNode*>& possibleNodes,
			      bool suppessErrmsg = false);

  virtual void setValue(const std::vector<double>&) {}
  virtual void setValue(const FaVec3& p) { P.setValue(p); }

  virtual void setTarget(const std::vector<int>& t);
  virtual void setTarget(int nid, int dummy = 0);
  virtual bool getTarget(int& nid, int& type) const;
  virtual size_t getTargetCount() const { return target.size(); }

  virtual int getLoad(std::vector<FaVec3>& p, int& type) const;

  virtual void calculateChecksum(FFaCheckSum* cs, int csMask) const;

  FFlField<FaVec3> P;

protected:
  const char myType; // 'F'=force, 'M'=moment

private:
  mutable size_t target_counter;
  std::vector< FFlReference<FFlNode> > target;
};


class FFlCFORCE : public FFlCLOAD
{
public:
  FFlCFORCE(int ID = 0) : FFlCLOAD(ID,'F') {}
  FFlCFORCE(const FFlCFORCE& obj) : FFlCLOAD(obj) {}
  virtual ~FFlCFORCE() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_LOAD_FACTORY_INIT(FFlCFORCE);
  FFL_FE_ATTRIBUTE_SPEC(FFlCFORCE);
  FFL_TYPE_INFO(FFlCFORCE);
};


class FFlCMOMENT : public FFlCLOAD
{
public:
  FFlCMOMENT(int ID = 0) : FFlCLOAD(ID,'M') {}
  FFlCMOMENT(const FFlCMOMENT& obj) : FFlCLOAD(obj) {}
  virtual ~FFlCMOMENT() {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  static void init();

  FFL_LOAD_FACTORY_INIT(FFlCMOMENT);
  FFL_FE_ATTRIBUTE_SPEC(FFlCMOMENT);
  FFL_TYPE_INFO(FFlCMOMENT);
};


typedef FFaSingelton<FFlTypeInfoSpec,FFlCFORCE>  FFlCFORCETypeInfoSpec;
typedef FFaSingelton<FFlTypeInfoSpec,FFlCMOMENT> FFlCMOMENTTypeInfoSpec;

typedef FFaSingelton<FFlFEAttributeSpec,FFlCFORCE>  FFlCFORCEAttributeSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlCMOMENT> FFlCMOMENTAttributeSpec;

#endif
