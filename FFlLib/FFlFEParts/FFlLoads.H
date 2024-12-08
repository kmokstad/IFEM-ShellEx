// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_LOADS_H
#define FFL_LOADS_H

#include "FFlLib/FFlLoadBase.H"
#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


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
  FFL_TYPE_INFO(FFlCMOMENT);
};


class FFlPLOAD : public FFlLoadBase
{
  using DoubleVec = std::vector<double>;

protected:
  FFlPLOAD(int ID);
  FFlPLOAD(const FFlPLOAD& obj);
  FFlPLOAD& operator=(const FFlPLOAD&) = delete;

public:
  virtual ~FFlPLOAD() {}

  virtual void setValue(const DoubleVec& p) { P.setValue(p); }
  virtual void setValue(const FaVec3&) {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  virtual int getLoad(std::vector<FaVec3>& p, int& face) const;

protected:
  using FFlLoadBase::getTarget;
  virtual const FFlElementBase* getTarget(int& face) const = 0;

private:
  FFlField<DoubleVec> P;
};


class FFlSURFLOAD : public FFlPLOAD
{
public:
  FFlSURFLOAD(int ID = 0) : FFlPLOAD(ID) { target_counter = 0; }
  FFlSURFLOAD(const FFlSURFLOAD& obj);
  virtual ~FFlSURFLOAD() {}

  virtual bool resolveElmRef(const std::vector<FFlElementBase*>& possibleElms,
			     bool suppessErrmsg = false);

  virtual void setTarget(const std::vector<int>& t);
  virtual bool getTarget(int& eid, int& face) const;
  virtual size_t getTargetCount() const { return target.size(); }

  virtual void calculateChecksum(FFaCheckSum* cs, int csMask) const;

  static void init();

  FFL_LOAD_FACTORY_INIT(FFlSURFLOAD);
  FFL_FE_ATTRIBUTE_SPEC(FFlSURFLOAD);
  FFL_TYPE_INFO(FFlSURFLOAD);

protected:
  virtual const FFlElementBase* getTarget(int& face) const;

private:
  mutable unsigned int target_counter;
  std::vector< FFlReference<FFlElementBase> > target;
};


class FFlFACELOAD : public FFlPLOAD
{
public:
  FFlFACELOAD(int ID = 0) : FFlPLOAD(ID) { target_counter = 0; }
  FFlFACELOAD(const FFlFACELOAD& obj);
  virtual ~FFlFACELOAD() {}

  virtual bool resolveElmRef(const std::vector<FFlElementBase*>& possibleElms,
			     bool suppessErrmsg = false);

  virtual void setTarget(int eid, int face = 0);
  virtual void setTarget(const std::vector<int>& t);
  virtual bool getTarget(int& eid, int& face) const;
  virtual size_t getTargetCount() const { return target.size(); }

  void setFace(int faceNum);

  virtual void calculateChecksum(FFaCheckSum* cs, int csMask) const;

  static void init();

  FFL_LOAD_FACTORY_INIT(FFlFACELOAD);
  FFL_FE_ATTRIBUTE_SPEC(FFlFACELOAD);
  FFL_TYPE_INFO(FFlFACELOAD);

protected:
  virtual const FFlElementBase* getTarget(int& face) const;

private:
  mutable unsigned int target_counter;
  std::vector< std::pair<FFlReference<FFlElementBase>,char> > target;
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
