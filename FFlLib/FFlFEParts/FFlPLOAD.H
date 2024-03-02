// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PLOAD_H
#define FFL_PLOAD_H

#include "FFlLib/FFlLoadBase.H"
#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFlLib/FFlField.H"

#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlPLOAD : public FFlLoadBase
{
protected:
  FFlPLOAD(int ID);
  FFlPLOAD(const FFlPLOAD& obj);
  FFlPLOAD& operator=(const FFlPLOAD&) = delete;

public:
  virtual ~FFlPLOAD() {}

  virtual void setValue(const std::vector<double>& p) { P.setValue(p); }
  virtual void setValue(const FaVec3&) {}

  virtual void convertUnits(const FFaUnitCalculator* convCal);

  virtual int getLoad(std::vector<FaVec3>& p, int& face) const;

protected:
  virtual const FFlElementBase* getTarget(int& face) const = 0;

private:
  FFlField< std::vector<double> > P;
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


typedef FFaSingelton<FFlTypeInfoSpec,FFlSURFLOAD> FFlSURFLOADTypeInfoSpec;
typedef FFaSingelton<FFlTypeInfoSpec,FFlFACELOAD> FFlFACELOADTypeInfoSpec;

typedef FFaSingelton<FFlFEAttributeSpec,FFlSURFLOAD> FFlSURFLOADAttributeSpec;
typedef FFaSingelton<FFlFEAttributeSpec,FFlFACELOAD> FFlFACELOADAttributeSpec;

#endif
