// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_RGD_H
#define FFL_RGD_H

#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif
#include "FFaLib/FFaPatterns/FFaSingelton.H"


class FFlRGDTopSpec : public FFlFEElementTopSpec
{
public:
  FFlRGDTopSpec() : FFlFEElementTopSpec() {}
  FFlRGDTopSpec(const FFlFEElementTopSpec& obj) : FFlFEElementTopSpec(obj) {}
  virtual ~FFlRGDTopSpec() {}

  // The master node always has 6 DOFs.
  // For the slave nodes it depends on the other elements.
  virtual int getNodeDOFs(int node) const { return node == 1 ? myNodeDOFs : 0; }
  virtual bool isSlaveNode(int node) const { return node>1 && !allowSlvAttach; }

  static bool allowSlvAttach; // May be set by command-line handler
};


class FFlRGD : public FFlElementBase
{
public:
  FFlRGD(int ID);
  FFlRGD(const FFlRGD& obj);
  virtual ~FFlRGD();

  static void init();

  virtual bool setNode(const int topologicalPos, FFlNode* nodeRef);
  virtual bool setNode(const int topologicalPos, int nodeRef);
  virtual bool setNodes(const std::vector<int>& nodeRefs, int offset);
  virtual bool setNodes(const std::vector<FFlNode*>& nodeRefs, int offset);

  void setMasterNode(FFlNode* nodeRef);
  void setMasterNode(int nodeRef);

  void addSlaveNode(FFlNode* nodeRef);
  void addSlaveNode(int nodeRef);
  void addSlaveNodes(const std::vector<int>& nodeRefs);
  void addSlaveNodes(const std::vector<FFlNode*>& nodeRefs);

  FFlNode* getMasterNode() const;
  void getSlaveNodes(std::vector<FFlNode*>& nodeRefs) const;

  // returns a unique topology spec for the RGD element
  virtual FFlFEElementTopSpec* getFEElementTopSpec() const;

  FFL_FE_ELEMENT_FACTORY_INIT(FFlRGD);
  FFL_FE_ATTRIBUTE_SPEC(FFlRGD);
  FFL_TYPE_INFO(FFlRGD);

protected:
  FFlFEElementTopSpec* myRGDElemTopSpec;

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};


typedef FFaSingelton<FFlFEAttributeSpec,FFlRGD> FFlRGDAttributeSpec;
typedef FFaSingelton<FFlTypeInfoSpec,FFlRGD>    FFlRGDTypeInfoSpec;

#endif
