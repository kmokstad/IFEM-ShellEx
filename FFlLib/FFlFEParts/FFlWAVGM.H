// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_WAVGM_H
#define FFL_WAVGM_H

#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEElementTopSpec.H"
#include "FFlLib/FFlFEAttributeSpec.H"
#include "FFlLib/FFlTypeInfoSpec.H"
#ifdef FT_USE_MEMPOOL
#include "FFaLib/FFaPatterns/FFaMemPool.H"
#endif

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


class FFlWAVGMTopSpec : public FFlFEElementTopSpec
{
public:
  FFlWAVGMTopSpec() : FFlFEElementTopSpec() {}
  FFlWAVGMTopSpec(const FFlFEElementTopSpec& obj) : FFlFEElementTopSpec(obj) {}
  virtual ~FFlWAVGMTopSpec() {}

  // The first element node is a slave and cannot be external
  virtual bool isSlaveNode(int node) const { return node == 1; }
  virtual bool isReferenceNode(int node) const { return node == 1; }
};


class FFlWAVGM : public FFlElementBase
{
public:
  FFlWAVGM(int ID);
  FFlWAVGM(const FFlWAVGM& obj);
  virtual ~FFlWAVGM();

  static void init();

  virtual bool setNode(const int topologicalPos, FFlNode* nodeRef);
  virtual bool setNode(const int topologicalPos, int nodeRef);
  virtual bool setNodes(const std::vector<int>& nodeRefs,
                        size_t offset, bool = false);
  virtual bool setNodes(const std::vector<FFlNode*>& nodeRefs,
                        size_t offset, bool = false);

  void setSlaveNode(FFlNode* nodeRef);
  void setSlaveNode(int nodeRef);

  void addMasterNode(FFlNode* nodeRef);
  void addMasterNode(int nodeRef);
  void addMasterNodes(const std::vector<int>& nodeRefs);
  void addMasterNodes(const std::vector<FFlNode*>& nodeRefs);
  bool removeMasterNodes(std::vector<int>& nodeRefs);

  FFlNode* getSlaveNode() const;
  void getMasterNodes(std::vector<FFlNode*>& nodeRefs) const;

  // returns a unique topology spec for the WAVGM element
  virtual FFlFEElementTopSpec* getFEElementTopSpec() const;

  FFL_FE_ELEMENT_FACTORY_INIT(FFlWAVGM);
  FFL_FE_ATTRIBUTE_SPEC(FFlWAVGM);
  FFL_TYPE_INFO(FFlWAVGM);

protected:
  FFlFEElementTopSpec* myWAVGMElemTopSpec;

#ifdef FT_USE_MEMPOOL
  FFA_MAKE_MEMPOOL;
#endif
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
