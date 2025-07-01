// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_FEELEMENT_TOP_SPEC_H
#define FFL_FEELEMENT_TOP_SPEC_H

#include <vector>
#include <utility>


#define FFL_FE_ELEMENT_TOP_SPEC(Class)					\
  virtual FFlFEElementTopSpec* getFEElementTopSpec() const		\
  { return FFaSingelton<FFlFEElementTopSpec,Class>::instance(); }


using EdgeType = std::pair<int,int>;


struct FaceType
{
  void addEdge(int i1, int i2) { myEdges.emplace_back(i1,i2); }
  std::vector<EdgeType> myEdges;
};


class FFlFEElementTopSpec
{
public:
  FFlFEElementTopSpec() : myExplEdgePattern(0xf0f0), myNodeDOFs(0)
  {
    myNodeCount = myExpandedNodeCount = 0;
    iAmShellFaces = allSlaveNodes = false;
  }
  virtual ~FFlFEElementTopSpec() {}

  void addFace(const FaceType& f)         { myFaces.push_back(f); }
  void addExpandedFace(const FaceType& f) { myExpandedFaces.push_back(f); }
  void addExplicitEdge(int i1, int i2) { myExplicitEdges.emplace_back(i1,i2); }

  int getNodeCount() const           { return myNodeCount; }
  int getExpandedNodeCount() const   { return myExpandedNodeCount; }
  virtual int getNodeDOFs(int) const { return myNodeDOFs; }

  virtual bool isSlaveNode(int) const     { return allSlaveNodes; }
  virtual bool isReferenceNode(int) const { return false; }

  bool isShellFaces() const { return iAmShellFaces; }

  void setNodeCount(int count)         { myNodeCount = count; }
  void setExpandedNodeCount(int count) { myExpandedNodeCount = count; }
  void setNodeDOFs(int count)          { myNodeDOFs = count; }

  void setShellFaces(bool trueOrFalse)  { iAmShellFaces = trueOrFalse; }
  void setSlaveStatus(bool trueOrFalse) { allSlaveNodes = trueOrFalse; }

  void setTopology(int faceCount, int* faceConnections, bool expanded = false);
  void setExpandedTopology(int faceCount, int* faceConnections);

  // Get local element node ids for each vertex in face
  bool getFaceTopology(short int faceNumber, bool isExpandedFace,
                       bool switchNormal, short int idxOffset,
                       std::vector<int>& topology) const;

  // Get the number of vertices in the given local face
  int getFaceVertexCount(short int faceNumber, bool isExpandedFace) const;

  // Get local face number defined by two local node ids
  short int getFaceNum(int node1, int node2) const;

  void dump() const;

  std::vector<FaceType> myFaces;
  std::vector<FaceType> myExpandedFaces;
  std::vector<EdgeType> myExplicitEdges;

  unsigned short myExplEdgePattern;

protected:
  int  myNodeCount;
  int  myExpandedNodeCount;
  char myNodeDOFs    : 7;
  bool iAmShellFaces : 1;
  bool allSlaveNodes : 1;
};

#endif
