// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_LINK_HANDLER_H
#define FFL_LINK_HANDLER_H

#include <map>
#include <set>
#include <vector>
#include <string>

#include "FFlTypeInfoSpec.H"
#include "FFlLinkCSMask.H"

class FFlGroup;
class FFlNode;
class FFlVertex;
class FFlPartBase;
class FFlNamedPartBase;
class FFlElementBase;
class FFlLoadBase;
class FFlAttributeBase;
class FFlVisualBase;
class FFlVDetail;
class FFlVAppearance;
class FFlFEResultBase;
class FFlrFELinkResult;
class FFlFaceGenerator;
class FFlConnectorItems;
class FFaCompoundGeometry;
class FFaUnitCalculator;
class FFaProfiler;
class FFaCheckSum;
class FFaTensor3;
class FaMat34;
class FaVec3;

typedef std::set<FFlElementBase*>        ElementsSet;
typedef std::vector<FFlElementBase*>     ElementsVec;
typedef ElementsVec::const_iterator      ElementsCIter;
typedef std::map<std::string,size_t>     ElmTypeCount;

typedef std::vector<FFlNode*>            NodesVec;
typedef NodesVec::const_iterator         NodesCIter;

#ifdef FT_USE_VERTEX
typedef std::vector<FaVec3*>             VertexVec;
typedef VertexVec::const_iterator        VertexCIter;
#endif

typedef std::map<int,FFlGroup*>          GroupMap;
typedef GroupMap::const_iterator         GroupCIter;

typedef std::vector<FFlLoadBase*>        LoadsVec;
typedef LoadsVec::const_iterator         LoadsCIter;

typedef std::map<int,FFlAttributeBase*>    AttributeMap;
typedef std::map<std::string,AttributeMap> AttributeTypeMap;
typedef AttributeTypeMap::const_iterator   AttributeTypeCIter;

#ifdef FT_USE_VISUALS
typedef std::vector<FFlVisualBase*>        VisualsVec;
typedef VisualsVec::const_iterator         VisualsCIter;
#endif

typedef std::pair<FFlElementBase*,int> FFlrElement;
typedef std::vector<FFlrElement>       FFlrElementVec;
typedef std::vector<FFlrElementVec>    FFlrVxToElmMap;


class FFlLinkHandler
{
public:
  FFlLinkHandler(size_t maxNodes = 0, size_t maxElms = 0);
  FFlLinkHandler(const FFlLinkHandler& otherLink);
  FFlLinkHandler(const FFlGroup& fromGroup);
  ~FFlLinkHandler();

  FFlLinkHandler& operator=(const FFlLinkHandler&) = delete;

  // Adding entities
  bool addElement(FFlElementBase* anElem, bool sortOnInsert = false);
  bool addNode(FFlNode* aNode, bool sortOnInsert = false);
  bool addGroup(FFlGroup* group, bool silence = false);
  void addLoad(FFlLoadBase* load, bool sortOnInsert = false);
  bool addAttribute(FFlAttributeBase* attr, bool silence = false);
  bool addAttribute(FFlAttributeBase* attr, bool silence, const std::string& name);
  int  addUniqueAttribute(FFlAttributeBase* newAtt, bool silence = false);
  bool removeAttribute(const std::string& typeName, int ID, bool silence = false);
#ifdef FT_USE_VISUALS
  void addVisual(FFlVisualBase* visual, bool sortOnInsert = false);
  void setRunningIdxOnAppearances();
#endif

  int  sortElementsAndNodes(bool deleteDuplicates = false) const;

  // Getting count of entities

  enum { FFL_ALL = 0, FFL_FEM = 1, FFL_STRC = 2 };

  const ElmTypeCount& getElmTypeCount() const;
  int getElementCount(const std::string& type) const;
  int getElementCount(int types = FFL_ALL, bool checkCF = false) const;
  int getNodeCount(int types = FFL_ALL) const;
  int getDofCount(bool includeExternalDofs = true) const;
  int getAttributeCount(const std::string& type) const;
  int buildFiniteElementVec(bool allFElements = true) const;

  // Looking up entities by ID (external or internal)

  FFlElementBase*   getElement(int ID, bool internalID = false) const;
  FFlElementBase*   getFiniteElement(int iel) const;
  FFlNode*          getNode(int ID) const;
  FFlNode*          getFENode(int inod) const;
  FFlGroup*         getGroup(int ID) const;
  FFlAttributeBase* getAttribute(const std::string& type, int ID) const;
  FFlVDetail*       getDetail(int ID) const;
  FFlVAppearance*   getAppearance(int ID) const;

  bool getRefNodes(std::vector<FFlNode*>& refNodes) const;
  bool getLoads(int ID, LoadsVec& loads) const;
  void getLoadCases(std::set<int>& IDs) const;

  // External-to-internal node and (finite) element number mapping
  int getIntNodeID(int extID) const;
  int getIntElementID(int extID) const;

  int getNewElmID() const;
  int getNewNodeID() const;
  int getNewGroupID() const;
  int getNewAttribID(const std::string& type) const;
#ifdef FT_USE_VISUALS
  int getNewVisualID() const;
#endif

  // Begin iterators to the containers

  ElementsCIter elementsBegin()  const { return myElements.begin(); }
  ElementsCIter fElementsBegin() const { return myFElements.begin(); }
  NodesCIter    nodesBegin()     const { return myNodes.begin(); }
#ifdef FT_USE_VERTEX
  VertexCIter   verticesBegin()  const { return myVertices.begin(); }
#endif
  GroupCIter    groupsBegin()    const { return myGroupMap.begin(); }
  LoadsCIter    loadsBegin()     const { return myLoads.begin(); }
#ifdef FT_USE_VISUALS
  VisualsCIter  visualsBegin()   const { return myVisuals.begin(); }
#endif
  AttributeTypeCIter attributeTypesBegin() const { return myAttributes.begin(); }

  // End iterators to the containers

  ElementsCIter elementsEnd()  const { return myElements.end(); }
  ElementsCIter fElementsEnd() const { return myFElements.end(); }
  NodesCIter    nodesEnd()     const { return myNodes.end(); }
#ifdef FT_USE_VERTEX
  VertexCIter   verticesEnd()  const { return myVertices.end(); }
#endif
  GroupCIter    groupsEnd()    const { return myGroupMap.end(); }
  LoadsCIter    loadsEnd()     const { return myLoads.end(); }
#ifdef FT_USE_VISUALS
  VisualsCIter  visualsEnd()   const { return myVisuals.end(); }
#endif
  AttributeTypeCIter attributeTypesEnd() const { return myAttributes.end(); }

  const AttributeMap& getAttributes(const std::string& name) const;

  // Special vertex management access :

#ifdef FT_USE_VERTEX
  size_t getVertexCount() const { return myVertices.size(); }
  const VertexVec& getVertexes() const { return myVertices; }
  const FaVec3& getVertexPos(size_t idx) const { return *myVertices[idx]; }
  FFlVertex* getVertex(size_t idx) const;
  void addVertex(FFlVertex* aVertex);

  // For each vertex, a vector of elements with the local index of the vertex
  const FFlrVxToElmMap& getVxToElementMapping();
#endif

  void getAllInternalCoordSys(std::vector<FaMat34>& mxes) const;

  // Utilities :

  bool getExtents(FaVec3& max, FaVec3& min) const;
  double getMeanElementSize() const;

  void getMassProperties(double& M, FaVec3& Xcg, FFaTensor3& I) const;

  NodesCIter findFreeNodeAtPoint(const FaVec3& point, double tol,
				 int dofFilter = -1) const;
  NodesCIter findClosestNode(const FaVec3& point) const;

  FFlNode* createAttachableNode(FFlNode* fromNode, const FaVec3& nodePos,
				FFlConnectorItems& cItems);

  ElementsCIter findClosestElement(const FaVec3& point,
				   const std::vector<FFlTypeInfoSpec::Cathegory>&) const;
  FFlElementBase* findClosestElement(const FaVec3& point, const FFlGroup& group) const;
  FFlElementBase* findClosestElement(const FaVec3& point, FFlGroup* group = NULL) const;

  FFlElementBase* findPoint(const FFlGroup& group, const FaVec3& point, double* xi) const;
  FFlElementBase* findPoint(const FaVec3& point, double* xi, int groupID = 0) const;

  // Strain coat creation (these methods are defined in FFlStrainCoatCreator.C)

  bool makeStrainCoat(FFlFaceGenerator* geometry,
                      FFlNamedPartBase* group = NULL);
  bool assignFatigueProperty(int stdIndx, int curveIndx, double sCF,
                             FFlNamedPartBase* group = NULL);

  static bool isStrainCoat(FFlElementBase* elm);

  // Result management :

  void deleteResults();
  bool hasResults() const { return myResults != NULL; }
  FFlrFELinkResult* getResults();

  void convertUnits(const FFaUnitCalculator* convCal);

  // Admin :

  bool resolve(bool subdivParabolic = false, bool fromSESAM = false);
  bool verify(bool fixNegElms = false);
  void deleteGeometry();
  bool hasGeometry() const { return this->getNodeCount() > 0; }
  bool isTooLarge() const { return tooLarge; }

  unsigned int calculateChecksum(int csType = FFl::CS_NOEXTINFO,
                                 bool rndOff = true) const;
  void calculateChecksum(FFaCheckSum* cs, bool rndOff = true) const;
  void calculateChecksum(FFaCheckSum* cs, int csType, bool rndOff) const;

  void initiateCalculationFlag(bool status = false);
  bool updateCalculationFlag(int groupId, bool status = true);
  bool updateCalculationFlag(FFlPartBase* part, bool status = true);
  bool updateCalculationFlag(const std::string& attType, int id, bool status = true);

#ifdef FT_USE_VISUALS
  // Visual settings
  bool setVisDetail(const FFlVDetail* det); // on all elements
  bool setVisDetail(FFlPartBase* part, const FFlVDetail* det);
  bool setVisDetail(const std::vector<FFlPartBase*>& parts, const FFlVDetail* det);
  bool setVisAppearance(FFlPartBase* partSpec, const FFlVAppearance* app);
  void updateGroupVisibilityStatus();

  // Predefined details
  FFlVDetail* getPredefDetail(int detailType);
  FFlVDetail* getOnDetail();
  FFlVDetail* getOffDetail();
#endif

  // For external Nastran files
  void addComponentModes(int nGen) { nGenDofs += nGen; }
  void addOP2file(const std::string& fname) { myOP2files.push_back(fname); }
  void clearOP2files() { myOP2files.clear(); }
  const std::vector<std::string>& getOP2files() const { return myOP2files; }
  size_t getNumberOfGenDofs() const { return nGenDofs; }

  void dump() const;

  FFlNode* createNodeAtPoint(const FaVec3& nodePos, int DOFs = 0);

  int createConnector(const FFaCompoundGeometry& compound,
                      const FaVec3& nodePos, int spiderType,
                      FFlConnectorItems& cItems);

  int deleteConnectorProperties(const std::vector<int>& propertiesID);
  int deleteConnectorElements(const std::vector<int>& elementsID);
  int deleteConnectorNodes(const std::vector<int>& nodesID);

protected:
  typedef NodesVec::iterator    NodesIter;
  typedef ElementsVec::iterator ElementsIter;

  NodesIter    getNodeIter(int ID) const;
  ElementsIter getElementIter(int ID) const;

  void countElements() const;

  bool areBUSHconnected(FFlNode* n1, FFlNode* n2) const;
  int  buildBUSHelementSet() const;

  //! \brief Subdivides a parabolic element into linear elements.
  int splitElement(FFlElementBase* elm);

  // Sorting quickies :

  int  sortElements(bool deleteDuplicates = false) const;
  int  sortNodes(bool deleteDuplicates = false) const;
  void sortLoads() const;
#ifdef FT_USE_VISUALS
  void sortVisuals() const;
#endif

  // Model size limits :

  bool   tooLarge;  // set by addNode and addElement
  size_t nodeLimit; // max number of nodes allowed
  size_t elmLimit;  // max number of elements allowed

  // Containers :

  mutable ElementsVec   myElements;
  mutable bool          areElementsSorted;
  mutable NodesVec      myNodes;
  mutable bool          areNodesSorted;
  mutable LoadsVec      myLoads;
  mutable bool          areLoadsSorted;
#ifdef FT_USE_VISUALS
  mutable VisualsVec    myVisuals;
  mutable bool          areVisualsSorted;
#endif

  GroupMap              myGroupMap;
  AttributeTypeMap      myAttributes;
#ifdef FT_USE_VERTEX
  VertexVec             myVertices;
  FFlrVxToElmMap        myVxMapping;
#endif

  bool                  isResolved;    // Flag to avoid resolving more than once
  bool                  hasLooseNodes; // true if unconnected nodes are present
  mutable std::map<int,int> ext2intNode; // external-to-internal node mapping
  mutable NodesVec      myFEnodes;     // pointers to Finite Element nodes only
  mutable ElmTypeCount  myNumElements; // the number of elements of each type
  mutable ElementsVec   myFElements;   // pointers to the Finite Elements only
  mutable ElementsSet   myBushElements;// pointers to the BUSH elements only
  std::vector<std::string> myOP2files; // files with externally reduced matrices
  size_t                nGenDofs;      // as read from Nastran bulk data file

  FFlFEResultBase*      myResults;
  FFaProfiler*          myProfiler;
};

#endif
