// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_NASTRAN_READER_H
#define FFL_NASTRAN_READER_H

#include "FFlReaderBase.H"
#include "FFaLib/FFaAlgebra/FFaMat34.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include <string>
#include <vector>
#include <map>
#include <set>


#define DEF_PROCESSOR(rec) bool process_##rec (std::vector<std::string>& entry);

class FFlAttributeBase;
class FFlElementBase;
class FFlLoadBase;
class FFlGroup;
class FFaCheckSum;


class FFlNastranReader : public FFlReaderBase
{
public:

  enum FieldFormat
  {
    UNDEFINED   =  0,
    CONT_FIELD  = -2,
    FREE_FIELD  = -1,
    SMALL_FIELD =  8,
    LARGE_FIELD = 16
  };

  enum coordSysType
  {
    rectangular,
    cylindrical,
    spherical
  };

  struct BulkEntry
  {
    std::string name;
    FieldFormat ffmt;
    std::vector<std::string> fields;
    std::string cont;
  };

private:

  struct CORD
  {
    bool isComputed;
    coordSysType type;
    int RID;
    int G[3];
    FaVec3 Origo;
    FaVec3 Zaxis;
    FaVec3 XZpnt;
    FaMat34 Tmat;
    CORD() : isComputed(false), type(rectangular) { RID = G[0]=G[1]=G[2] = 0; }
  };

  struct GRDSET
  {
    int CP;
    int CD;
    int PS;
    int SEID;
    GRDSET() { CP = CD = PS = SEID = 0; }
  };

  struct BEAMOR
  {
    int PID;
    int G0;
    FaVec3 X;
    bool isBAR;
    bool empty[5];
    BEAMOR(bool b = false) : isBAR(b) { PID = G0 = 0; }
  };

  typedef std::map<int,int>  IntMap;
  typedef std::pair<int,int> IntPair;

  std::map<std::string,int> ignoredBulk; // Number of unsupported entries
  std::map<std::string,int> sxErrorBulk; // Number of syntax errors

  std::vector<BulkEntry> ucEntries; // List of uncompleted bulk entries
  std::set<int> beamPIDnsm;  // Set of PIDs for beam elements with NSM property
  std::set<int> shellPIDnsm; // Set PIDs for shell elements with NSM property

  std::map<int,CORD*> cordSys; // coordinate system ID --> coordinate system
  IntMap nodeStat;  // node ID --> dof status (exteral/internal)
  IntMap nodeCPID;  // node ID --> coordinate system ID for positions
  IntMap nodeCDID;  // node ID --> coordinate system ID for displacements
  IntMap massCID;   // element ID --> coordinate system ID for mass matrix
  IntMap solidPID;  // element ID --> property ID for solid elements
  IntMap shellPID;  // element ID --> property ID for shell elements
  IntMap sprPID;    // element ID --> property ID for spring elements
  IntMap sprComp;   // element ID --> component number for spring elements
  IntMap weldGS;    // element ID --> location node ID for weld elements
  std::map<int,double> sprK;   // element ID --> spring stiffness
  std::map<int,double> propK;  // property ID --> spring stiffness
  std::map<int,IntMap> propMID;// element category -> property ID -> material ID
  std::map<int,BEAMOR*> bOri;  // element ID --> Bar/Beam-orientation vectors
  std::map<int,FaVec3*> massX; // element ID --> offset distance for mass matrix
  GRDSET* gridDefault;    // Default options for GRID entries
  BEAMOR* barDefault;     // Default property values for CBAR entries
  BEAMOR* beamDefault;    // Default property values for CBEAM entries
  std::set<int> PTHICKs;  // Set of PIDs that are of PTHICK type
  std::set<int> PCOMPs;   // Set of PIDs that are of PCOMP type

  std::vector<FFlElementBase*> myWelds; // Temporary vector of all weld elements
  std::vector<std::map<int,FFlElementBase*> > weld; // Weld connectivities

  typedef std::map<FFlLoadBase*,IntPair> LoadFaceMap;
  LoadFaceMap                loadFace; // Pressure load --> element face nodes
  std::map<FFlLoadBase*,int> loadCID; // Load --> coord. system ID for direction

  FFaCheckSum*                             attChkSum;
  std::map<unsigned int,FFlAttributeBase*> uniqueAtts;

  static int nWarnings, nNotes; // Accumulates the number of warnings and notes

  mutable std::pair<int,std::string> lastComment;

  int lineCounter;
  bool sizeOK;

public:

  FFlNastranReader (FFlLinkHandler* readIntoLink, const int startHere);
  FFlNastranReader (const FFlNastranReader&) = delete;
  virtual ~FFlNastranReader ();

  FFlNastranReader& operator= (const FFlNastranReader&) = delete;

  static void init ();
  static void identifierCB (const std::string& fileName, int& positiveIdent);
  static void readerCB     (const std::string& fileName, FFlLinkHandler* link);

  static bool extractNameFromComment (std::string& commentLine,
                                      bool first = false);

protected:

  bool read (const std::string& fileName, bool includedFile = false);
  bool read (std::istream& is);
  bool resolve (bool stillOk);

  void processAssignFile (const std::string& line);

  // Detection and parsing of Nastran sets
  bool setsArePresent (const std::string& fileName);
  bool processSet     (const std::string& fileName, const int startBulk);
  FFlGroup* processThisSet (std::string& setLine,
			    const int startL, const int stopL);

  // Read next bulk-entry from the input stream
  bool getNextEntry (std::istream& is, BulkEntry& entry);
  // Read data fields for the given entry from the input stream
  bool getFields (std::istream& is, BulkEntry& entry);

  // Read next (name- or data-) field from the input stream
  int getNextField (std::istream& is, std::string& field,
		    const FieldFormat size = UNDEFINED);

  // Parse the given entry and put the data into the FFlLinkHandler object
  bool processThisEntry (BulkEntry& entry);

  // Methods for resolving references after finished parsing
  // These methods are defined in FFlNastranResolver.C

  bool resolveCoordinates ();
  bool transformNode (IntMap::iterator& nit);
  bool transformPoint (FaVec3& X, const int CID, bool orientationOnly = false);
  bool computeTmatrix (const int CID, CORD& cs);
  bool getTmatrixAtPt (const int CID, CORD& cs, const FaVec3& X, FaMat33& T);
  bool resolveAttributes ();

  int  resolveBeamAttributes   (FFlElementBase* theBeam, bool& ok);
  bool resolveWeldElement      (FFlElementBase* theElem, int& newEID,
				const int PID);
  bool getElementAxis          (FFlElementBase* theElem,
				const int n1, const int n2, FaVec3& axis);
  void resolveSpringAttributes (FFlElementBase* theSpring,
				const double K, const int C, bool& ok);
  void resolveBushAttributes   (FFlElementBase* theSpring,
			        const double S, const int CID, bool& ok);

  bool transformMassMatrix  (FFlElementBase* theMass, const int CID, FaVec3* X);
  bool transformSymmMatrix6 (double mat[6][6], const FaVec3& X,
			     const int CID, const int node = 0);

  bool transformVec3 (FaVec3& v, const FaVec3& X, const int CID);

  bool resolveLoads ();
  bool resolveLoadFace (const FFlLoadBase* load, const IntPair& nodes) const;
  bool resolveLoadDirection (const FFlLoadBase* load, int CID);

  // The remaining methods are defined in FFlNastranProcessor.C

  bool processThisEntry (const std::string& name,
			 std::vector<std::string>& entry);

  bool insertBeamPropMat (const char* bulk, int PID, int MID);

  int addUniqueAttribute (FFlAttributeBase* newAtt);

  // Alphabetical list of bulk entries that will be processed by this reader:
  DEF_PROCESSOR( ASET   ) // Specification of external nodes
  DEF_PROCESSOR( ASET1  ) // Specification of external nodes, alternate form
  DEF_PROCESSOR( BAROR  ) // Default orientation and properties for CBAR
  DEF_PROCESSOR( BEAMOR ) // Default orientation and properties for CBEAM
  DEF_PROCESSOR( CBAR   ) // Prismatic beam element (beam without tapering)
  DEF_PROCESSOR( CBEAM  ) // General beam element
  DEF_PROCESSOR( CBUSH  ) // Generalized spring connection
  DEF_PROCESSOR( CHEXA  ) // Hexahedron solid element
  DEF_PROCESSOR( CELAS1 ) // Scalar spring element
  DEF_PROCESSOR( CELAS2 ) // Scalar spring element without reference to PELAS
  DEF_PROCESSOR( CONROD ) // Rod element with only axial and torsional stiffness
  DEF_PROCESSOR( CONM1  ) // Concentrated mass element, 6x6 mass matrix
  DEF_PROCESSOR( CONM2  ) // Concentrated mass element
  DEF_PROCESSOR( CORD1C ) // Cylindrical coordinate system
  DEF_PROCESSOR( CORD1R ) // Rectangular coordinate system
  DEF_PROCESSOR( CORD1S ) // Spherical coordinate system
  DEF_PROCESSOR( CORD2C ) // Cylindrical coordinate system
  DEF_PROCESSOR( CORD2R ) // Rectangular coordinate system
  DEF_PROCESSOR( CORD2S ) // Spherical coordinate system
  DEF_PROCESSOR( CPENTA ) // Pentahedron solid element
  DEF_PROCESSOR( CQUAD4 ) // 4-noded shell element
  DEF_PROCESSOR( CQUAD8 ) // 8-noded shell element
  DEF_PROCESSOR( CROD   ) // Rod element with only axial and torsional stiffness
  DEF_PROCESSOR( CTETRA ) // Tetrahedron solid element
  DEF_PROCESSOR( CTRIA3 ) // 3-noded shell element
  DEF_PROCESSOR( CTRIA6 ) // 6-noded shell element
  DEF_PROCESSOR( CWELD  ) // Spot weld element
  DEF_PROCESSOR( FORCE  ) // Concentrated force at grid point
  DEF_PROCESSOR( GRDSET ) // Default coordinate system IDs
  DEF_PROCESSOR( GRID   ) // Nodal Grid point
  DEF_PROCESSOR( INCLUDE) // Include another bulk-data file
  DEF_PROCESSOR( MAT1   ) // Material properties
  DEF_PROCESSOR( MAT2   ) // Material properties
  DEF_PROCESSOR( MAT8   ) // Material properties
  DEF_PROCESSOR( MAT9   ) // Material properties
  DEF_PROCESSOR( MOMENT ) // Concentrated moment at grid point
  DEF_PROCESSOR( PBAR   ) // Properties for CBAR
  DEF_PROCESSOR( PBARL  ) // Properties for CBAR with cross section dimensions
  DEF_PROCESSOR( PBEAM  ) // Properties for CBEAM
  DEF_PROCESSOR( PBEAML ) // Properties for CBEAM with cross section dimensions
  DEF_PROCESSOR( PBUSH  ) // Properties for CBUSH
  DEF_PROCESSOR( PCOMP  ) // Properties for composite shell
  DEF_PROCESSOR( PELAS  ) // Properties for CELAS1
  DEF_PROCESSOR( PLOAD2 ) // Uniform normal pressure load on a surface element
  DEF_PROCESSOR( PLOAD4 ) // Pressure load on surfaces and solid element faces
  DEF_PROCESSOR( PROD   ) // Properties for CROD
  DEF_PROCESSOR( PSHELL ) // Properties for shell elements
  DEF_PROCESSOR( PSOLID ) // Properties for solid elements
  DEF_PROCESSOR( PWELD  ) // Properties for spot weld element
  DEF_PROCESSOR( QSET1  ) // Specification of component modes
  DEF_PROCESSOR( RBAR   ) // Rigid bar element
  DEF_PROCESSOR( RBE2   ) // Rigid body with any number of dependent nodes
  DEF_PROCESSOR( RBE3   ) // Motion of a reference node as a weighted average
  DEF_PROCESSOR( SPC    ) // Single-point constraints
  DEF_PROCESSOR( SPC1   ) // Single-point constraints, alternate form
};

#endif
