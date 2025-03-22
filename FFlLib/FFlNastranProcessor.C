// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlIOAdaptors/FFlNastranReader.H"
#include "FFlLib/FFlIOAdaptors/FFlReaders.H"
#include "FFlLib/FFlIOAdaptors/FFlCrossSection.H"
#include "FFlLib/FFlLinkHandler.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlLoadBase.H"
#include "FFlLib/FFlGroup.H"
#include "FFlLib/FFlFEParts/FFlPBEAMSECTION.H"
#include "FFlLib/FFlFEParts/FFlPBEAMECCENT.H"
#include "FFlLib/FFlFEParts/FFlPBEAMPIN.H"
#include "FFlLib/FFlFEParts/FFlPTHICK.H"
#include "FFlLib/FFlFEParts/FFlPCOMP.H"
#include "FFlLib/FFlFEParts/FFlPMASS.H"
#include "FFlLib/FFlFEParts/FFlPMAT.H"
#include "FFlLib/FFlFEParts/FFlPNSM.H"
#include "FFlLib/FFlFEParts/FFlPRGD.H"
#include "FFlLib/FFlFEParts/FFlPRBAR.H"
#include "FFlLib/FFlFEParts/FFlPWAVGM.H"
#include "FFlLib/FFlFEParts/FFlPBUSHCOEFF.H"
#include "FFlLib/FFlFEParts/FFlPBUSHECCENT.H"
#include "FFlLib/FFlFEParts/FFlPORIENT.H"
#include "FFaLib/FFaAlgebra/FFaMath.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"
#ifdef FFL_TIMER
#include "FFaLib/FFaProfiler/FFaProfiler.H"
#endif
#include <cstring>

#ifdef FFL_TIMER
#define START_TIMER(func) myProfiler->startTimer(func);
#define STOPP_TIMER(func) myProfiler->stopTimer(func);
#define PROCC_TIMER(func) myProfiler->stopTimer("process_"+std::string(func));
#else
#define START_TIMER(func)
#define STOPP_TIMER(func)
#define PROCC_TIMER(func)
#endif

#define CONVERT_ENTRY(name,ok) if (!(ok)) {\
  if (++sxErrorBulk[name] == 11)\
    ListUI <<"\n            ...\n";\
  else if (sxErrorBulk[name] < 11) {\
    ListUI <<"\n *** Error: Syntax error in "<< name <<" entry: ";\
    for (size_t i = 0; i < entry.size(); i++)\
      ListUI << entry[i] << (i+1 < entry.size() ? ", " : "\n"); }\
  PROCC_TIMER(name)\
  return false; }

#define fieldValue(field,value) FFlFieldBase::parseNumericField(value,field)

// This macro first tries to parse the field as an integer value.
// Only if that fails, it then tries to parse it as a real value.
#define fieldValue2(field,ival,rval) \
  (FFlFieldBase::parseNumericField(ival,field,true) ? true : \
   FFlFieldBase::parseNumericField(rval,field))

#define CREATE_ATTRIBUTE(type,name,id) \
  dynamic_cast<type*>(AttributeFactory::instance()->create(name,id))

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


////////////////////////////////////////////////////////////////////////////////

static FFlAttributeBase* createBeamSection (int PID, FFlCrossSection& data,
                                            std::pair<int,std::string>& comment)
{
  FFlPBEAMSECTION* myAtt = CREATE_ATTRIBUTE(FFlPBEAMSECTION,"PBEAMSECTION",PID);
  myAtt->crossSectionArea = round(data.A,10);
  myAtt->phi = round(data.findMainAxes(),10);
  myAtt->Iy  = round(data.Izz,10);
  myAtt->Iz  = round(data.Iyy,10);
  myAtt->It  = round(data.J,10);
  myAtt->Kxy = round(data.K1,10);
  myAtt->Kxz = round(data.K2,10);
  myAtt->Sy  = round(data.S1,10);
  myAtt->Sz  = round(data.S2,10);

  if (comment.first > 0)
    // Set property name from the first comment line before this PBEAM entry
    if (FFlNastranReader::extractNameFromComment(comment.second,true))
      myAtt->setName(comment.second);

#ifdef FFL_DEBUG
  std::cout << data;
  myAtt->print();
#endif
  return myAtt;
}

////////////////////////////////////////////////////////////////////////////////

static FFlAttributeBase* createNSM (int PID, double NSM, bool isShell = false)
{
  FFlPNSM* myAtt = CREATE_ATTRIBUTE(FFlPNSM,"PNSM",PID);
  myAtt->NSM     = round(NSM,10);
  myAtt->isShell = isShell ? 1 : 0;

#ifdef FFL_DEBUG
  myAtt->print();
#endif
  return myAtt;
}

////////////////////////////////////////////////////////////////////////////////

static FFlElementBase* createElement (const std::string& elmType, int elmId,
                                      std::vector<int>& elmNodes,
                                      int propId = 0, int coordId = 0)
{
#if FFL_DEBUG > 2
  std::cout << elmType <<" element "<< elmId;
  if (propId > 0) std::cout <<", property "<< propId;
  std::cout <<", Nodes:";
  for (int node : elmNodes) std::cout <<" "<< node;
  std::cout << std::endl;
#endif

  FFlElementBase* theElm = ElementFactory::instance()->create(elmType,elmId);
  if (!theElm)
  {
    ListUI <<"\n *** Error: Failure creating element "<< elmId
	   <<" of type "<< elmType <<".\n";
    return theElm;
  }

  theElm->setNodes(elmNodes);
  if (propId > 0)
  {
    if (theElm->getCathegory() == FFlTypeInfoSpec::SHELL_ELM)
      theElm->setAttribute("PTHICK",propId);
    else if (elmType == "BEAM2")
      theElm->setAttribute("PBEAMSECTION",propId);
    else if (elmType == "CMASS")
      theElm->setAttribute("PMASS",propId);
    else if (elmType == "SPRING" || elmType == "RSPRING")
      theElm->setAttribute("PSPRING",propId);
    else if (elmType == "BUSH")
      theElm->setAttribute("PBUSHCOEFF",propId);
    else if (elmType == "RGD")
      theElm->setAttribute("PRGD",propId);
    else if (elmType == "RBAR")
      theElm->setAttribute("PRBAR",propId);
    else if (elmType == "WAVGM")
      theElm->setAttribute("PWAVGM",propId);
  }

  if (coordId > 0)
    if (theElm->getCathegory() == FFlTypeInfoSpec::SHELL_ELM)
      theElm->setAttribute("PCOORDSYS",coordId);

  return theElm;
}

////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::processThisEntry (const std::string& name,
					 std::vector<std::string>& entry)
{
  // Check the most likely entry names first
  if (name == "GRID")   return process_GRID   (entry);
  if (name == "CQUAD4") return process_CQUAD4 (entry);
  if (name == "CTRIA3") return process_CTRIA3 (entry);
  if (name == "CTETRA") return process_CTETRA (entry);
  if (name == "CHEXA")  return process_CHEXA  (entry);
  if (name == "CPENTA") return process_CPENTA (entry);
  if (name == "CQUADR") return process_CQUAD4 (entry); // NB: CQUADR --> CQUAD4
  if (name == "CTRIAR") return process_CTRIA3 (entry); // NB: CTRIAR --> CTRIA3
  if (name == "CQUAD8") return process_CQUAD8 (entry);
  if (name == "CTRIA6") return process_CTRIA6 (entry);
  if (name == "RBE2")   return process_RBE2   (entry);
  if (name == "RBE3")   return process_RBE3   (entry);
  if (name == "RBAR")   return process_RBAR   (entry);
  if (name == "MPC")    return process_MPC    (entry);
  if (name == "CWELD")  return process_CWELD  (entry);
  if (name == "CBEAM")  return process_CBEAM  (entry);
  if (name == "CBAR")   return process_CBAR   (entry);
  if (name == "CROD")   return process_CROD   (entry);
  if (name == "CONROD") return process_CONROD (entry);
  if (name == "CELAS1") return process_CELAS1 (entry);
  if (name == "CELAS2") return process_CELAS2 (entry);
  if (name == "CBUSH")  return process_CBUSH  (entry);
  if (name == "CONM1")  return process_CONM1  (entry);
  if (name == "CONM2")  return process_CONM2  (entry);
  if (name == "PSHELL") return process_PSHELL (entry);
  if (name == "PCOMP")  return process_PCOMP  (entry);
  if (name == "PWELD")  return process_PWELD  (entry);
  if (name == "PBEAM")  return process_PBEAM  (entry);
  if (name == "PBAR")   return process_PBAR   (entry);
  if (name == "PROD")   return process_PROD   (entry);
  if (name == "PBEAML") return process_PBEAML (entry);
  if (name == "PBARL")  return process_PBARL  (entry);
  if (name == "PELAS")  return process_PELAS  (entry);
  if (name == "PBUSH")  return process_PBUSH  (entry);
  if (name == "PSOLID") return process_PSOLID (entry);
  if (name == "PLOAD2") return process_PLOAD2 (entry);
  if (name == "PLOAD4") return process_PLOAD4 (entry);
  if (name == "FORCE")  return process_FORCE  (entry);
  if (name == "MOMENT") return process_MOMENT (entry);
  if (name == "MAT1")   return process_MAT1   (entry);
  if (name == "MAT2")   return process_MAT2   (entry);
  if (name == "MAT8")   return process_MAT8   (entry);
  if (name == "MAT9")   return process_MAT9   (entry);
  if (name == "CORD1R") return process_CORD1R (entry);
  if (name == "CORD2R") return process_CORD2R (entry);
  if (name == "CORD1C") return process_CORD1C (entry);
  if (name == "CORD2C") return process_CORD2C (entry);
  if (name == "CORD1S") return process_CORD1S (entry);
  if (name == "CORD2S") return process_CORD2S (entry);
  if (name == "INCLUDE")return process_INCLUDE(entry);
  if (name == "ASET")   return process_ASET   (entry);
  if (name == "ASET1")  return process_ASET1  (entry);
  if (name == "QSET1")  return process_QSET1  (entry);
  if (name == "SPC")    return process_SPC    (entry);
  if (name == "SPC1")   return process_SPC1   (entry);
  if (name == "SET1")   return process_SET1   (entry);
  if (name == "GRDSET") return process_GRDSET (entry);
  if (name == "BEAMOR") return process_BEAMOR (entry);
  if (name == "BAROR")  return process_BAROR  (entry);
  if (name == "SEQGP")  return true; // Silently ignore this entry
  if (name == "PARAM")  return true; // Silently ignore this entry
  if (name == "param")  return true; // Silently ignore this entry
  if (name == "EIGRL")  return true; // Silently ignore this entry
  if (name == "SPOINT") return true; // Silently ignore this entry
  if (name == "BEGINBU")return true; // Silently ignore this entry
  if (name == "TextInp")return true; // Silently ignore this entry

  // Only print up to 5 warnings of each kind
  int nMsg = ++ignoredBulk[name];
  if (nMsg > 5) return true;

  nWarnings++;
  ListUI <<"\n  ** Warning: Unknown bulk-entry ignored: "<< name;
  if (entry.size() > 4)
  {
    for (size_t i = 0; i < 4; i++)
      ListUI <<", "<< entry[i];
    ListUI <<", ...";
  }
  else
    for (size_t i = 0; i < entry.size(); i++)
      ListUI <<", "<< entry[i];
  ListUI <<" (Line: "<< lineCounter <<").\n";

  if (nMsg == 5)
    ListUI <<"              (subsequent incidences of such entries"
           <<", if any, are silently ignored).\n";

  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// ASET ///
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_ASET (std::vector<std::string>& entry)
{
  START_TIMER("process_ASET")

  int node, dofs;
  for (size_t i = 0; i+1 < entry.size(); i+=2)
  {
    node = dofs = 0;
    CONVERT_ENTRY ("ASET",
		   fieldValue(entry[i],node) &&
		   fieldValue(entry[i+1],dofs));

    if (node > 0) nodeStat[node] = dofs;
  }

  STOPP_TIMER("process_ASET")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// ASET1 //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_ASET1 (std::vector<std::string>& entry)
{
  START_TIMER("process_ASET1")

  int node1 = 0, node2 = 0, dofs = 0;
  if (!entry.empty())
    CONVERT_ENTRY ("ASET1",fieldValue(entry.front(),dofs));

  for (size_t i = 1; i < entry.size(); i++)
    if (entry[i] == "THRU")
      node1 = node2;
    else
    {
      node2 = 0;
      CONVERT_ENTRY ("ASET1",fieldValue(entry[i],node2));
      if (node1 > 0)
      {
        for (int node = node1+1; node <= node2; node++)
          nodeStat[node] = dofs;
        node1 = 0;
      }
      else if (node2 > 0)
        nodeStat[node2] = dofs;
    }

  STOPP_TIMER("process_ASET1")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// BAROR //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_BAROR (std::vector<std::string>& entry)
{
  START_TIMER("process_BAROR")

  if (barDefault)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: More than one BAROR entries were encountered,"
           <<" only the first one is used (Line: "<< lineCounter <<").\n";
    STOPP_TIMER("process_BAROR")
    return true;
  }

  barDefault = new BEAMOR(true);
  barDefault->G0 = -999;

  if (entry.size() < 7) entry.resize(7,"");

  CONVERT_ENTRY ("BAROR",   entry[0].empty() &&
		 fieldValue(entry[1],barDefault->PID) &&
		            entry[2].empty() &&
		            entry[3].empty() &&
		 fieldValue(entry[4],barDefault->X[0]) &&
		 fieldValue(entry[5],barDefault->X[1]) &&
		 fieldValue(entry[6],barDefault->X[2]));

  if (!entry[4].empty() && entry[5].empty() && entry[6].empty())
    if (entry[4].find('.') == std::string::npos)
      barDefault->G0 = barDefault->X[0];

  barDefault->empty[0] = entry[1].empty();
  barDefault->empty[1] = entry[4].empty();
  barDefault->empty[2] = entry[5].empty();
  barDefault->empty[3] = entry[6].empty();
  barDefault->empty[4] = (barDefault->G0 > 0 ? false : true);

#ifdef FFL_DEBUG
  std::cout <<"Default bar-orientation, PID = "<< barDefault->PID;
  if (barDefault->G0 >= 0)
    std::cout <<": G0 = "<< barDefault->G0 << std::endl;
  else
    std::cout <<": X = "<< barDefault->X << std::endl;
#endif

  STOPP_TIMER("process_BAROR")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// BEAMOR /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_BEAMOR (std::vector<std::string>& entry)
{
  START_TIMER("process_BEAMOR")

  if (beamDefault)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: More than one BEAMOR entries were encountered,"
           <<" only the first one is used (Line: "<< lineCounter <<").\n";
    STOPP_TIMER("process_BEAMOR")
    return true;
  }

  beamDefault = new BEAMOR;
  beamDefault->G0 = -999;

  if (entry.size() < 7) entry.resize(7,"");

  CONVERT_ENTRY ("BEAMOR",  entry[0].empty() &&
		 fieldValue(entry[1],beamDefault->PID) &&
		            entry[2].empty() &&
		            entry[3].empty() &&
		 fieldValue(entry[4],beamDefault->X[0]) &&
		 fieldValue(entry[5],beamDefault->X[1]) &&
		 fieldValue(entry[6],beamDefault->X[2]));

  if (!entry[4].empty() && entry[5].empty() && entry[6].empty())
    if (entry[4].find('.') == std::string::npos)
      beamDefault->G0 = beamDefault->X[0];

  beamDefault->empty[0] = entry[1].empty();
  beamDefault->empty[1] = entry[4].empty();
  beamDefault->empty[2] = entry[5].empty();
  beamDefault->empty[3] = entry[6].empty();
  beamDefault->empty[4] = beamDefault->G0 > 0 ? false : true;

#ifdef FFL_DEBUG
  std::cout <<"Default beam-orientation, PID = "<< beamDefault->PID;
  if (beamDefault->G0 >= 0)
    std::cout <<": G0 = "<< beamDefault->G0 << std::endl;
  else
    std::cout <<": X = "<< beamDefault->X << std::endl;
#endif

  STOPP_TIMER("process_BEAMOR")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CBAR ///
////////////////////////////////////////////////////////////////////////////////

// Auxilliary function to order the digits in an integer in ascending order

static int sortDOFs (int C)
{
  bool on[6] = { false, false, false, false, false, false };
  int digit;
  while (C != 0)
  {
    digit = C%10;
    if (digit > 0 && digit < 7) on[digit-1] = true;
    C /= 10;
  }

  for (digit = 1; digit < 7; digit++)
    if (on[digit-1]) C = 10*C + digit;

  return C;
}

bool FFlNastranReader::process_CBAR (std::vector<std::string>& entry)
{
  START_TIMER("process_CBAR")

  int    EID, PID = 0, PA = 0, PB = 0;
  FaVec3 X, WA, WB;

  std::vector<int> G(2,0);

  if (entry.size() < 16) entry.resize(16,"");

  CONVERT_ENTRY ("CBAR",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],PID) &&
		 fieldValue(entry[2],G[0]) &&
		 fieldValue(entry[3],G[1]) &&
		 fieldValue(entry[4],X[0]) &&
		 fieldValue(entry[5],X[1]) &&
		 fieldValue(entry[6],X[2]) &&
	 	            entry[7].empty() &&
		 fieldValue(entry[8],PA) &&
		 fieldValue(entry[9],PB) &&
		 fieldValue(entry[10],WA[0]) &&
		 fieldValue(entry[11],WA[1]) &&
		 fieldValue(entry[12],WA[2]) &&
		 fieldValue(entry[13],WB[0]) &&
		 fieldValue(entry[14],WB[1]) &&
		 fieldValue(entry[15],WB[2]));

  if (entry.front().empty())
    EID = myLink->getNewElmID();

  // Store the property and beam orientation data temporarily in the bOri map
  // and must be resolved later after all coordinate systems have been read
  BEAMOR* bo = new BEAMOR(true);
  bo->PID = PID;
  if (!entry[4].empty() && entry[5].empty() && entry[6].empty())
    if (entry[4].find('.') == std::string::npos)
      bo->G0 = X[0];

  bo->X = X;
  bo->empty[0] = entry[1].empty();
  bo->empty[1] = entry[4].empty();
  bo->empty[2] = entry[5].empty();
  bo->empty[3] = entry[6].empty();
  bo->empty[4] = bo->G0 > 0 ? false : true;
  bOri[EID] = bo;

  FFlElementBase* theBeam = createElement("BEAM2",EID,G);

  if (PA != 0 || PB != 0)
  {
    // This beam element has pin-flags (local DOFs to be released)
    FFlPBEAMPIN* myAtt = CREATE_ATTRIBUTE(FFlPBEAMPIN,"PBEAMPIN",EID);
    myAtt->PA = sortDOFs(PA);
    myAtt->PB = sortDOFs(PB);
    theBeam->setAttribute("PBEAMPIN",myLink->addUniqueAttribute(myAtt));
  }

  if (WA.sqrLength() > 0.0 || WB.sqrLength() > 0.0)
  {
    // This beam element has eccentricities
    FFlPBEAMECCENT* myAtt = CREATE_ATTRIBUTE(FFlPBEAMECCENT,"PBEAMECCENT",EID);
    myAtt->node1Offset = WA;
    myAtt->node2Offset = WB;
#if FFL_DEBUG > 1
    myAtt->print();
#endif
    myLink->addAttribute(myAtt);
    theBeam->setAttribute("PBEAMECCENT",EID);
  }

  sizeOK = myLink->addElement(theBeam);

  STOPP_TIMER("process_CBAR")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CBEAM //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CBEAM (std::vector<std::string>& entry)
{
  START_TIMER("process_CBEAM")

  int    EID, PID = 0, PA = 0, PB = 0;
  FaVec3 X, WA, WB;

  std::vector<int> G(2,0);

  if (entry.size() < 16) entry.resize(16,"");

  CONVERT_ENTRY ("CBEAM",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],PID) &&
		 fieldValue(entry[2],G[0]) &&
		 fieldValue(entry[3],G[1]) &&
		 fieldValue(entry[4],X[0]) &&
		 fieldValue(entry[5],X[1]) &&
		 fieldValue(entry[6],X[2]) &&
		 fieldValue(entry[8],PA) &&
		 fieldValue(entry[9],PB) &&
		 fieldValue(entry[10],WA[0]) &&
		 fieldValue(entry[11],WA[1]) &&
		 fieldValue(entry[12],WA[2]) &&
		 fieldValue(entry[13],WB[0]) &&
		 fieldValue(entry[14],WB[1]) &&
		 fieldValue(entry[15],WB[2]));

  if (entry.front().empty())
    EID = myLink->getNewElmID();

  // Store the property and beam orientation data temporarily in the bOri map
  // and must be resolved later after all coordinate systems have been read
  BEAMOR* bo = new BEAMOR;
  bo->PID = PID;
  if (!entry[4].empty() && entry[5].empty() && entry[6].empty())
    if (entry[4].find('.') == std::string::npos)
      bo->G0 = X[0];

  if (entry[7] == "BGG")
    bo->basic = true; // Orientation vector in the Basic coordinate system
  else if (!entry[7].empty() && entry[7] != "GGG")
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: CBEAM element "<< EID
           <<" has offset vector flag \""<< entry[7] <<"\"."
           <<"\n              This is not implemented yet, \"GGG\" is used.\n";
  }

  bo->X = X;
  bo->empty[0] = entry[1].empty();
  bo->empty[1] = entry[4].empty();
  bo->empty[2] = entry[5].empty();
  bo->empty[3] = entry[6].empty();
  bo->empty[4] = bo->G0 > 0 ? false : true;
  bOri[EID] = bo;

  FFlElementBase* theBeam = createElement("BEAM2",EID,G);

  if (PA != 0 || PB != 0)
  {
    // This beam element has pin-flags (local DOFs to be released)
    FFlPBEAMPIN* myAtt = CREATE_ATTRIBUTE(FFlPBEAMPIN,"PBEAMPIN",EID);
    myAtt->PA = sortDOFs(PA);
    myAtt->PB = sortDOFs(PB);
    theBeam->setAttribute("PBEAMPIN",myLink->addUniqueAttribute(myAtt));
  }

  if (WA.sqrLength() > 0.0 || WB.sqrLength() > 0.0)
  {
    // This beam element has eccentricities
    FFlPBEAMECCENT* myAtt = CREATE_ATTRIBUTE(FFlPBEAMECCENT,"PBEAMECCENT",EID);
    myAtt->node1Offset = WA;
    myAtt->node2Offset = WB;
#if FFL_DEBUG > 1
    myAtt->print();
#endif
    myLink->addAttribute(myAtt);
    theBeam->setAttribute("PBEAMECCENT",EID);
  }

  sizeOK = myLink->addElement(theBeam);

  STOPP_TIMER("process_CBEAM")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CBUSH //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CBUSH (std::vector<std::string>& entry)
{
  START_TIMER("process_CBUSH")

  int    EID, PID = 0, OCID = -1, CID = -1;
  double S = 0.5;
  FaVec3 X, Sv;

  std::vector<int> G(2,0);

  if (entry.size() < 13) entry.resize(13,"");

  CONVERT_ENTRY ("CBUSH",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],PID) &&
		 fieldValue(entry[2],G[0]) &&
		 fieldValue(entry[3],G[1]) &&
		 fieldValue(entry[4],X[0]) &&
		 fieldValue(entry[5],X[1]) &&
		 fieldValue(entry[6],X[2]) &&
		 fieldValue(entry[7],CID) &&
		 fieldValue(entry[8],S) &&
		 fieldValue(entry[9],OCID) &&
		 fieldValue(entry[10],Sv[0]) &&
		 fieldValue(entry[11],Sv[1]) &&
		 fieldValue(entry[12],Sv[2]));

  if (entry.front().empty())
    EID = myLink->getNewElmID();

  if (G[0] == 0 || G[1] == 0 || G[0] == G[1])
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: CBUSH element "<< EID <<" is grounded."
           <<"\n              This is not supported so this element is ignored."
           <<"\n              Use a Free Joint on the system level instead.\n";

    STOPP_TIMER("process_CBUSH")
    return true;
  }

  FFlElementBase* theBush = createElement("BUSH",EID,G,PID);

  if (CID >= 0)
  {
    sprComp[EID] = CID; // Store the coordinate system ID
    if (CID > 0) theBush->setAttribute("PCOORDSYS",CID);
  }
  else if (!entry[4].empty() || !entry[5].empty() || !entry[6].empty())
  {
    // This bush element has an orientation vector
    FFlPORIENT* myOr = CREATE_ATTRIBUTE(FFlPORIENT,"PORIENT",EID);
    if (entry[5].empty() && entry[6].empty())
      if (entry[4].find('.') == std::string::npos)
        CID = X[0];
    if (CID > 0)
      sprComp[EID] = -CID; // Store node for later computation of orientation
    else
      myOr->directionVector = X*1.1; // Explicitly defined orientation vector
    myLink->addAttribute(myOr);
    theBush->setAttribute("PORIENT",EID);
  }
  else
    sprComp[EID] = -1; // Neither a local system nor an orientation vector

  if (OCID >= 0)
  {
    // Explicitly defined offset vector
    FFlPBUSHECCENT* myEcc = CREATE_ATTRIBUTE(FFlPBUSHECCENT,"PBUSHECCENT",EID);
    myEcc->offset = OCID == 0 ? Sv.round(10) : Sv;
#if FFL_DEBUG > 1
    myEcc->print();
#endif
    myLink->addAttribute(myEcc);
    theBush->setAttribute("PBUSHECCENT",EID);
    if (OCID > 0) sprPID[EID] = OCID; // Store coordinate system id
  }
  else if (S < 0.0 || S > 1.0)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Invalid location parameter S = "<< S
           <<" for CBUSH element "<< EID
           <<".\n              Reset to default value 0.5\n";
  }
  else if (S != 0.5)
    sprK[EID] = S; // Relative offset along line from node 1 to node 2

  sizeOK = myLink->addElement(theBush);

  STOPP_TIMER("process_CBUSH")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CELAS1 /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CELAS1 (std::vector<std::string>& entry)
{
  START_TIMER("process_CELAS1")

  int EID, PID, C1 = 0, C2 = 0;

  std::vector<int> G(2,0);

  if (entry.size() < 6) entry.resize(6,"");

  CONVERT_ENTRY ("CELAS1",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],PID) &&
		 fieldValue(entry[2],G[0]) &&
		 fieldValue(entry[3],C1) &&
		 fieldValue(entry[4],G[1]) &&
		 fieldValue(entry[5],C2));

  if (entry.front().empty())
    EID = myLink->getNewElmID();

  if (G[0] == 0 || G[1] == 0)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: CELAS1 element "<< EID <<" is grounded."
           <<"\n              This is not supported so this element is ignored."
           <<"\n              Use a Free Joint on the system level instead.\n";

    STOPP_TIMER("process_CELAS1")
    return true;
  }

  if (!entry[1].empty()) sprPID[EID] = PID;
  sprComp[EID] = C1;

  if (C1 != C2)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: CELAS1 element "<< EID <<" does not have the same"
           <<" component numbers at its two nodes: "<< C1 <<", "<< C2 <<".\n"
           <<"\n              The value for node 1 will be used for both.\n";
  }

  bool hasTranslation = false;
  bool hasRotation = false;
  for (int i = C1; i > 0; i /= 10)
    if (i%10 <= 3)
      hasTranslation = true;
    else if (i%10 <= 6)
      hasRotation = true;

  if (hasTranslation && hasRotation)
  {
    hasRotation = false;
    nWarnings++;
    ListUI <<"\n  ** Warning: CELAS1 element "<< EID
           <<" has both translational and rotational components: "<< C1
           <<"\n              The rotational components will be ignored."
           <<"\n              They have to be specified through a separate"
           <<" CELAS1 element, if needed.\n";
  }

  myLink->addAttribute(AttributeFactory::instance()->create("PSPRING",EID));

  if (hasRotation)
    sizeOK = myLink->addElement(createElement("RSPRING",EID,G,EID));
  else
    sizeOK = myLink->addElement(createElement("SPRING",EID,G,EID));

  STOPP_TIMER("process_CELAS1")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CELAS2 /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CELAS2 (std::vector<std::string>& entry)
{
  START_TIMER("process_CELAS2")

  int    EID, C1 = 0, C2 = 0;
  double K = 0.0;

  std::vector<int> G(2,0);

  if (entry.size() < 6) entry.resize(6,"");

  CONVERT_ENTRY ("CELAS2",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],K) &&
		 fieldValue(entry[2],G[0]) &&
		 fieldValue(entry[3],C1) &&
		 fieldValue(entry[4],G[1]) &&
		 fieldValue(entry[5],C2));

  if (entry.front().empty())
    EID = myLink->getNewElmID();

  if (G[0] == 0 || G[1] == 0)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: CELAS2 element "<< EID <<" is grounded."
           <<"\n              This is not supported so this element is ignored."
           <<"\n              Use a Free Joint on the system level instead.\n";

    STOPP_TIMER("process_CELAS2")
    return true;
  }

  sprComp[EID] = C1;
  sprK[EID] = K;

  if (C1 != C2)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: CELAS2 element "<< EID <<" does not have the same"
           <<" component numbers at its two nodes: "<< C1 <<", "<< C2 <<".\n"
           <<"\n              The value for node 1 will be used for both.\n";
  }

  bool hasTranslation = false;
  bool hasRotation = false;
  for (int i = C1; i > 0; i /= 10)
    if (i%10 <= 3)
      hasTranslation = true;
    else if (i%10 <= 6)
      hasRotation = true;

  if (hasTranslation && hasRotation)
  {
    hasRotation = false;
    nWarnings++;
    ListUI <<"\n  ** Warning: CELAS2 element "<< EID
           <<" has both translational and rotational components: "<< C1
           <<"\n              The rotational components will be ignored."
           <<"\n              They have to be specified through a separate"
           <<" CELAS2 element, if needed.\n";
  }

  myLink->addAttribute(AttributeFactory::instance()->create("PSPRING",EID));

  if (hasRotation)
    sizeOK = myLink->addElement(createElement("RSPRING",EID,G,EID));
  else
    sizeOK = myLink->addElement(createElement("SPRING",EID,G,EID));

  STOPP_TIMER("process_CELAS2")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CHEXA //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CHEXA (std::vector<std::string>& entry)
{
  START_TIMER("process_CHEXA")

  int EID, PID = 0;

  std::vector<int> G(20,0);

  if (entry.size() < 22) entry.resize(22,"");

  CONVERT_ENTRY ("CHEXA",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],PID) &&
		 fieldValue(entry[2],G[0]) &&
		 fieldValue(entry[3],G[1]) &&
		 fieldValue(entry[4],G[2]) &&
		 fieldValue(entry[5],G[3]) &&
		 fieldValue(entry[6],G[4]) &&
		 fieldValue(entry[7],G[5]) &&
		 fieldValue(entry[8],G[6]) &&
		 fieldValue(entry[9],G[7]) &&
		 fieldValue(entry[10],G[8]) &&
		 fieldValue(entry[11],G[9]) &&
		 fieldValue(entry[12],G[10]) &&
		 fieldValue(entry[13],G[11]) &&
		 fieldValue(entry[14],G[12]) &&
		 fieldValue(entry[15],G[13]) &&
		 fieldValue(entry[16],G[14]) &&
		 fieldValue(entry[17],G[15]) &&
		 fieldValue(entry[18],G[16]) &&
		 fieldValue(entry[19],G[17]) &&
		 fieldValue(entry[20],G[18]) &&
		 fieldValue(entry[21],G[19]) &&
		 G[0] && G[1] && G[2] && G[3] &&
		 G[4] && G[5] && G[6] && G[7]);

  if (entry.front().empty())
    EID = myLink->getNewElmID();

  if (G[ 8] && G[ 9] && G[10] && G[11] && G[12] && G[13] &&
      G[14] && G[15] && G[16] && G[17] && G[18] && G[19])
  {
    // 20-noded hexahedron     1  2  3  4  5  6  7  8  9 10
    const int nodeperm[20] = { 1, 3, 5, 7,13,15,17,19, 2, 4,
    //                        11 12 13 14 15 16 17 18 19 20
                               6, 8, 9,10,11,12,14,16,18,20 };
    std::vector<int> tmp(20,0);
    for (int i = 0; i < 20; i++) tmp[nodeperm[i]-1] = G[i];
    solidPID[EID] = PID;
    sizeOK = myLink->addElement(createElement("HEX20",EID,tmp));
  }
  else if (G[ 8] || G[ 9] || G[10] || G[11] || G[12] || G[13] ||
           G[14] || G[15] || G[16] || G[17] || G[18] || G[19])
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: CHEXA element "<< EID
           <<" with invalid number of nodes ignored.\n             ";
    for (int i = 0; i < 20; i++)
      if (G[i]) ListUI <<" G"<< i+1 <<"="<< G[i];
    ListUI <<"\n";
  }
  else
  {
    // 8-noded hexahedron
    G.resize(8);
    solidPID[EID] = PID;
    sizeOK = myLink->addElement(createElement("HEX8",EID,G));
  }

  STOPP_TIMER("process_CHEXA")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CONM1 //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CONM1 (std::vector<std::string>& entry)
{
  START_TIMER("process_CONM1")

  int EID, CID = 0;

  std::vector<int>    G(1,0);
  std::vector<double> M(21,0.0);

  if (entry.size() < 24) entry.resize(24,"");

  CONVERT_ENTRY ("CONM1",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],G[0]) &&
		 fieldValue(entry[2],CID) &&
		 fieldValue(entry[3],M[0]) &&
		 fieldValue(entry[4],M[1]) &&
		 fieldValue(entry[5],M[2]) &&
		 fieldValue(entry[6],M[3]) &&
		 fieldValue(entry[7],M[4]) &&
		 fieldValue(entry[8],M[5]) &&
		 fieldValue(entry[9],M[6]) &&
		 fieldValue(entry[10],M[7]) &&
		 fieldValue(entry[11],M[8]) &&
		 fieldValue(entry[12],M[9]) &&
		 fieldValue(entry[13],M[10]) &&
		 fieldValue(entry[14],M[11]) &&
		 fieldValue(entry[15],M[12]) &&
		 fieldValue(entry[16],M[13]) &&
		 fieldValue(entry[17],M[14]) &&
		 fieldValue(entry[18],M[15]) &&
		 fieldValue(entry[19],M[16]) &&
		 fieldValue(entry[20],M[17]) &&
		 fieldValue(entry[21],M[18]) &&
		 fieldValue(entry[22],M[19]) &&
		 fieldValue(entry[23],M[20]));

  if (entry.front().empty())
    EID = myLink->getNewElmID();
  if (CID > 0)
    massCID[EID] = CID;

  FFlPMASS* myAtt = CREATE_ATTRIBUTE(FFlPMASS,"PMASS",EID);
  std::vector<double>& Mvec = myAtt->M.data();

  int nField = 21;
  while (nField > 0 && M[nField-1] == 0.0) nField--;
  Mvec.reserve(nField);
  for (int i = 0; i < nField; i++)
    Mvec.push_back(round(M[i],10));

#ifdef FFL_DEBUG
  myAtt->print();
#endif
  myLink->addAttribute(myAtt);
  sizeOK = myLink->addElement(createElement("CMASS",EID,G,EID));

  STOPP_TIMER("process_CONM1")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CONM2 //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CONM2 (std::vector<std::string>& entry)
{
  START_TIMER("process_CONM2")

  int    EID, CID = 0;
  double M = 0.0;
  FaVec3 X;

  std::vector<int>    G(1,0);
  std::vector<double> I(6,0.0);

  if (entry.size() < 14) entry.resize(14,"");

  CONVERT_ENTRY ("CONM2",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],G[0]) &&
		 fieldValue(entry[2],CID) &&
		 fieldValue(entry[3],M) &&
		 fieldValue(entry[4],X[0]) &&
		 fieldValue(entry[5],X[1]) &&
		 fieldValue(entry[6],X[2]) &&
		            entry[7].empty() &&
		 fieldValue(entry[8],I[0]) &&
		 fieldValue(entry[9],I[1]) &&
		 fieldValue(entry[10],I[2]) &&
		 fieldValue(entry[11],I[3]) &&
		 fieldValue(entry[12],I[4]) &&
		 fieldValue(entry[13],I[5]));

  if (entry.front().empty())
    EID = myLink->getNewElmID();
  if (CID > 0 || CID == -1)
    massCID[EID] = CID;

  FFlPMASS* myAtt = CREATE_ATTRIBUTE(FFlPMASS,"PMASS",EID);
  std::vector<double>& Mvec = myAtt->M.data();

  // Store mass value on the diagonal
  Mvec.reserve(6);
  Mvec.push_back(M);
  Mvec.push_back(0.0);
  Mvec.push_back(M);
  Mvec.push_back(0.0);
  Mvec.push_back(0.0);
  Mvec.push_back(M);

  for (int i = 0; i < 6; i++)
    if (I[i] != 0.0)
    {
      // Store the mass moments of inertia
      Mvec.resize(21,0.0);
      Mvec[9]  =  I[0];
      Mvec[13] = -I[1];
      Mvec[14] =  I[2];
      Mvec[18] = -I[3];
      Mvec[19] = -I[4];
      Mvec[20] =  I[5];
      break;
    }

  if (X.sqrLength() > 0.0 || CID == -1)
    massX[EID] = new FaVec3(X);
  else if (CID == 0)
    for (double& m : Mvec)
      m = round(m,10);

#ifdef FFL_DEBUG
  myAtt->print();
#endif
  myLink->addAttribute(myAtt);
  sizeOK = myLink->addElement(createElement("CMASS",EID,G,EID));

  STOPP_TIMER("process_CONM2")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CONROD /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CONROD (std::vector<std::string>& entry)
{
  START_TIMER("process_CONROD")

  int             EID, MID = 0;
  FFlCrossSection params;

  std::vector<int> G(2,0);

  if (entry.size() < 8) entry.resize(8,"");

  CONVERT_ENTRY ("CONROD",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],G[0]) &&
		 fieldValue(entry[2],G[1]) &&
		 fieldValue(entry[3],MID) &&
		 fieldValue(entry[4],params.A) &&
		 fieldValue(entry[5],params.J) &&
		 fieldValue(entry[7],params.NSM));

  if (entry.front().empty())
    EID = myLink->getNewElmID();

#ifdef FFL_DEBUG
  std::cout <<"Rod property, ID = "<< EID <<" --> material ID = "<< MID
	    << std::endl;
#endif

  this->insertBeamPropMat("PROD",EID,MID);
  myLink->addAttribute(createBeamSection(EID,params,lastComment));

  if (params.NSM != 0.0)
  {
    beamPIDnsm.insert(EID);
    myLink->addAttribute(createNSM(EID,params.NSM));
  }

  sizeOK = myLink->addElement(createElement("BEAM2",EID,G,EID));

  STOPP_TIMER("process_CONROD")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CORD1C /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CORD1C (std::vector<std::string>& entry)
{
  START_TIMER("process_CORD1C")

  int CID, G1, G2, G3;

  if (entry.size() < 4)
    entry.resize(4,"");
  else if (entry.size() > 4 && entry.size() < 8)
    entry.resize(8,"");

  for (size_t i = 0; i+3 < entry.size(); i++)
  {
    CONVERT_ENTRY ("CORD1C",
		   fieldValue(entry[i+0],CID) &&
		   fieldValue(entry[i+1],G1) &&
		   fieldValue(entry[i+2],G2) &&
		   fieldValue(entry[i+3],G3));

#ifdef FFL_DEBUG
    std::cout <<"Cylindrical coordinate system, CID = "<< CID
	      <<": G1 = "<< G1 <<" G2 = "<< G2 <<" G3 = "<< G3 << std::endl;
#endif

    CORD* theCord = new CORD;
    theCord->type = cylindrical;
    theCord->G[0] = G1;
    theCord->G[1] = G2;
    theCord->G[2] = G3;
    cordSys[CID] = theCord;
  }

  STOPP_TIMER("process_CORD1C")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CORD1R /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CORD1R (std::vector<std::string>& entry)
{
  START_TIMER("process_CORD1R")

  int CID, G1, G2, G3;

  if (entry.size() < 4)
    entry.resize(4,"");
  else if (entry.size() > 4 && entry.size() < 8)
    entry.resize(8,"");

  for (size_t i = 0; i+3 < entry.size(); i++)
  {
    CONVERT_ENTRY ("CORD1R",
		   fieldValue(entry[i+0],CID) &&
		   fieldValue(entry[i+1],G1) &&
		   fieldValue(entry[i+2],G2) &&
		   fieldValue(entry[i+3],G3));

#ifdef FFL_DEBUG
    std::cout <<"Rectangular coordinate system, CID = "<< CID
	      <<": G1 = "<< G1 <<" G2 = "<< G2 <<" G3 = "<< G3 << std::endl;
#endif

    CORD* theCord = new CORD;
    theCord->type = rectangular;
    theCord->G[0] = G1;
    theCord->G[1] = G2;
    theCord->G[2] = G3;
    cordSys[CID] = theCord;
  }

  STOPP_TIMER("process_CORD1R")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CORD1S /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CORD1S (std::vector<std::string>& entry)
{
  START_TIMER("process_CORD1S")

  int CID, G1, G2, G3;

  if (entry.size() < 4)
    entry.resize(4,"");
  else if (entry.size() > 4 && entry.size() < 8)
    entry.resize(8,"");

  for (size_t i = 0; i+3 < entry.size(); i++)
  {
    CONVERT_ENTRY ("CORD1S",
		   fieldValue(entry[i+0],CID) &&
		   fieldValue(entry[i+1],G1) &&
		   fieldValue(entry[i+2],G2) &&
		   fieldValue(entry[i+3],G3));

#ifdef FFL_DEBUG
    std::cout <<"Spherical coordinate system, CID = "<< CID
	      <<": G1 = "<< G1 <<" G2 = "<< G2 <<" G3 = "<< G3 << std::endl;
#endif

    CORD* theCord = new CORD;
    theCord->type = spherical;
    theCord->G[0] = G1;
    theCord->G[1] = G2;
    theCord->G[2] = G3;
    cordSys[CID] = theCord;
  }

  STOPP_TIMER("process_CORD1S")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CORD2C /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CORD2C (std::vector<std::string>& entry)
{
  START_TIMER("process_CORD2C")

  int    CID, RID = 0;
  FaVec3 A, B, C;

  if (entry.size() < 11) entry.resize(11,"");

  CONVERT_ENTRY ("CORD2C",
		 fieldValue(entry[0],CID) &&
		 fieldValue(entry[1],RID) &&
		 fieldValue(entry[2],A[0]) &&
		 fieldValue(entry[3],A[1]) &&
		 fieldValue(entry[4],A[2]) &&
		 fieldValue(entry[5],B[0]) &&
		 fieldValue(entry[6],B[1]) &&
		 fieldValue(entry[7],B[2]) &&
		 fieldValue(entry[8],C[0]) &&
		 fieldValue(entry[9],C[1]) &&
		 fieldValue(entry[10],C[2]));

#ifdef FFL_DEBUG
  std::cout <<"Cylindrical coordinate system, CID = "<< CID <<", RID = "<< RID
	    <<"\n  Origo = "<< A <<"\n  Zaxis = "<< B <<"\n  XZpnt = "<< C
	    << std::endl;
#endif

  CORD* theCord = new CORD;
  theCord->type = cylindrical;
  theCord->RID  = RID;
  theCord->G[0] = -999;
  theCord->Origo = A;
  theCord->Zaxis = B;
  theCord->XZpnt = C;
  cordSys[CID] = theCord;

  STOPP_TIMER("process_CORD2C")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CORD2R /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CORD2R (std::vector<std::string>& entry)
{
  START_TIMER("process_CORD2R")

  int CID, RID = 0;
  FaVec3 A, B, C;

  if (entry.size() < 11) entry.resize(11,"");

  CONVERT_ENTRY ("CORD2R",
		 fieldValue(entry[0],CID) &&
		 fieldValue(entry[1],RID) &&
		 fieldValue(entry[2],A[0]) &&
		 fieldValue(entry[3],A[1]) &&
		 fieldValue(entry[4],A[2]) &&
		 fieldValue(entry[5],B[0]) &&
		 fieldValue(entry[6],B[1]) &&
		 fieldValue(entry[7],B[2]) &&
		 fieldValue(entry[8],C[0]) &&
		 fieldValue(entry[9],C[1]) &&
		 fieldValue(entry[10],C[2]));

#ifdef FFL_DEBUG
  std::cout <<"Rectangular coordinate system, CID = "<< CID <<", RID = "<< RID
	    <<"\n  Origo = "<< A <<"\n  Zaxis = "<< B <<"\n  XZpnt = "<< C
	    << std::endl;
#endif

  CORD* theCord = new CORD;
  theCord->type = rectangular;
  theCord->RID  = RID;
  theCord->G[0] = -999;
  theCord->Origo = A;
  theCord->Zaxis = B;
  theCord->XZpnt = C;
  cordSys[CID] = theCord;

  STOPP_TIMER("process_CORD2R")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CORD2S /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CORD2S (std::vector<std::string>& entry)
{
  START_TIMER("process_CORD2S")

  int CID, RID = 0;
  FaVec3 A, B, C;

  if (entry.size() < 11) entry.resize(11,"");

  CONVERT_ENTRY ("CORD2S",
		 fieldValue(entry[0],CID) &&
		 fieldValue(entry[1],RID) &&
		 fieldValue(entry[2],A[0]) &&
		 fieldValue(entry[3],A[1]) &&
		 fieldValue(entry[4],A[2]) &&
		 fieldValue(entry[5],B[0]) &&
		 fieldValue(entry[6],B[1]) &&
		 fieldValue(entry[7],B[2]) &&
		 fieldValue(entry[8],C[0]) &&
		 fieldValue(entry[9],C[1]) &&
		 fieldValue(entry[10],C[2]));

#ifdef FFL_DEBUG
  std::cout <<"Spherical coordinate system, CID = "<< CID <<", RID = "<< RID
	    <<"\n  Origo = "<< A <<"\n  Zaxis = "<< B <<"\n  XZpnt = "<< C
	    << std::endl;
#endif

  CORD* theCord = new CORD;
  theCord->type = spherical;
  theCord->RID  = RID;
  theCord->G[0] = -999;
  theCord->Origo = A;
  theCord->Zaxis = B;
  theCord->XZpnt = C;
  cordSys[CID] = theCord;

  STOPP_TIMER("process_CORD2S")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CPENTA /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CPENTA (std::vector<std::string>& entry)
{
  START_TIMER("process_CPENTA")

  int EID, PID = 0;

  std::vector<int> G(15,0);

  if (entry.size() < 17) entry.resize(17,"");

  CONVERT_ENTRY ("CPENTA",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],PID) &&
		 fieldValue(entry[2],G[0]) &&
		 fieldValue(entry[3],G[1]) &&
		 fieldValue(entry[4],G[2]) &&
		 fieldValue(entry[5],G[3]) &&
		 fieldValue(entry[6],G[4]) &&
		 fieldValue(entry[7],G[5]) &&
		 fieldValue(entry[8],G[6]) &&
		 fieldValue(entry[9],G[7]) &&
		 fieldValue(entry[10],G[8]) &&
		 fieldValue(entry[11],G[9]) &&
		 fieldValue(entry[12],G[10]) &&
		 fieldValue(entry[13],G[11]) &&
		 fieldValue(entry[14],G[12]) &&
		 fieldValue(entry[15],G[13]) &&
		 fieldValue(entry[16],G[14]) &&
		 G[0] && G[1] && G[2] && G[3] && G[4] && G[5]);

  if (entry.front().empty())
    EID = myLink->getNewElmID();

  if (G[ 6] && G[ 7] && G[ 8] && G[ 9] && G[10] &&
      G[11] && G[12] && G[13] && G[14])
  {
    // 15-noded pentahedron    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
    const int nodeperm[15] = { 1, 3, 5,10,12,14, 2, 4, 6, 7, 8, 9,11,13,15 };

    std::vector<int> tmp(15,0);
    for (int i = 0; i < 15; i++) tmp[nodeperm[i]-1] = G[i];
    solidPID[EID] = PID;
    sizeOK = myLink->addElement(createElement("WEDG15",EID,tmp));
  }
  else if (G[ 6] || G[ 7] || G[ 8] || G[ 9] || G[10] ||
           G[11] || G[12] || G[13] || G[14])
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: CPENTA element "<< EID
           <<" with invalid number of nodes ignored.\n             ";
    for (int i = 0; i < 15; i++)
      if (G[i]) ListUI <<" G"<< i+1 <<"="<< G[i];
    ListUI <<"\n";
  }
  else
  {
    // 6-noded pentahedron
    G.resize(6);
    solidPID[EID] = PID;
    sizeOK = myLink->addElement(createElement("WEDG6",EID,G));
  }

  STOPP_TIMER("process_CPENTA")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CQUAD4 /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CQUAD4 (std::vector<std::string>& entry)
{
  START_TIMER("process_CQUAD4")

  int    EID, PID, MCID = 0;
  double THETA, ZOFFS;

  std::vector<int>    G(4,0);
  std::vector<double> T(4,0.0);

  if (entry.size() < 14) entry.resize(14,"");

  CONVERT_ENTRY ("CQUAD4",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],PID) &&
		 fieldValue(entry[2],G[0]) &&
		 fieldValue(entry[3],G[1]) &&
		 fieldValue(entry[4],G[2]) &&
		 fieldValue(entry[5],G[3]) &&
		 fieldValue2(entry[6],MCID,THETA) && // either MCID or THETA
		 fieldValue(entry[7],ZOFFS) &&
	                    entry[8].empty() &&
	                    entry[9].empty() &&
		 fieldValue(entry[10],T[0]) &&
		 fieldValue(entry[11],T[1]) &&
		 fieldValue(entry[12],T[2]) &&
		 fieldValue(entry[13],T[3]));

  if (entry.front().empty())
    EID = myLink->getNewElmID();
  if (entry[1].empty())
    PID = EID;
  else
    shellPID[EID] = PID;

  double thickness = (T[0]+T[1]+T[2]+T[3])/4.0;
  if (thickness > 0.0)
  {
    PID = EID;
    FFlPTHICK* myAtt = CREATE_ATTRIBUTE(FFlPTHICK,"PTHICK",PID);
    myAtt->thickness = round(thickness,10);
#if FFL_DEBUG > 1
    myAtt->print();
#endif
    myLink->addAttribute(myAtt);
    PTHICKs.insert(PID);
  }
  else
    PID = 0;

  sizeOK = myLink->addElement(createElement("QUAD4",EID,G,PID,MCID));

  STOPP_TIMER("process_CQUAD4")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CQUAD8 /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CQUAD8 (std::vector<std::string>& entry)
{
  if (FFlReaders::convertToLinear == 1)
  {
    static bool firstCall = true;
    if (firstCall)
    {
      firstCall = false;
      nNotes++;
      ListUI <<"\n   * Note: Bulk input contains CQUAD8 shell elements.\n"
             <<"           These are all converted to CQUAD4 elements.\n";
    }

    if (entry.size() < 14) entry.resize(14,"");

    if (entry.size() < 15)
      entry[6].erase();
    else
      entry[6] = entry[14];

    if (entry.size() < 16)
      entry[7].erase();
    else
      entry[7] = entry[15];

    entry[8].erase();
    entry[9].erase();

    return process_CQUAD4(entry);
  }

  START_TIMER("process_CQUAD8")

  int    EID, PID, MCID = 0;
  double THETA, ZOFFS;

  std::vector<int>    G(8,0);
  std::vector<double> T(4,0.0);

  if (entry.size() < 17) entry.resize(17,"");

  CONVERT_ENTRY ("CQUAD8",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],PID) &&
		 fieldValue(entry[2],G[0]) &&
		 fieldValue(entry[3],G[2]) &&
		 fieldValue(entry[4],G[4]) &&
		 fieldValue(entry[5],G[6]) &&	
		 fieldValue(entry[6],G[1]) &&
		 fieldValue(entry[7],G[3]) &&
		 fieldValue(entry[8],G[5]) &&
		 fieldValue(entry[9],G[7]) &&
		 fieldValue(entry[10],T[0]) &&
		 fieldValue(entry[11],T[1]) &&
		 fieldValue(entry[12],T[2]) &&
		 fieldValue(entry[13],T[3]) &&
		 fieldValue2(entry[14],MCID,THETA) && // either MCID or THETA
		 fieldValue(entry[15],ZOFFS) &&
		 entry[16].empty());

  if (entry.front().empty())
    EID = myLink->getNewElmID();
  if (entry[1].empty())
    PID = EID;
  else
    shellPID[EID] = PID;

  double thickness = (T[0]+T[1]+T[2]+T[3])/4.0;
  if (thickness > 0.0)
  {
    PID = EID;
    FFlPTHICK* myAtt = CREATE_ATTRIBUTE(FFlPTHICK,"PTHICK",PID);
    myAtt->thickness = round(thickness,10);
#if FFL_DEBUG > 1
    myAtt->print();
#endif
    myLink->addAttribute(myAtt);
    PTHICKs.insert(PID);
  }
  else
    PID = 0;

  sizeOK = myLink->addElement(createElement("QUAD8",EID,G,PID,MCID));

  STOPP_TIMER("process_CQUAD8")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CROD ///
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CROD (std::vector<std::string>& entry)
{
  START_TIMER("process_CROD")

  int EID, PID = 0;

  std::vector<int> G(2,0);

  if (entry.size() < 4) entry.resize(4,"");

  CONVERT_ENTRY ("CROD",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],PID) &&
		 fieldValue(entry[2],G[0]) &&
		 fieldValue(entry[3],G[1]));

  if (entry.front().empty())
    EID = myLink->getNewElmID();

  sizeOK = myLink->addElement(createElement("BEAM2",EID,G,PID));

  STOPP_TIMER("process_CROD")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CTETRA /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CTETRA (std::vector<std::string>& entry)
{
  START_TIMER("process_CTETRA")

  int EID, PID = 0;

  std::vector<int> G(10,0);

  if (entry.size() < 12) entry.resize(12,"");

  CONVERT_ENTRY ("CTETRA",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],PID) &&
		 fieldValue(entry[2],G[0]) &&
		 fieldValue(entry[3],G[1]) &&
		 fieldValue(entry[4],G[2]) &&
		 fieldValue(entry[5],G[3]) &&
		 fieldValue(entry[6],G[4]) &&
		 fieldValue(entry[7],G[5]) &&
		 fieldValue(entry[8],G[6]) &&
		 fieldValue(entry[9],G[7]) &&
		 fieldValue(entry[10],G[8]) &&
		 fieldValue(entry[11],G[9]) &&
		 G[0] && G[1] && G[2] && G[3]);

  if (entry.front().empty())
    EID = myLink->getNewElmID();

  if (G[4] && G[5] && G[6] && G[7] && G[8] && G[9])
  {
    // 10-noded tetrahedron    1  2  3  4  5  6  7  8  9 10
    const int nodeperm[10] = { 1, 3, 5,10, 2, 4, 6, 7, 8, 9 };

    std::vector<int> tmp(10,0);
    for (int i = 0; i < 10; i++) tmp[nodeperm[i]-1] = G[i];
    solidPID[EID] = PID;
    sizeOK = myLink->addElement(createElement("TET10",EID,tmp));
  }
  else if (G[4] || G[5] || G[6] || G[7] || G[8] || G[9])
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: CTETRA element "<< EID
           <<" with invalid number of nodes ignored.\n             ";
    for (int i = 0; i < 10; i++)
      if (G[i]) ListUI <<" G"<< i+1 <<"="<< G[i];
    ListUI <<"\n";
  }
  else
  {
    // 4-noded tetrahedron
    G.resize(4);
    solidPID[EID] = PID;
    sizeOK = myLink->addElement(createElement("TET4",EID,G));
  }

  STOPP_TIMER("process_CTETRA")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CTRIA3 /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CTRIA3 (std::vector<std::string>& entry)
{
  START_TIMER("process_CTRIA3")

  int    EID, PID, MCID = 0;
  double THETA, ZOFFS;

  std::vector<int>    G(3,0);
  std::vector<double> T(3,0.0);

  if (entry.size() < 13) entry.resize(13,"");

  CONVERT_ENTRY ("CTRIA3",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],PID) &&
		 fieldValue(entry[2],G[0]) &&
		 fieldValue(entry[3],G[1]) &&
		 fieldValue(entry[4],G[2]) &&
		 fieldValue2(entry[5],MCID,THETA) && // either MCID or THETA
		 fieldValue(entry[6],ZOFFS) &&
		            entry[7].empty() &&
		            entry[8].empty() &&
		            entry[9].empty() &&
		 fieldValue(entry[10],T[0]) &&
		 fieldValue(entry[11],T[1]) &&
		 fieldValue(entry[12],T[2]));

  if (entry.front().empty())
    EID = myLink->getNewElmID();
  if (entry[1].empty())
    PID = EID;
  else
    shellPID[EID] = PID;

  double thickness = (T[0]+T[1]+T[2])/3.0;
  if (thickness > 0.0)
  {
    PID = EID;
    FFlPTHICK* myAtt = CREATE_ATTRIBUTE(FFlPTHICK,"PTHICK",PID);
    myAtt->thickness = round(thickness,10);
#if FFL_DEBUG > 1
    myAtt->print();
#endif
    myLink->addAttribute(myAtt);
    PTHICKs.insert(PID);
  }
  else
    PID = 0;

  sizeOK = myLink->addElement(createElement("TRI3",EID,G,PID,MCID));

  STOPP_TIMER("process_CTRIA3")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CTRIA6 /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CTRIA6 (std::vector<std::string>& entry)
{
  if (FFlReaders::convertToLinear == 1)
  {
    static bool firstCall = true;
    if (firstCall)
    {
      firstCall = false;
      nNotes++;
      ListUI <<"\n   * Note: Bulk input contains CTRIA6 shell elements.\n"
             <<"           These are all converted to CTRIA3 elements.\n";
    }

    if (entry.size() < 13) entry.resize(13,"");

    entry[5] = entry[8];
    entry[6] = entry[9];
    entry[7].erase();
    entry[8].erase();
    entry[9].erase();

    return process_CTRIA3(entry);
  }

  START_TIMER("process_CTRIA6")

  int    EID, PID, MCID = 0;
  double THETA, ZOFFS;

  std::vector<int>    G(6,0);
  std::vector<double> T(3,0.0);

  if (entry.size() < 14) entry.resize(14,"");

  CONVERT_ENTRY ("CTRIA6",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],PID) &&
		 fieldValue(entry[2],G[0]) &&
		 fieldValue(entry[3],G[2]) &&
		 fieldValue(entry[4],G[4]) &&
		 fieldValue(entry[5],G[1]) &&
		 fieldValue(entry[6],G[3]) &&
		 fieldValue(entry[7],G[5]) &&
		 fieldValue2(entry[8],MCID,THETA) && // either MCID or THETA
		 fieldValue(entry[9],ZOFFS) &&
		 fieldValue(entry[10],T[0]) &&
		 fieldValue(entry[11],T[1]) &&
		 fieldValue(entry[12],T[2]) &&
		            entry[13].empty());

  if (entry.front().empty())
    EID = myLink->getNewElmID();
  if (entry[1].empty())
    PID = EID;
  else
    shellPID[EID] = PID;

  double thickness = (T[0]+T[1]+T[2])/3.0;
  if (thickness > 0.0)
  {
    PID = EID;
    FFlPTHICK* myAtt = CREATE_ATTRIBUTE(FFlPTHICK,"PTHICK",PID);
    myAtt->thickness = round(thickness,10);
#if FFL_DEBUG > 1
    myAtt->print();
#endif
    myLink->addAttribute(myAtt);
    PTHICKs.insert(PID);
  }
  else
    PID = 0;

  sizeOK = myLink->addElement(createElement("TRI6",EID,G,PID,MCID));

  STOPP_TIMER("process_CTRIA6")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// CWELD //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_CWELD (std::vector<std::string>& entry)
{
  START_TIMER("process_CWELD")

  int EWID, PWID = 0, GS = 0;

  std::vector<int> G(2,0);

  if (entry.size() < 4) entry.resize(4,"");

  CONVERT_ENTRY ("CWELD",
		 fieldValue(entry[0],EWID) &&
		 fieldValue(entry[1],PWID) &&
		 fieldValue(entry[2],GS));

  if (entry.front().empty())
    EWID = myLink->getNewElmID();

  if (entry[3] == "GRIDID")
  {
    std::vector<int> GA(1,0), GB(1,0);
    if (entry.size() < 24) entry.resize(24,"");

    if (entry[6] == "QQ" || entry[6] == "QT" || entry[6] == "Q")
      if (entry[12].empty())
      {
	// 4-noded patch A
	GA.resize(5,0);
	CONVERT_ENTRY ("CWELD",
		       fieldValue(entry[8] ,GA[1]) &&
		       fieldValue(entry[9] ,GA[2]) &&
		       fieldValue(entry[10],GA[3]) &&
		       fieldValue(entry[11],GA[4]));
      }
      else
      {
	// 8-noded patch A
	GA.resize(9,0);
	CONVERT_ENTRY ("CWELD",
		       fieldValue(entry[8] ,GA[1]) &&
		       fieldValue(entry[9] ,GA[2]) &&
		       fieldValue(entry[10],GA[3]) &&
		       fieldValue(entry[11],GA[4]) &&
		       fieldValue(entry[12],GA[5]) &&
		       fieldValue(entry[13],GA[6]) &&
		       fieldValue(entry[14],GA[7]) &&
		       fieldValue(entry[15],GA[8]));
      }

    else if (entry[6] == "TQ" || entry[6] == "TT" || entry[6] == "T")
      if (entry[11].empty())
      {
	// 3-noded patch A
	GA.resize(4,0);
	CONVERT_ENTRY ("CWELD",
		       fieldValue(entry[8] ,GA[1]) &&
		       fieldValue(entry[9] ,GA[2]) &&
		       fieldValue(entry[10],GA[3]));
      }
      else
      {
	// 6-noded patch A
	GA.resize(7,0);
	CONVERT_ENTRY ("CWELD",
		       fieldValue(entry[8] ,GA[1]) &&
		       fieldValue(entry[9] ,GA[2]) &&
		       fieldValue(entry[10],GA[3]) &&
		       fieldValue(entry[11],GA[4]) &&
		       fieldValue(entry[12],GA[5]) &&
		       fieldValue(entry[13],GA[6]));
      }

    else
    {
      nWarnings++;
      ListUI <<"\n  ** Warning: Invalid SPTYP code \""<< entry[6]
             <<"\" for CWELD element "<< EWID
             <<".\n              This element is ignored.\n";
      STOPP_TIMER("process_CWELD")
      return true;
    }

    if (entry[6] == "QQ" || entry[6] == "TQ")
    {
      if (entry[20].empty())
      {
	// 4-noded patch B
	GB.resize(5,0);
	CONVERT_ENTRY ("CWELD",
		       fieldValue(entry[16],GB[1]) &&
		       fieldValue(entry[17],GB[2]) &&
		       fieldValue(entry[18],GB[3]) &&
		       fieldValue(entry[19],GB[4]));
      }
      else
      {
	// 8-noded patch B
	GB.resize(9,0);
	CONVERT_ENTRY ("CWELD",
		       fieldValue(entry[16],GB[1]) &&
		       fieldValue(entry[17],GB[2]) &&
		       fieldValue(entry[18],GB[3]) &&
		       fieldValue(entry[19],GB[4]) &&
		       fieldValue(entry[20],GB[5]) &&
		       fieldValue(entry[21],GB[6]) &&
		       fieldValue(entry[22],GB[7]) &&
		       fieldValue(entry[23],GB[8]));
      }
    }
    else if (entry[6] == "QT" || entry[6] == "TT")
    {
      if (entry[19].empty())
      {
	// 3-noded patch B
	GB.resize(4,0);
	CONVERT_ENTRY ("CWELD",
		       fieldValue(entry[16],GB[1]) &&
		       fieldValue(entry[17],GB[2]) &&
		       fieldValue(entry[18],GB[3]));
      }
      else
      {
	// 6-noded patch B
	GB.resize(7,0);
	CONVERT_ENTRY ("CWELD",
		       fieldValue(entry[16],GB[1]) &&
		       fieldValue(entry[17],GB[2]) &&
		       fieldValue(entry[18],GB[3]) &&
		       fieldValue(entry[19],GB[4]) &&
		       fieldValue(entry[20],GB[5]) &&
		       fieldValue(entry[21],GB[6]));
      }
    }

    weld.resize(2);

    // Add a property-less WAVGM element for each surface patch.
    // The element IDs are temporarily set to zero since we don't
    // know which number will be used by other elements yet.
    FFlElementBase* pch = createElement("WAVGM",0,GA);
    myLink->addElement(pch);
    weld[0][EWID] = pch;
    if (GB.size() > 1)
    {
      pch = createElement("WAVGM",0,GB);
      myLink->addElement(pch);
      weld[1][EWID] = pch;
    }
    else if (GS > 0)
      G[1] = GS; // This is a point-to-patch connection
  }

  else if (entry[3] == "ELEMID")
  {
    if (entry.size() < 10) entry.resize(10,"");

    CONVERT_ENTRY ("CWELD",
		   fieldValue(entry[8],G[0]) &&
		   fieldValue(entry[9],G[1]));

    // Negate the shell element ID numbers to indicate that these must be
    // resolved into surface patches later, onto which the point GS will be
    // projected to obtain the actual nodal points of the weld connector.
    G[0] = -G[0];
    G[1] = -G[1];
  }

  else if (entry[3] == "ALIGN")
  {
    if (entry.size() < 6) entry.resize(6,"");

    GS = 0;
    CONVERT_ENTRY ("CWELD",
		   fieldValue(entry[4],G[0]) &&
		   fieldValue(entry[5],G[1]));
  }

  else
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Invalid connection type \""<< entry[3]
           <<"\" for CWELD element "<< EWID
           <<".\n              This element is ignored.\n";
    STOPP_TIMER("process_CWELD")
    return true;
  }

  if (GS > 0) weldGS[EWID] = GS;
  FFlElementBase* theWeld = createElement("BEAM2",EWID,G,PWID);
  myWelds.push_back(theWeld);
  sizeOK = myLink->addElement(theWeld);

  STOPP_TIMER("process_CWELD")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// FORCE //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_FORCE (std::vector<std::string>& entry)
{
  START_TIMER("process_FORCE")

  int    SID = 0, G = 0, CID = 0;
  double F = 0.0;
  FaVec3 N;

  if (entry.size() < 7) entry.resize(7,"");

  CONVERT_ENTRY ("FORCE",
		 fieldValue(entry[0],SID) &&
		 fieldValue(entry[1],G) &&
		 fieldValue(entry[2],CID) &&
		 fieldValue(entry[3],F) &&
		 fieldValue(entry[4],N[0]) &&
		 fieldValue(entry[5],N[1]) &&
		 fieldValue(entry[6],N[2]));

  N *= F;
#ifdef FFL_DEBUG
  std::cout <<"Concentrated force, SID = "<< SID <<" --> F = "<< N
	    <<", G = "<< G << std::endl;
#endif

  FFlLoadBase* theLoad = LoadFactory::instance()->create("CFORCE",SID);
  theLoad->setValue(N);
  theLoad->setTarget(G);
  myLink->addLoad(theLoad);

  if (CID > 0) loadCID[theLoad] = CID;

  STOPP_TIMER("process_FORCE")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// GRDSET /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_GRDSET (std::vector<std::string>& entry)
{
  START_TIMER("process_GRDSET")

  if (gridDefault)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: More than one GRDSET entries were encountered,"
           <<" only the first one is used (Line: "<< lineCounter <<").\n";
    STOPP_TIMER("process_GRDSET")
    return true;
  }

  gridDefault = new GRDSET;

  if (entry.size() < 8) entry.resize(8,"");

  CONVERT_ENTRY ("GRDSET",  entry[0].empty() &&
		 fieldValue(entry[1],gridDefault->CP) &&
		            entry[2].empty() &&
		            entry[3].empty() &&
		            entry[4].empty() &&
		 fieldValue(entry[5],gridDefault->CD) &&
		 fieldValue(entry[6],gridDefault->PS) &&
		 fieldValue(entry[7],gridDefault->SEID));

#ifdef FFL_DEBUG
  std::cout <<"Default grid-point options: CP = "<< gridDefault->CP
	    << std::endl;
#endif

  STOPP_TIMER("process_GRDSET")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// GRID ///
////////////////////////////////////////////////////////////////////////////////

static int convertDOF (int C)
{
  std::set<int> digits;
  for (; C > 0; C /= 10)
    if (C%10 > 0)
      digits.insert(C%10);

  int status = 0;
  for (int digit : digits)
    switch (digit) {
    case 1: status += 1; break;
    case 2: status += 2; break;
    case 3: status += 4; break;
    case 4: status += 8; break;
    case 5: status += 16; break;
    case 6: status += 32; break;
    }

  return status;
}

bool FFlNastranReader::process_GRID (std::vector<std::string>& entry)
{
  START_TIMER("process_GRID")

  int    ID, CP, CD, PS = 0, SEID = 0;
  double X1 = 0.0, X2 = 0.0, X3 = 0.0;

  if (entry.size() < 8) entry.resize(8,"");

  CONVERT_ENTRY ("GRID",
		 fieldValue(entry[0],ID) &&
		 fieldValue(entry[1],CP) &&
		 fieldValue(entry[2],X1) &&
		 fieldValue(entry[3],X2) &&
		 fieldValue(entry[4],X3) &&
		 fieldValue(entry[5],CD) &&
		 fieldValue(entry[6],PS) &&
		 fieldValue(entry[7],SEID));

#if FFL_DEBUG > 2
  std::cout <<"Grid point, ID = "<< ID;
  if (!entry[1].empty()) std::cout <<", CP = "<< CP;
  std::cout <<", X = "<< X1 <<" "<< X2 <<" "<< X3;
  if (!entry[6].empty()) std::cout <<", PS = "<< PS;
  std::cout << std::endl;
#endif

  if (!entry[1].empty()) nodeCPID[ID] = CP;
  if (!entry[5].empty()) nodeCDID[ID] = CD;

  sizeOK = myLink->addNode(new FFlNode(ID,X1,X2,X3,-convertDOF(PS)));

  STOPP_TIMER("process_GRID")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// INCLUDE
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_INCLUDE (std::vector<std::string>& entry)
{
  std::string& fname = entry.front();
  if (entry.empty()) return false;
  if (fname.empty()) return false;

  const int i = fname.length()-1;
  if (fname[i] == '\'' || fname[i] == '"') fname.erase(i,1);
  if (fname[0] == '\'' || fname[0] == '"') fname.erase(0,1);

  // Save some counter variables on stack before FFlNastranReader::read
  // is invoked recursively on the included file
  std::map<std::string,int> saveIgnored(ignoredBulk);
  std::map<std::string,int> saveSxError(sxErrorBulk);
  int saveCounter = lineCounter;

  ListUI <<"\nReading included file \""<< fname <<"\"\n";
  lineCounter = 0;
  bool ok = read(fname,true);
  ListUI <<"\nDone reading included file \""<< fname <<"\"\n";

  // Restore the counter variables associated with current file from stack
  ignoredBulk = std::move(saveIgnored);
  sxErrorBulk = std::move(saveSxError);
  lineCounter = saveCounter;

  return ok;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// MAT1 ///
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_MAT1 (std::vector<std::string>& entry)
{
  START_TIMER("process_MAT1")

  int    MID = 0;
  double E = 0.0, G = 0.0, NU = 0.0, RHO = 0.0;

  if (entry.size() < 5) entry.resize(5,"");

  CONVERT_ENTRY ("MAT1",
		 fieldValue(entry[0],MID) &&
		 fieldValue(entry[1],E) &&
		 fieldValue(entry[2],G) &&
		 fieldValue(entry[3],NU) &&
		 fieldValue(entry[4],RHO));

  if (E <= 0.0)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Material "<< MID
           <<" has a zero or negative elasticity modulus (E).\n";
    if (G > 0.0 && NU >= 0.0 && NU < 0.5)
    {
      E = G*(2.0+NU+NU);
      ListUI <<"              Resetting to "<< E <<" = 2*G*(1+nu).\n";
    }
    else
      ListUI <<"              This may cause a singular stiffness matrix.\n";
  }
  if (RHO <= 0.0)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Material "<< MID
           <<" has a zero or negative mass density (rho).\n"
           <<"              This may result in a singular mass matrix.\n";
  }

  FFlPMAT* myAtt = CREATE_ATTRIBUTE(FFlPMAT,"PMAT",MID);
  myAtt->youngsModule    = round(E,10);
  myAtt->shearModule     = round(G,10);
  myAtt->poissonsRatio   = round(NU,10);
  myAtt->materialDensity = round(RHO,10);

  if (lastComment.first > 0)
    if (extractNameFromComment(lastComment.second))
      myAtt->setName(lastComment.second);

#ifdef FFL_DEBUG
  myAtt->print();
#endif
  myLink->addAttribute(myAtt);

  STOPP_TIMER("process_MAT1")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// MAT2 ///
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_MAT2 (std::vector<std::string>& entry)
{
  START_TIMER("process_MAT2")

  int    MID = 0;
  double RHO = 0.0;

  std::vector<double> C(6,0.0);

  if (entry.size() < 8) entry.resize(8,"");

  CONVERT_ENTRY ("MAT2",
		 fieldValue(entry[0],MID) &&
		 fieldValue(entry[7],RHO));
  for (size_t i = 0; i < C.size(); i++)
    CONVERT_ENTRY ("MAT2",fieldValue(entry[1+i],C[i]));

  if (C[0] <= 0.0 || C[3] <= 0.0 || C[5] <= 0.0)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Material "<< MID
           <<" has a zero or negative diagonal element Cii \n"
           <<"              This may result in a singular stiffness matrix.\n";
  }
  if (RHO <= 0.0)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Material "<< MID
           <<" has a zero or negative mass density (rho).\n"
           <<"              This may result in a singular mass matrix.\n";
  }

  FFlPMAT2D* myAtt = CREATE_ATTRIBUTE(FFlPMAT2D,"PMAT2D",MID);
  for (size_t i = 0; i < C.size(); i++)
    myAtt->C[i] = round(C[i],10);
  myAtt->materialDensity = round(RHO,10);

  if (lastComment.first > 0)
    if (extractNameFromComment(lastComment.second))
      myAtt->setName(lastComment.second);

#ifdef FFL_DEBUG
  myAtt->print();
#endif
  myLink->addAttribute(myAtt);

  STOPP_TIMER("process_MAT2")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// MAT8 ///
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_MAT8 (std::vector<std::string>& entry)
{
  START_TIMER("process_MAT8")

  int    MID = 0;
  double E1 = 0.0, E2 = 0.0, NU12 = 0.0;
  double G12 = 0.0, G1Z = 0.0, G2Z = 0.0, RHO = 0.0;

  if (entry.size() < 8) entry.resize(8,"");

  CONVERT_ENTRY ("MAT8",
		 fieldValue(entry[0],MID) &&
		 fieldValue(entry[1],E1) &&
		 fieldValue(entry[2],E2) &&
		 fieldValue(entry[3],NU12) &&
		 fieldValue(entry[4],G12) &&
		 fieldValue(entry[5],G1Z) &&
		 fieldValue(entry[6],G2Z) &&
		 fieldValue(entry[7],RHO));

  if (E1 <= 0.0 || E2 <= 0.0)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Material "<< MID
           <<" has a zero or negative elasticity modulus (E).\n"
           <<"              This may result in a singular stiffness matrix.\n";
  }
  if (RHO <= 0.0)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Material "<< MID
           <<" has a zero or negative mass density (rho).\n"
           <<"              This may result in a singular mass matrix.\n";
  }

  FFlPMATSHELL* myAtt = CREATE_ATTRIBUTE(FFlPMATSHELL,"PMATSHELL",MID);
  myAtt->E1   = round(E1,10);
  myAtt->E2   = round(E2,10);
  myAtt->NU12 = round(NU12,10);
  myAtt->G12  = round(G12,10);
  myAtt->G1Z  = round(G1Z,10);
  myAtt->G2Z  = round(G2Z,10);
  myAtt->materialDensity = round(RHO,10);

  if (lastComment.first > 0)
    if (extractNameFromComment(lastComment.second))
      myAtt->setName(lastComment.second);

#ifdef FFL_DEBUG
  myAtt->print();
#endif
  myLink->addAttribute(myAtt);

  STOPP_TIMER("process_MAT8")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// MAT9 ///
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_MAT9 (std::vector<std::string>& entry)
{
  START_TIMER("process_MAT9")

  int    MID = 0;
  double RHO = 0.0;

  std::vector<double> C(21,0.0);

  if (entry.size() < 23) entry.resize(23,"");

  CONVERT_ENTRY ("MAT9",
		 fieldValue(entry[ 0],MID) &&
		 fieldValue(entry[22],RHO));
  for (size_t i = 0; i < C.size(); i++)
    CONVERT_ENTRY ("MAT2",fieldValue(entry[1+i],C[i]));

  if (C[ 0] <= 0.0 || C[ 6] <= 0.0 || C[11] <= 0.0 ||
      C[15] <= 0.0 || C[18] <= 0.0 || C[20] <= 0.0)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Material "<< MID
           <<" has a zero or negative diagonal element Cii \n"
           <<"              This may result in a singular stiffness matrix.\n";
  }
  if (RHO <= 0.0)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Material "<< MID
           <<" has a zero or negative mass density (rho).\n"
           <<"              This may result in a singular mass matrix.\n";
  }

  FFlPMAT3D* myAtt = CREATE_ATTRIBUTE(FFlPMAT3D,"PMAT3D",MID);
  for (size_t i = 0; i < C.size(); i++)
    myAtt->C[i] = round(C[i],10);
  myAtt->materialDensity = round(RHO,10);

  if (lastComment.first > 0)
    if (extractNameFromComment(lastComment.second))
      myAtt->setName(lastComment.second);

#ifdef FFL_DEBUG
  myAtt->print();
#endif
  myLink->addAttribute(myAtt);

  STOPP_TIMER("process_MAT9")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// PCOMP //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_PCOMP (std::vector<std::string>& entry)
{
  START_TIMER("process_PCOMP")

  int    PID = 0;
  double Z0 = 0.0, NSMass = 0.0;

  if (entry.size() < 3) entry.resize(3,"");

  CONVERT_ENTRY ("PCOMP",
		 fieldValue(entry[0],PID) &&
		 fieldValue(entry[1],Z0) &&
		 fieldValue(entry[2],NSMass));

  if (entry.size() < 10)
  {
    ListUI <<"\n *** Error: Composite property PCOMP "<< PID
	   <<" is insufficiently defined.\n";
    STOPP_TIMER("process_PCOMP")
    return false;
  }

  double T_total = 0.0;
  FFlPlyVec plyVec;
  plyVec.reserve((entry.size()-8)/4);
  for (size_t idx = 8; idx+2 < entry.size(); idx += 4)
  {
    FFlPly ply;
    CONVERT_ENTRY ("PCOMP",
		   fieldValue(entry[idx  ],ply.MID) &&
		   fieldValue(entry[idx+1],ply.T) &&
		   fieldValue(entry[idx+2],ply.theta));
    T_total += ply.T;
    plyVec.push_back(ply);
  }

  FFlPCOMP* myAtt = CREATE_ATTRIBUTE(FFlPCOMP,"PCOMP",PID);
  myAtt->plySet.setValue(plyVec);
  myAtt->Z0.setValue(round(Z0 == 0.0 ? -0.5*T_total : Z0, 10));

  myLink->addAttribute(myAtt);
  PCOMPs.insert(PID);

  if (NSMass != 0.0)
  {
    shellPIDnsm.insert(PID);
    myLink->addAttribute(createNSM(PID,NSMass,true));
  }

  STOPP_TIMER("process_PCOMP")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// MOMENT /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_MOMENT (std::vector<std::string>& entry)
{
  START_TIMER("process_MOMENT")

  int    SID = 0, G = 0, CID = 0;
  double M = 0.0;
  FaVec3 N;

  if (entry.size() < 7) entry.resize(7,"");

  CONVERT_ENTRY ("MOMENT",
		 fieldValue(entry[0],SID) &&
		 fieldValue(entry[1],G) &&
		 fieldValue(entry[2],CID) &&
		 fieldValue(entry[3],M) &&
		 fieldValue(entry[4],N[0]) &&
		 fieldValue(entry[5],N[1]) &&
		 fieldValue(entry[6],N[2]));

  N *= M;
#ifdef FFL_DEBUG
  std::cout <<"Concentrated moment, SID = "<< SID <<" --> M = "<< N
	    <<", G = "<< G << std::endl;
#endif

  FFlLoadBase* theLoad = LoadFactory::instance()->create("CMOMENT",SID);
  theLoad->setValue(N);
  theLoad->setTarget(G);
  myLink->addLoad(theLoad);

  if (CID > 0) loadCID[theLoad] = CID;

  STOPP_TIMER("process_MOMENT")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// MPC ////
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_MPC (std::vector<std::string>& entry)
{
  START_TIMER("process_MPC")

  if (entry.size() < 7) entry.resize(7,"");

  int    SID  = 0;
  size_t nMst = 1 + (entry.size()-6)/4;
  std::vector<int>    G(1+nMst,0), C(1+nMst,0);
  std::vector<double> A(1+nMst,0.0);

  CONVERT_ENTRY ("MPC",
		 fieldValue(entry[0],SID) &&
		 fieldValue(entry[1],G[0]) &&
		 fieldValue(entry[2],C[0]) &&
		 fieldValue(entry[3],A[0]) &&
		 fieldValue(entry[4],G[1]) &&
		 fieldValue(entry[5],C[1]) &&
		 fieldValue(entry[6],A[1]));

  size_t i = 8;
  for (size_t j = 2; i+2 < entry.size(); j++)
  {
    CONVERT_ENTRY ("MPC",
		   fieldValue(entry[i+1],G[j]) &&
		   fieldValue(entry[i+2],C[j]) &&
		   (i+3 >= entry.size() || fieldValue(entry[i+3],A[j])));
    i += j%2 ? 5 : 3;
  }

#ifdef FFL_DEBUG
  std::cout <<"Multi-point constraint, SID = "<< SID <<": ";
  for (i = 0; i < G.size(); i++)
  {
    if (i == 0)
      std::cout << A.front();
    else if (A[i] < 0.0)
      std::cout <<" - "<< -A[i];
    else if (A[i] > 0.0)
      std::cout <<" + "<< A[i];
    else
      continue;
    std::cout <<"*("<< G[i] <<","<< C[i] <<")";
  }
  std::cout <<" = 0"<< std::endl;
#endif

  if (fabs(A.front()) < 1.0e-12)
  {
    ListUI <<"\n *** Error: A1 ("<< A.front() <<") must be non-zero"
           <<" for MPC "<< SID <<".\n";
    STOPP_TIMER("process_MPC")
    return false;
  }

  FFl::DepDOFs& masters = myMPCs[G.front()][C.front()];
  masters.reserve(G.size()-1);
  for (i = 1; i < G.size(); i++)
    masters.push_back(FFl::DepDOF(G[i],C[i],-A[i]/A.front()));

  STOPP_TIMER("process_MPC")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// PBAR ///
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_PBAR (std::vector<std::string>& entry)
{
  START_TIMER("process_PBAR")

  int             PID = 0, MID = 0;
  double          C1, C2, D1, D2, E1, E2, F1, F2;
  FFlCrossSection params;

  if (entry.size() < 19) entry.resize(19,"");

  CONVERT_ENTRY ("PBAR",
		 fieldValue(entry[0],PID) &&
		 fieldValue(entry[1],MID) &&
		 fieldValue(entry[2],params.A) &&
		 fieldValue(entry[3],params.Izz) &&
		 fieldValue(entry[4],params.Iyy) &&
		 fieldValue(entry[5],params.J) &&
		 fieldValue(entry[6],params.NSM) &&
		            entry[7].empty() &&
		 fieldValue(entry[8],C1) &&
		 fieldValue(entry[9],C2) &&
		 fieldValue(entry[10],D1) &&
		 fieldValue(entry[11],D2) &&
		 fieldValue(entry[12],E1) &&
		 fieldValue(entry[13],E2) &&
		 fieldValue(entry[14],F1) &&
		 fieldValue(entry[15],F2) &&
		 fieldValue(entry[16],params.K1) &&
		 fieldValue(entry[17],params.K2) &&
		 fieldValue(entry[18],params.Izy));

#ifdef FFL_DEBUG
  std::cout <<"Beam property, ID = "<< PID <<" --> material ID = "<< MID
	    << std::endl;
#endif

  this->insertBeamPropMat("PBAR",PID,MID);
  myLink->addAttribute(createBeamSection(PID,params,lastComment));

  if (params.NSM != 0.0)
  {
    beamPIDnsm.insert(PID);
    myLink->addAttribute(createNSM(PID,params.NSM));
  }

  STOPP_TIMER("process_PBAR")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// PBARL //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_PBARL (std::vector<std::string>& entry)
{
  START_TIMER("process_PBARL")

  int    PID = 0, MID = 0;
  size_t numDim = 0;
  double NSMass = 0.0;

  // The "L" type is not listed in the MSC manual. Should it be?
  static std::vector<std::string> Types4({
    "HAT", "CHAN", "CHAN1", "BOX", "CHAN2", "CROSS",
    "T", "T1", "T2", "Z", "H", "I1"
  });

  // ======================= Read & check =======================

  if (entry.size() < 4) entry.resize(4,"");

  // Read IDs and check cross-section section group/type
  CONVERT_ENTRY ("PBARL",
		 fieldValue(entry[0],PID) &&
		 fieldValue(entry[1],MID));

  if (!entry[2].empty() && entry[2] != "MSCBML0")
  {
    ListUI <<"\n *** Error: Beam property "<< PID <<" cross-section group "
           << entry[2] <<" is not supported. Use standard MSCBML0.\n";
    STOPP_TIMER("process_PBARL")
    return false;
  }

  std::string Type = entry[3];
  if (Type.empty())
  {
    ListUI <<"\n *** Error: Cross section type not specified"
           <<" for beam property "<< PID <<".\n";
    STOPP_TIMER("process_PBARL")
    return false;
  }
  else if (Type == "ROD")
    numDim = 1;
  else if (Type == "TUBE" || Type == "BAR")
    numDim = 2;
  else if (Type == "HEXA")
    numDim = 3;
  else if (Type == "BOX1" || Type == "I")
    numDim = 6;
  else if (std::find(Types4.begin(),Types4.end(),Type) != Types4.end())
    numDim = 4;

  if (entry.size() < 9+numDim) entry.resize(9+numDim,"");

  std::vector<double> Dim(numDim,0.0);
  for (size_t i = 8; i < 8+numDim; i++)
    CONVERT_ENTRY ("PBARL", fieldValue(entry[i],Dim[i-8]) && !entry[i].empty());
  CONVERT_ENTRY ("PBARL", fieldValue(entry[8+numDim],NSMass));

  // ======================= Process data =======================

#ifdef FFL_DEBUG
  std::cout <<"Beam property, ID = "<< PID <<" --> material ID = "<< MID
	    << std::endl;
#endif

  FFlCrossSection params(Type,Dim);
  if (params.A > 0.0)
  {
    this->insertBeamPropMat("PBARL",PID,MID);
    myLink->addAttribute(createBeamSection(PID,params,lastComment));
  }
  else
    ListUI <<"            Error occurred when processing PBARL "<< PID <<".\n";

  if (NSMass != 0.0)
  {
    beamPIDnsm.insert(PID);
    myLink->addAttribute(createNSM(PID,NSMass));
  }

  STOPP_TIMER("process_PBARL")
  return params.A > 0.0;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// PBEAM //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_PBEAM (std::vector<std::string>& entry)
{
  START_TIMER("process_PBEAM")

  // Lambda function checking for presence of the Stress Output option field
  auto&& isSOFIELD = [](const std::string& field) -> bool
  {
    return field == "YES" || field == "YESA" || field == "NO";
  };

  // Lambda function for averaging a cross section property
  auto&& averageBS = [](double& A, double B, const char* s) -> short int
  {
    if (fabs(A-B) <= 1.0e-9*(fabs(A)+fabs(B)))
      return 1;

    ListUI <<"           "<< s <<"(A) = "<< A <<" "<< s <<"(B) = "<< B;
    A = 0.5*(A+B);
    ListUI <<" --> "<< s <<" = "<< A <<"\n";
    return 0;
  };

  int             PID = 0, MID = 0;
  FFlCrossSection params(1.0);

  if (entry.size() < 17) entry.resize(17,"");

  CONVERT_ENTRY ("PBEAM",
		 fieldValue(entry[0],PID) &&
		 fieldValue(entry[1],MID) &&
		 fieldValue(entry[2],params.A) &&
		 fieldValue(entry[3],params.Izz) &&
		 fieldValue(entry[4],params.Iyy) &&
		 fieldValue(entry[5],params.Izy) &&
		 fieldValue(entry[6],params.J) &&
		 fieldValue(entry[7],params.NSM));

  // Check if this is a tapered beam, and
  // find the index of the first SO-field (stress output option), if any
  size_t iLast = 16;
  char hasTapering = 0;
  if (isSOFIELD(entry[16]))
    hasTapering = 'y';
  else if (isSOFIELD(entry[8]))
  {
    hasTapering = 'y';
    iLast = 8; // First continuation is omitted
    nNotes++;
    ListUI <<"\n   * Note: Stress output option \""<< entry[8]
           <<"\" was found in field 9 of a PBEAM bulk entry.\n           "
           <<"This is non-standard. Assuming fields 9-18 are zero.\n           "
           <<"Please verify that corresponding PBEAMSECTION entry "<< PID
           <<"\n           in the corresponding ftl-file is correct.\n";
  }
  else if (entry.size() <= 17)
    iLast = 0; // No tapering

  if (hasTapering) // Find the properties at end B
    while (iLast < entry.size() && isSOFIELD(entry[iLast]))
    {
      double endB[11]; endB[0] = 0.0;
      memcpy(endB+1,&params.A,10*sizeof(double));
      for (size_t i = 1; i < 8 && iLast+i < entry.size(); i++)
        if (!entry[iLast+i].empty())
          CONVERT_ENTRY ("PBEAM",fieldValue(entry[iLast+i],endB[i-1]));

      iLast += entry[iLast] == "YES" ? 16 : 8;
      if (endB[0] > 0.999 && hasTapering == 'y')
        hasTapering = (averageBS(params.A  ,endB[1],"A") &
                       averageBS(params.Izz,endB[2],"I1") &
                       averageBS(params.Iyy,endB[3],"I2") &
                       averageBS(params.Izy,endB[4],"I12") &
                       averageBS(params.J  ,endB[5],"J") &
                       averageBS(params.NSM,endB[6],"NSM")) ? 'N' : 'Y';
    }

  if (hasTapering == 'Y')
  {
    nNotes++;
    ListUI <<"   * Note: Beam property "<< PID <<" has tapering.\n           "
           <<"The properties specified at the two end points are averaged "
           <<"(see above).\n";
  }
  else if (hasTapering == 'y')
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Beam property "<< PID <<" has tapering,\n"
           <<"                but properties at end B were not found.\n";
  }
#ifdef FFL_DEBUG
  else if (hasTapering == 'N')
    std::cout <<"Beam property with ID = "<< PID <<" has tapering,\n"
              <<"but all properties at both ends are equal."<< std::endl;
#endif

  if (iLast > 0 && iLast < entry.size())
  {
    double S1, S2, NSIA, NSIB, CWIA, CWIB, M1A, M2A, M1B, M2B;
    double N1A = 0.0, N2A = 0.0, N1B = 0.0, N2B = 0.0;

    // The last two continuations
    if (entry.size() < iLast+16) entry.resize(iLast+16,"");
    CONVERT_ENTRY ("PBEAM",
		   fieldValue(entry[iLast+0],params.K1) &&
		   fieldValue(entry[iLast+1],params.K2) &&
		   fieldValue(entry[iLast+2],S1) &&
		   fieldValue(entry[iLast+3],S2) &&
		   fieldValue(entry[iLast+4],NSIA) &&
		   fieldValue(entry[iLast+5],NSIB) &&
		   fieldValue(entry[iLast+6],CWIA) &&
		   fieldValue(entry[iLast+7],CWIB) &&
		   fieldValue(entry[iLast+8],M1A) &&
		   fieldValue(entry[iLast+9],M2A) &&
		   fieldValue(entry[iLast+10],M1B) &&
		   fieldValue(entry[iLast+11],M2B) &&
		   fieldValue(entry[iLast+12],N1A) &&
		   fieldValue(entry[iLast+13],N2A) &&
		   fieldValue(entry[iLast+14],N1B) &&
		   fieldValue(entry[iLast+15],N2B));

    // Shear center offset with respect to the neutral axis (averaged)
    params.S1 = -0.5*(N1A+N1B);
    params.S2 = -0.5*(N2A+N2B);
  }

#ifdef FFL_DEBUG
  std::cout <<"Beam property, ID = "<< PID <<" --> material ID = "<< MID
	    << std::endl;
#endif

  this->insertBeamPropMat("PBEAM",PID,MID);
  myLink->addAttribute(createBeamSection(PID,params,lastComment));

  if (params.NSM != 0.0)
  {
    beamPIDnsm.insert(PID);
    myLink->addAttribute(createNSM(PID,params.NSM));
  }

  STOPP_TIMER("process_PBEAM")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// PBEAML /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_PBEAML (std::vector<std::string>& entry)
{
  START_TIMER("process_PBEAML")

  int    PID = 0, MID = 0;
  size_t i, numDim = 0;
  double NSMass = 0.0;

  static std::vector<std::string> Types4({
    "HAT", "CHAN", "CHAN1", "BOX", "CHAN2", "CROSS",
    "T", "T1", "T2", "L", "Z", "H", "I1"
  });

  // ======================= Read & check =======================

  if (entry.size() < 4) entry.resize(4,"");

  // Read IDs and check cross-section section group/type
  CONVERT_ENTRY ("PBEAML",
		 fieldValue(entry[0],PID) &&
		 fieldValue(entry[1],MID));

  if (!entry[2].empty() && entry[2] != "MSCBML0")
  {
    ListUI <<"\n *** Error: Beam property "<< PID <<" cross-section group "
           << entry[2] <<" is not supported. Use standard MSCBML0.\n";
    STOPP_TIMER("process_PBEAML")
    return false;
  }

  std::string Type = entry[3];
  if (Type.empty())
  {
    ListUI <<"\n *** Error: Cross section type not specified"
           <<" for beam property "<< PID <<".\n";
    STOPP_TIMER("process_PBEAML")
    return false;
  }
  else if (Type == "ROD")
    numDim = 1;
  else if (Type == "TUBE" || Type == "BAR")
    numDim = 2;
  else if (Type == "HEXA")
    numDim = 3;
  else if (Type == "BOX1" || Type == "I")
    numDim = 6;
  else if (std::find(Types4.begin(),Types4.end(),Type) != Types4.end())
    numDim = 4;

  // Read end A data
  if (entry.size() < 9+numDim) entry.resize(9+numDim,"");

  std::vector<double> Dim(numDim,0.0);
  for (i = 8; i < 8+numDim; i++)
    CONVERT_ENTRY ("PBEAML",fieldValue(entry[i],Dim[i-8]) && !entry[i].empty());
  CONVERT_ENTRY ("PBEAML",fieldValue(entry[8+numDim],NSMass));

  // Find starting point of end B data
  size_t endBStart = 9+numDim;
  size_t inc       = 3+numDim;
  while (endBStart+inc < entry.size()) endBStart += inc;

  // Read end B data - substitute by data of end A if field is empty
  if (entry.size() < endBStart+inc) entry.resize(endBStart+inc,"");

  std::vector<double> DimB(Dim);
  double NsmB = NSMass;
  for (i = 0; i < numDim; i++)
    CONVERT_ENTRY ("PBEAML", fieldValue(entry[endBStart+2+i],DimB[i]));
  CONVERT_ENTRY ("PBEAML",fieldValue(entry[endBStart+2+numDim],NsmB));

  // ======================= Process data =======================

  if (Dim != DimB || NSMass != NsmB || endBStart > numDim+9)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Beam property "<< PID <<" has tapering.\n"
           <<"              In the current version, the properties specified at"
           <<" end A\n              are used throughout each beam element.\n";
  }

#ifdef FFL_DEBUG
  std::cout <<"Beam property, ID = "<< PID <<" --> material ID = "<< MID
	    << std::endl;
#endif

  FFlCrossSection params(Type,Dim);
  if (params.A > 0.0)
  {
    this->insertBeamPropMat("PBEAML",PID,MID);
    myLink->addAttribute(createBeamSection(PID,params,lastComment));
  }
  else
    ListUI <<"            Error occurred when processing PBEAML "<< PID <<".\n";

  if (NSMass != 0.0)
  {
    beamPIDnsm.insert(PID);
    myLink->addAttribute(createNSM(PID,NSMass));
  }

  STOPP_TIMER("process_PBEAML")
  return params.A > 0.0;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// PBUSH //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_PBUSH (std::vector<std::string>& entry)
{
  START_TIMER("process_PBUSH")

  size_t i;
  int    PID = 0;

  std::vector<double> K(6,0.0);

  if (entry.size() < 4)
    entry.resize(4,"");
  else if (entry.size() < 8)
    entry.resize(8,"");

  for (i = 0; i < entry.size(); i += 8)
    if (entry[i+1] == "K")
    {
      if (entry.size() < i+8) entry.resize(i+8,"");
      CONVERT_ENTRY ("PBUSH",
		     fieldValue(entry[0],PID) &&
		     fieldValue(entry[i+2],K[0]) &&
		     fieldValue(entry[i+3],K[1]) &&
		     fieldValue(entry[i+4],K[2]) &&
		     fieldValue(entry[i+5],K[3]) &&
		     fieldValue(entry[i+6],K[4]) &&
		     fieldValue(entry[i+7],K[5]));
    }
    else if (entry[i+1] == "B")
    {
      nWarnings++;
      ListUI <<"\n  ** Warning: Bushing property "<< PID
             <<" has force-per-velocity damping (ignored).\n";
    }
    else if (entry[i+1] == "GE")
    {
      nWarnings++;
      ListUI <<"\n  ** Warning: Bushing property "<< PID
             <<" has structural damping (ignored).\n";
    }

#ifdef FFL_DEBUG
  std::cout <<"Bushing property, ID = "<< PID <<" --> K =";
  for (double k : K) std::cout <<" "<< k;
  std::cout << std::endl;
#endif

  FFlPBUSHCOEFF* myAtt = CREATE_ATTRIBUTE(FFlPBUSHCOEFF,"PBUSHCOEFF",PID);
  for (i = 0; i < 6; i++) myAtt->K[i] = round(K[i],10);

  if (lastComment.first > 0)
    if (extractNameFromComment(lastComment.second))
      myAtt->setName(lastComment.second);

  myLink->addAttribute(myAtt);

  STOPP_TIMER("process_PBUSH")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// PELAS //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_PELAS (std::vector<std::string>& entry)
{
  START_TIMER("process_PELAS")

  int    PID = 0;
  double K = 0.0, GE = 0.0;

  if (entry.size() < 4)
    entry.resize(4,"");
  else if (entry.size() < 8)
    entry.resize(8,"");

  for (size_t i = 0; i+3 < entry.size(); i += 4)
  {
    CONVERT_ENTRY ("PELAS",
		   fieldValue(entry[i+0],PID) &&
		   fieldValue(entry[i+1],K) &&
		   fieldValue(entry[i+2],GE));

#ifdef FFL_DEBUG
    std::cout <<"Spring property, ID = "<< PID <<" --> K = "<< K << std::endl;
#endif
    if (GE != 0.0)
    {
      nWarnings++;
      ListUI <<"\n  ** Warning: Spring property "<< PID
             <<" has structural damping (ignored).\n";
    }

    propK[PID] = K;
  }

  STOPP_TIMER("process_PELAS")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// PLOAD2 /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_PLOAD2 (std::vector<std::string>& entry)
{
  START_TIMER("process_PLOAD2")

  int    SID = 0;
  double P = 0.0;

  std::vector<int> EID;

  if (entry.size() < 3) entry.resize(3,"");

  CONVERT_ENTRY ("PLOAD2",
		 fieldValue(entry[0],SID) &&
		 fieldValue(entry[1],P));

  if (entry.size() > 4 && entry[3] == "THRU")
  {
    int EID1 = 0, EID2 = 0;
    CONVERT_ENTRY ("PLOAD2",
		   fieldValue(entry[2],EID1) &&
		   fieldValue(entry[4],EID2) &&
		   EID2 > EID1);

    EID.reserve(EID2-EID1+1);
    for (int e = EID1; e <= EID2; e++) EID.push_back(e);
  }
  else
  {
    EID.resize(entry.size()-2,0);
    for (size_t i = 2; i < entry.size(); i++)
      CONVERT_ENTRY ("PLOAD2",fieldValue(entry[i],EID[i-2]));
  }

#ifdef FFL_DEBUG
  std::cout <<"Pressure load, SID = "<< SID <<" --> P = "<< P <<", EID";
  for (int e : EID) std::cout <<" "<< e;
  std::cout << std::endl;
#endif

  FFlLoadBase* theLoad = LoadFactory::instance()->create("SURFLOAD",SID);
  theLoad->setValue({round(P,10)});
  theLoad->setTarget(EID);
  myLink->addLoad(theLoad);

  STOPP_TIMER("process_PLOAD2")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// PLOAD4 /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_PLOAD4 (std::vector<std::string>& entry)
{
  START_TIMER("process_PLOAD4")

  int    SID = 0, CID = 0, G1 = 0, G34 = 0;
  FaVec3 N;

  std::vector<int>    EID(1,0);
  std::vector<double> P(1,0.0);

  if (entry.size() < 3) entry.resize(3,"");

  CONVERT_ENTRY ("PLOAD4",
		 fieldValue(entry[0],SID) &&
		 fieldValue(entry[1],EID[0]) &&
		 fieldValue(entry[2],P[0]));

  if (entry.size() > 3)
  {
    P.resize(entry.size() > 5 ? 4 : entry.size()-2,P[0]);
    for (size_t i = 1; i < P.size(); i++)
      CONVERT_ENTRY ("PLOAD4",fieldValue(entry[2+i],P[i]));
    while (P.size() > 1 && P.back() == P.front()) P.pop_back();
  }

  if (entry.size() > 7)
  {
    if (entry[6] == "THRU")
    {
      int EID2 = 0;
      CONVERT_ENTRY ("PLOAD4",
		     fieldValue(entry[7],EID2) &&
		     EID2 > EID.front());
      EID.reserve(EID2-EID.front()+1);
      for (int e = EID.front()+1; e <= EID2; e++) EID.push_back(e);
    }
    else
      CONVERT_ENTRY ("PLOAD4",
		     fieldValue(entry[6],G1) &&
		     fieldValue(entry[7],G34));
  }

  if (entry.size() > 9)
  {
    if (entry.size() < 12) entry.resize(12,"");

    CONVERT_ENTRY ("PLOAD4",
		   fieldValue(entry[8],CID) &&
		   fieldValue(entry[9],N[0]) &&
		   fieldValue(entry[10],N[1]) &&
		   fieldValue(entry[11],N[2]));
  }

#ifdef FFL_DEBUG
  std::cout <<"Pressure load, SID = "<< SID <<" --> P = "<< P[0] <<", EID";
  for (int e : EID) std::cout <<" "<< e;
  std::cout << std::endl;
#endif

  FFlLoadBase* theLoad;
  if (G1 > 0)
  {
    // Solid face load
    theLoad = LoadFactory::instance()->create("FACELOAD",SID);
    loadFace[theLoad] = std::make_pair(G1,G34);
  }
  else
    // Shell surface load
    theLoad = LoadFactory::instance()->create("SURFLOAD",SID);

  for (double& p : P) round(p,10);
  theLoad->setValue(P);
  theLoad->setTarget(EID);
  myLink->addLoad(theLoad);

  if (!N.isZero())
  {
    // This load has an orientation vector
    FFlPORIENT* myOr = CREATE_ATTRIBUTE(FFlPORIENT,"PORIENT",EID.front());
    myOr->directionVector = N.normalize().round(10);

    if (CID > 0)
    {
      nWarnings++;
      ListUI <<"\n  ** Warning: Pressure load "<< SID
             <<" on element "<< EID.front()
             <<"\n              has a direction vector specified in local"
	     <<" coordinate system "<< CID
	     <<"\n              This is not implemented yet"
	     <<" (assuming global system instead).\n";
    }

    theLoad->setAttribute("PORIENT",myLink->addUniqueAttribute(myOr));
  }

  STOPP_TIMER("process_PLOAD4")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// PROD ///
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_PROD (std::vector<std::string>& entry)
{
  START_TIMER("process_PROD")

  int             PID = 0, MID = 0;
  FFlCrossSection params;

  if (entry.size() < 6) entry.resize(6,"");

  CONVERT_ENTRY ("PROD",
		 fieldValue(entry[0],PID) &&
		 fieldValue(entry[1],MID) &&
		 fieldValue(entry[2],params.A) &&
		 fieldValue(entry[3],params.J) &&
		 fieldValue(entry[5],params.NSM));

#ifdef FFL_DEBUG
  std::cout <<"Rod property, ID = "<< PID <<" --> material ID = "<< MID
	    << std::endl;
#endif

  this->insertBeamPropMat("PROD",PID,MID);
  myLink->addAttribute(createBeamSection(PID,params,lastComment));

  if (params.NSM != 0.0)
  {
    beamPIDnsm.insert(PID);
    myLink->addAttribute(createNSM(PID,params.NSM));
  }

  STOPP_TIMER("process_PROD")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// PSHELL /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_PSHELL (std::vector<std::string>& entry)
{
  START_TIMER("process_PSHELL")

  int    PID = 0, MID1 = 0, MID2 = 0, MID3, MID4;
  double T, Iratio, TSratio, Z1, Z2, NSMass = 0.0;

  if (entry.size() < 11) entry.resize(11,"");

  CONVERT_ENTRY ("PSHELL",
		 fieldValue(entry[0],PID) &&
		 fieldValue(entry[1],MID1) &&
		 fieldValue(entry[2],T) &&
		 fieldValue(entry[3],MID2) &&
		 fieldValue(entry[4],Iratio) &&
		 fieldValue(entry[5],MID3) &&
		 fieldValue(entry[6],TSratio) &&
		 fieldValue(entry[7],NSMass) &&
		 fieldValue(entry[8],Z1) &&
		 fieldValue(entry[9],Z2) &&
		 fieldValue(entry[10],MID4));

#ifdef FFL_DEBUG
  std::cout <<"Shell property, ID = "<< PID <<" --> material ID = "<< MID1
	    << std::endl;
#endif

  bool ok = true;
  if (MID1 > 0)
    propMID[FFlTypeInfoSpec::SHELL_ELM][PID] = MID1;
  else if (MID2 > 0)
    propMID[FFlTypeInfoSpec::SHELL_ELM][PID] = MID2;
  else
  {
    ok = false;
    ListUI <<"\n *** Error: PSHELL "<< PID
	   <<" lacks reference to a material property (MID1 and MID2).\n";
  }

  if (MID1*MID2 > 0 && MID2 != MID1)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: PSHELL "<< PID
	   <<" refer to different material properties\n"
	   <<"              for membrane (MID1=" << MID1
	   <<") and bending (MID2="<< MID2 <<")\n              "
	   <<"This is not supported, MID1 will be used also for bending.\n";
  }

  if (!entry[2].empty())
  {
    FFlPTHICK* myAtt = CREATE_ATTRIBUTE(FFlPTHICK,"PTHICK",PID);
    myAtt->thickness = round(T,10);
    if (myAtt->thickness.getValue() <= 0.0)
    {
      nWarnings++;
      ListUI <<"\n  ** Warning: Invalid thickness value, "
             << myAtt->thickness.getValue() <<", for PTHICK "<< PID <<".\n  "
             <<"            This may result in a singular stiffness matrix.\n";
    }

    if (lastComment.first > 0)
      if (extractNameFromComment(lastComment.second))
	myAtt->setName(lastComment.second);

#ifdef FFL_DEBUG
    myAtt->print();
#endif
    myLink->addAttribute(myAtt);
    PTHICKs.insert(PID);
  }

  if (NSMass != 0.0)
  {
    shellPIDnsm.insert(PID);
    myLink->addAttribute(createNSM(PID,NSMass,true));
  }

  STOPP_TIMER("process_PSHELL")
  return ok;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// PSOLID /
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_PSOLID (std::vector<std::string>& entry)
{
  START_TIMER("process_PSOLID")

  int PID = 0, MID = 0;

  if (entry.size() < 2) entry.resize(2,"");

  CONVERT_ENTRY ("PSOLID",
		 fieldValue(entry[0],PID) &&
		 fieldValue(entry[1],MID));

#ifdef FFL_DEBUG
  std::cout <<"Solid property, ID = "<< PID <<" --> material ID = "<< MID
	    << std::endl;
#endif

  propMID[FFlTypeInfoSpec::SOLID_ELM][PID] = MID;

  STOPP_TIMER("process_PSOLID")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// PWELD //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_PWELD (std::vector<std::string>& entry)
{
  START_TIMER("process_PWELD")

  int PID = 0, MID = 0;
  double D = 0.0;

  if (entry.size() < 8) entry.resize(8,"");

  CONVERT_ENTRY ("PWELD",
		 fieldValue(entry[0],PID) &&
		 fieldValue(entry[1],MID) &&
		 fieldValue(entry[2],D));

  if (entry[5] == "OFF")
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: MSET = \"OFF\" is not implemented,"
           <<" for CWELD property "<< PID
           <<".\n              MSET = \"ON\" is assumed instead.\n";
  }
  if (entry[7] == "SPOT")
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Connection type \"SPOT\" is not implemented,"
           <<" for CWELD property "<< PID <<" (ignored).\n";
  }

#ifdef FFL_DEBUG
  std::cout <<"Weld property, ID = "<< PID <<" --> material ID = "<< MID
	    << std::endl;
#endif

  this->insertBeamPropMat("PWELD",PID,MID);

  FFlCrossSection params("ROD",{0.5*D});
  myLink->addAttribute(createBeamSection(PID,params,lastComment));

  STOPP_TIMER("process_PWELD")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// QSET1 //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_QSET1 (std::vector<std::string>& entry)
{
  START_TIMER("process_QSET1")

  int node1 = 0, node2 = 0, dofs = 0;

  if (!entry.empty())
    CONVERT_ENTRY ("QSET1",fieldValue(entry.front(),dofs));

  for (size_t i = 1; i < entry.size(); i++)
    if (entry[i] == "THRU")
      node1 = node2;
    else
    {
      node2 = 0;
      CONVERT_ENTRY ("QSET1",fieldValue(entry[i],node2));
      if (node1 > 0)
      {
        myLink->addComponentModes(node2-node1+1);
        node1 = 0;
      }
      else if (node2 > 0)
        myLink->addComponentModes(1);
    }

  STOPP_TIMER("process_QSET1")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// RBAR ///
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_RBAR (std::vector<std::string>& entry)
{
  START_TIMER("process_RBAR")

  bool replaceByRGD = false;
  int  EID, PID, CNA = 0, CNB = 0, CMA = 0, CMB = 0;

  std::vector<int> G(2,0);

  if (entry.size() < 7) entry.resize(7,"");

  CONVERT_ENTRY ("RBAR",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],G[0]) &&
		 fieldValue(entry[2],G[1]) &&
		 fieldValue(entry[3],CNA) &&
		 fieldValue(entry[4],CNB) &&
		 fieldValue(entry[5],CMA) &&
		 fieldValue(entry[6],CMB));

  if (entry.front().empty())
    EID = myLink->getNewElmID();

  // Check if the RBAR can be represented by a two-noded RBE2 instead
  if (sortDOFs(CNA) == 123456 && CNB == 0 && CMA == 0)
  {
    replaceByRGD = true; // Node B is the dependent node
    CMB = sortDOFs(CMB);
  }
  else if (sortDOFs(CNB) == 123456 && CNA == 0 && CMB == 0)
  {
    replaceByRGD = true; // Node A is the dependent node (swap nodes)
    PID = G[0]; G[0] = G[1]; G[1] = PID;
    CMB = sortDOFs(CMA);
  }

  if (G[0] == G[1])
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Ignoring invalid RBAR element "<< EID
           <<", both ends are connected to the same node "<< G[0] <<".\n";
  }
  else if (replaceByRGD) // Create an RGD element for this RBAR
  {
    if (CMB <= 0 || CMB >= 123456)
      PID = 0; // All DOFs in the dependent node are coupled
    else
    {
      FFlPRGD* myAtt = CREATE_ATTRIBUTE(FFlPRGD,"PRGD",EID);
      myAtt->dependentDofs = CMB;
      PID = myLink->addUniqueAttribute(myAtt);
    }
    sizeOK = myLink->addElement(createElement("RGD",EID,G,PID));
  }
  else // This RBAR element has its dependent DOFs distributed on both nodes
  {
    FFlPRBAR* myAtt = CREATE_ATTRIBUTE(FFlPRBAR,"PRBAR",EID);
    myAtt->CNA = sortDOFs(CNA);
    myAtt->CNB = sortDOFs(CNB);
    myAtt->CMA = sortDOFs(CMA);
    myAtt->CMB = sortDOFs(CMB);
    PID = myLink->addUniqueAttribute(myAtt);
    sizeOK = myLink->addElement(createElement("RBAR",EID,G,PID));
  }

  STOPP_TIMER("process_RBAR")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// RBE2 ///
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_RBE2 (std::vector<std::string>& entry)
{
  START_TIMER("process_RBE2")

  int EID, PID, GM, CM = 0;

  if (entry.size() < 3)
    entry.resize(3,"");
  else if (entry.size() > 3 && entry.back().find('.') != std::string::npos)
    entry.pop_back(); // Ignore thermal expansion coefficient in the last field

  std::vector<int> G(1,0);
  G.reserve(entry.size()-2);

  CONVERT_ENTRY ("RBE2",
		 fieldValue(entry[0],EID) &&
		 fieldValue(entry[1],G[0]) &&
		 fieldValue(entry[2],CM));

  if (entry.front().empty())
    EID = myLink->getNewElmID();

  for (size_t i = 3; i < entry.size(); i++)
    if (!entry[i].empty())
    {
      CONVERT_ENTRY ("RBE2",fieldValue(entry[i],GM));
      if (GM == G[0])
      {
        nWarnings++;
        ListUI <<"\n  ** Warning: Ignoring node "<< GM
               <<" as dependent node for RBE2 element "<< EID
               <<" since it is also specified as the independent node.\n";
      }
      else
        G.push_back(GM);
    }

  if (G.size() < 2)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: One-noded RBE2 element "<< EID <<" (ignored).\n";
    STOPP_TIMER("process_RBE2")
    return true;
  }

  if (sortDOFs(CM) == 123456)
    PID = 0; // All DOFs are coupled
  else
  {
    FFlPRGD* myAtt = CREATE_ATTRIBUTE(FFlPRGD,"PRGD",EID);
    myAtt->dependentDofs = sortDOFs(CM);
    PID = myLink->addUniqueAttribute(myAtt);
  }
  sizeOK = myLink->addElement(createElement("RGD",EID,G,PID));

  STOPP_TIMER("process_RBE2")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// RBE3 ///
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_RBE3 (std::vector<std::string>& entry)
{
  // Lambda function checking if given target integer contains given digit
  auto&& isDigitIn = [](int target, const int digit)
  {
    while (target > 0)
      if (target%10 == digit)
        return true;
      else
        target /= 10;
    return false;
  };

  // Lambda function checking if given target integer contains a set of digits.
  auto&& isSubset = [&isDigitIn](int target, int value)
  {
    while (value > 0)
      if (!isDigitIn(target,value%10))
        return false;
      else
        value /= 10;
    return true;
  };

  // Lambda function removong a set of digits from given target integer.
  auto&& setMinusIfSubset = [&isSubset,&isDigitIn](int& target, int value)
  {
    if (!isSubset(target,value))
      return false;

    int digit, result = 0, expon = 1;
    while (target > 0)
    {
      digit = target%10;
      if (!isDigitIn(value,digit))
        result += digit*expon;
      expon *= 10;
      target /= 10;
    }

    target = result;
    return true;
  };

  struct NodeGroup
  {
    double WT = 0.0;
    int    C  = 0;
    std::vector<int> G;
  };

  START_TIMER("process_RBE3")

  int n, EID, REFC = 0;

  std::vector<int>       C, G(1,0);
  std::vector<NodeGroup> mGroup;
  NodeGroup              aGroup;

  if (entry.size() < 4) entry.resize(4,"");

  // Parse element ID and the reference grid point and components
  CONVERT_ENTRY ("RBE3",
		 fieldValue(entry[0],EID) &&
		            entry[1].empty() &&
		 fieldValue(entry[2],G[0]) &&
		 fieldValue(entry[3],REFC));

  if (entry.front().empty())
    EID = myLink->getNewElmID();

  if (REFC <= 0)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: RBE3 element "<< EID <<" has a zero REFC field.\n"
           <<"              This element will not define any constraints and is"
           <<" therefore ignored.\n";
    STOPP_TIMER("process_RBE3")
    return true;
  }

  // Parse the node group definition and store temporarily in aGroup
  for (size_t i = 4; i < entry.size(); i++)
    if (entry[i] == "UM")
    {
      nWarnings++;
      ListUI <<"\n  ** Warning: RBE3 element "<< EID <<" has the UM field.\n"
             <<"              This is currently not supported (ignored).\n";
      break;
    }
    else if (FFlFieldBase::parseNumericField(n,entry[i],true)) // is field int?
      aGroup.G.push_back(n); // A grid point number was read
    else if (i+1 < entry.size())
    {
      // The next field was not an integer, it must be a weight factor then
      if (aGroup.G.size() > 0)
      {
        // Store the previously read node group in the array first
        if (aGroup.WT != 0.0 && aGroup.C > 0)
          mGroup.push_back(aGroup);
        else
        {
          nWarnings++;
          ListUI <<"\n  ** Warning: RBE3 element "<< EID
                 <<" has a zero weight factor field\n             "
                 <<" and/or the associated component number field is zero.\n  "
                 <<"            This weight factor set is therefore ignored.\n";
        }
        aGroup.G.clear();
      }
      // Then parse weighting factor and associated component number
      CONVERT_ENTRY ("RBE3",
		     fieldValue(entry[i++],aGroup.WT) &&
		     fieldValue(entry[i],aGroup.C));
    }

  // Store the last read node group in the array
  if (aGroup.G.size() > 0)
  {
    if (aGroup.WT != 0.0 && aGroup.C > 0)
      mGroup.push_back(aGroup);
    else
    {
      nWarnings++;
      ListUI <<"\n  ** Warning: RBE3 element "<< EID
             <<" has a zero weight factor field\n             "
             <<" and/or the associated component number field is zero.\n"
             <<"              This weight factor set is therefore ignored.\n";
    }
  }

  // Find the element nodes by inserting unique numbers from the node groups
  // and the number of different component numbers for the weigthed grid points
  for (const NodeGroup& group : mGroup)
  {
    for (int node : group.G)
      if (std::find(G.begin(),G.end(),node) == G.end())
        G.push_back(node);

    if (std::find(C.begin(),C.end(),group.C) == C.end())
    {
      for (int& comp : C)
        if (setMinusIfSubset(comp,group.C)) break;
      C.push_back(group.C);
    }
  }

  size_t nCol = G.size()-1;
  size_t nRow = C.size();
  if (nCol*nRow < 1)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Empty RBE3 element "<< EID <<" (ignored).\n";
    STOPP_TIMER("process_RBE3")
    return true;
  }

  // Now construct the weigthing matrix from the temporary node group data
  std::vector<double> tmpW(nRow*nCol,0.0);
  for (const NodeGroup& group : mGroup)
    for (size_t i = 0; i < nRow; i++)
      if (isSubset(group.C,C[i]))
        for (size_t j = 0; j < nCol; j++)
          for (int node : group.G)
            if (node == G[1+j])
            {
              tmpW[i*nCol+j] = round(group.WT,10);
              break;
            }

  FFlPWAVGM* newAtt = CREATE_ATTRIBUTE(FFlPWAVGM,"PWAVGM",EID);
  newAtt->refC = sortDOFs(REFC);
  newAtt->weightMatrix.data().swap(tmpW);

  // Compute the component indices
  for (size_t i = 0; i < nRow; i++)
    for (size_t j = 0; j < 6; j++)
      if (isDigitIn(C[i],j+1))
        newAtt->indC[j] = i*nCol+1;

  FFlAttributeBase* myAtt = newAtt;
  int PID = myLink->addUniqueAttributeCS(myAtt);
  sizeOK = myLink->addElement(createElement("WAVGM",EID,G,PID));

  STOPP_TIMER("process_RBE3")
  return sizeOK;
}


////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// SET1 //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_SET1 (std::vector<std::string>& entry)
{
  START_TIMER("process_SET1")

  int ID1 = 0, ID2 = 0, SID = 0;

  if (!entry.empty())
    CONVERT_ENTRY ("SET1",fieldValue(entry.front(),SID));

  std::vector<int> IDs;
  IDs.reserve(entry.size()-1);
  for (size_t i = 1; i < entry.size(); i++)
    if (entry[i] == "THRU")
      ID1 = ID2;
    else
    {
      ID2 = 0;
      CONVERT_ENTRY ("SET1",fieldValue(entry[i],ID2));
      if (ID1 > 0)
      {
        for (int i = ID1+1; i <= ID2; i++)
          IDs.push_back(i);
        ID1 = 0;
      }
      else if (ID2 > 0)
        IDs.push_back(ID2);
    }

  int oldNotes = nNotes;
  if (SID > 0 && !IDs.empty())
  {
    lastGroup = new FFlGroup(SID,"Nastran SET");
    for (int eID : IDs)
      if (myLink->getElement(eID))
        lastGroup->addElement(eID);
      else if (nNotes++ < oldNotes+10)
        ListUI <<"\n   * Note: Ignoring non-existing element "<< eID
               <<" in Nastran SET "<< SID;
    if (nNotes > oldNotes) ListUI <<"\n";
    lastGroup->sortElements();
    myLink->addGroup(lastGroup);
  }

  STOPP_TIMER("process_SET1")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// SPC ///
////////////////////////////////////////////////////////////////////////////////

FFlNastranReader::IntMap::iterator FFlNastranReader::setDofFlag (int n, int flg)
{
  IntMap::iterator it = nodeStat.find(n);
  if (it == nodeStat.end())
    return nodeStat.insert({n,flg}).first;

  if (flg > 0)
    it->second |= flg;
  else if (flg < 0)
    it->second = -(-it->second | -flg);

  return it;
}


bool FFlNastranReader::process_SPC (std::vector<std::string>& entry)
{
  START_TIMER("process_SPC")

  int node, dofs;
  double D = 0.0;

  for (size_t i = 1; i+2 < entry.size(); i+=3)
  {
    node = dofs = 0;
    D = 0.0;
    CONVERT_ENTRY ("SPC",
		   fieldValue(entry[i],node) &&
		   fieldValue(entry[i+1],dofs) &&
		   fieldValue(entry[i+2],D));

    if (node > 0 && dofs > 0)
    {
      this->setDofFlag(node,-convertDOF(dofs));
      if (D != 0.0)
      {
        nWarnings++;
        ListUI <<"\n  ** Warning: SPC "<< node <<" "<< dofs
               <<" has a non-zero prescribed value "<< D
               <<" (ignored).\n";
      }
    }
  }

  STOPP_TIMER("process_SPC")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// SPC1 //
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::process_SPC1 (std::vector<std::string>& entry)
{
  START_TIMER("process_SPC1")

  int node1 = 0, node2 = 0, dofs = 0;

  if (entry.size() > 1)
    CONVERT_ENTRY ("SPC1",fieldValue(entry[1],dofs));

  if (dofs > 0)
  {
    int status = convertDOF(dofs);
    for (size_t i = 2; i < entry.size(); i++)
      if (entry[i] == "THRU")
        node1 = node2;
      else
      {
        node2 = 0;
        CONVERT_ENTRY ("SPC1",fieldValue(entry[i],node2));
        if (node1 > 0)
        {
          for (int node = node1+1; node <= node2; node++)
            this->setDofFlag(node,-status);
          node1 = 0;
        }
        else if (node2 > 0)
          this->setDofFlag(node2,-status);
      }
  }

  STOPP_TIMER("process_SPC1")
  return true;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::insertBeamPropMat (const char* bulk, int PID, int MID)
{
  IntMap& beamMID = propMID[FFlTypeInfoSpec::BEAM_ELM];

  // Check if this beam property number has been assigned a material earlier
  IntMap::const_iterator pit = beamMID.find(PID);
  if (pit == beamMID.end())
    beamMID[PID] = MID;
  else if (pit->second != MID)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Beam property "<< bulk <<" "<< PID
           <<" uses material "<< MID
           <<".\n              However, this property number has already been"
           <<" assigned material number "<< pit->second
           <<".\n              All beam properties referring to the same"
           <<" material need unique identification numbers"
           <<".\n              You need to edit the FE data file unless the"
           <<" materials "<< MID <<" and "<< pit->second <<" are identical.\n";
    return false;
  }

  return true;
}

#ifdef FF_NAMESPACE
} // namespace
#endif
