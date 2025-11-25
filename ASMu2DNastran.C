// $Id$
//==============================================================================
//!
//! \file ASMu2DNastran.C
//!
//! \date Feb 27 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of 2D %Lagrange FE models from Nastran Bulk Data File (BDF).
//!
//==============================================================================

#include "ASMu2DNastran.h"
#include "FiniteElement.h"
#include "ElasticBeam.h"
#include "BeamProperty.h"
#include "ElementBlock.h"
#include "MPC.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "Profiler.h"
#include "IFEM.h"
#include <sstream>

namespace ASM {
  char useBeam = 1; //!< If nonzero, include beam elements as a separate patch
  char skipMass = 0; //!< 1: ignore one-noded mass elements, 2: no VTF output
  bool skipRBE2 = false; //!< If \e true, ignore all RBE2 elements
  bool replRBE3 = false; //!< If \e true, convert RBE3 elements to RBE2 elements
  bool skipLoad = false; //!< If \e true, ignore all FE surface loads
  bool readSets = true; //!< If \e true, read pre-bulk Nastran SET definitions
  double Ktra = 0.0; //!< Translation stiffness for added springs in slave nodes
  double Krot = 0.0; //!< Rotation stiffness for added sprins in slave nodes
  IntVec fixRBE3; //!< List of RBE3 elements to be constrained
  IntVec ignoredElms; //!< List of elements to ignore
}

#ifdef HAS_FFLLIB
#include "FFlLinkHandler.H"
#include "FFlNastranReader.H"
#include "FFlElementBase.H"
#include "FFlLoadBase.H"
#include "FFlGroup.H"
#include "FFlFEParts/FFlNode.H"
#include "FFlFEParts/FFlPMAT.H"
#include "FFlFEParts/FFlPMASS.H"
#include "FFlFEParts/FFlPBUSHCOEFF.H"
#include "FFlFEParts/FFlPTHICK.H"
#include "FFlFEParts/FFlPWAVGM.H"
#include "FFlFEParts/FFlPBEAMSECTION.H"
#include "FFlFEParts/FFlPBEAMECCENT.H"
#include "FFlFEParts/FFlPORIENT.H"
#include <unordered_map>

namespace FFlNastran {
  extern std::string mainPath;
}


/*!
  \brief Class for reading Nastran bulk data into a FFlLinkHandler object.
*/

class MyNastranReader : public FFlNastranReader
{
  int nPreBulk; //!< Number of lines before BEGIN BULK

public:
  //! \brief The constructor forwards to the parent class constructor.
  MyNastranReader(FFlLinkHandler& fePart, int lCount)
    : FFlNastranReader(&fePart,lCount), nPreBulk(lCount) {}

  //! \brief Reads the FE model from the Nastran file stream.
  bool readFE(std::istream& is, std::istream& iset)
  {
    PROFILE("Nastran file parser");

    if (!this->resolve(this->read(is)))
      myLink->deleteGeometry(); // Parsing failure, delete all FE data
    else if (nWarnings+nNotes > 0)
      IFEM::cout <<"\n  ** Parsing FE data succeeded."
                 <<"\n     However, "<< nWarnings
                 <<" warning(s) and "<< nNotes <<" note(s) were reported.\n"
                 <<"     Review the messages and check the FE data file.\n"
                 << std::endl;
    if (!myLink->hasGeometry())
      return false;

    // Remove all solid elements (not yet supported),
    // and (optionally) all the mass and beam elements, and ...
    ElementsVec toBeErased;
    ElementsCIter it;
    for (it = myLink->elementsBegin(); it != myLink->elementsEnd(); ++it)
      if ((*it)->getCathegory() == FFlTypeInfoSpec::SOLID_ELM ||
          utl::findIndex(ASM::ignoredElms,(*it)->getID()) >= 0 ||
          (ASM::skipRBE2 && (*it)->getTypeName() == "RGD") ||
          (ASM::skipMass == 1 && (*it)->getTypeName() == "CMASS") ||
          (!ASM::useBeam && (*it)->getCathegory() == FFlTypeInfoSpec::BEAM_ELM))
        toBeErased.push_back(*it);
    if (!toBeErased.empty())
    {
      IFEM::cout <<"  ** Erasing "<< toBeErased.size()
                 <<" elements from the model."<< std::endl;
      myLink->removeElements(toBeErased);
    }

    // Now parse the element set definitions, if any
    lastComment = { 0, "" };
    return iset ? this->processAllSets(iset,nPreBulk) : true;
  }
};
#endif


ASMu2DNastran::ASMu2DNastran (unsigned char n, unsigned char n_f,
                              const std::string& path)
  : ASMu2DLag(n,n_f,'N')
{
  massMax = 1.0;
  beamPatch = nullptr;
  nGnod.fill(0);
#ifdef HAS_FFLLIB
  FFlNastran::mainPath = path;
#endif
}


ASMu2DNastran::~ASMu2DNastran ()
{
  for (const ElementBlockID& geo : myBlocks)
    delete geo.second;
}


bool ASMu2DNastran::read (std::istream& is)
{
  this->clear(true);

  nGnod.fill(0);
  spiders.clear();
  myProps.clear();
  myBprops.clear();
  myStiff.clear();
  myMass.clear();
  myLoads.clear();
  elmProp.clear();
  elmPres.clear();
  IntVec beamElms;
  IntVec beamNodes;
  IntMat beamMNPC;
  int nBel = 0;
  int nErr = 0;
  int iErr = 0;

  // Fast-forward until "BEGIN BULK"
  int lCount = 0;
  char cline[256];
  std::stringstream sets;
  while (is.getline(cline,255))
    if (!strncmp(cline,"BEGIN BULK",10))
      break;
    else
    {
      ++lCount;
      // Copy element SET definitions to a second stream,
      // since they have to be parsed after the FE model is loaded
      if (ASM::readSets && (!strncmp(cline,"SET ",4) || sets.tellp() > 0))
        sets << cline << '\n';
    }

  if (!is) return false; // No bulk data file

  if (lCount > 0)
    IFEM::cout <<"\tNastran bulk data starting at line "
               << lCount+1 << std::endl;

#ifdef HAS_FFLLIB
  bool ok = false;
  FFlLinkHandler fem;
  if (MyNastranReader reader(fem,lCount); reader.readFE(is,sets))
  {
    PROFILE("Resolve cross references");
    ok = fem.resolve();
  }
  if (ok)
    IFEM::cout <<"\nParsing Nastran bulk data file succceeded."<< std::endl;
  else
  {
    std::cerr <<"\n *** Parsing/resolving FE data failed.\n"
              <<"     The FE model is probably not consistent and has not been"
              <<" resolved completely.\n";
    return false;
  }

  for (int eId : ASM::fixRBE3)
    if (FFlElementBase* elm = fem.getElement(eId); elm)
      if (FFlNode* refNode = elm->getNode(1); refNode && refNode->isRefNode())
        if ((refNode = fem.createAttachableNode(refNode,refNode->getPos(),
                                                nullptr,ASM::Ktra,ASM::Krot)))
        {
          refNode->setExternal();
          IFEM::cout <<"   * Adding External node "<< refNode->getID()
                     <<" connected to RBE3 reference node "<< elm->getNodeID(1)
                     << std::endl;
        }

  nnod = fem.getNodeCount(FFlLinkHandler::FFL_FEM);
  nel  = fem.getElementCount(FFlTypeInfoSpec::SHELL_ELM);
  nBel = fem.getElementCount(FFlTypeInfoSpec::BEAM_ELM);
#ifdef INT_DEBUG
  fem.dump();
#else
  IFEM::cout <<"\nTotal number of nodes:          "<< nnod
             <<"\nNumber of shell elements:       "<< nel;
  if (nBel > 0)
    IFEM::cout <<"\nNumber of beam elements:        "<< nBel;
  IFEM::cout <<"\nNumber of constraint elements:  "
             << fem.getElementCount(FFlTypeInfoSpec::CONSTRAINT_ELM)
             <<"\nNumber of other elements:       "
             << fem.getElementCount(FFlTypeInfoSpec::OTHER_ELM)
             << std::endl;
#endif
  if (size_t allN = fem.getNodeCount(FFlLinkHandler::FFL_ALL); allN > nnod)
    IFEM::cout <<"\n  ** Warning: This model contains "<< allN-nnod
               <<" node(s) without any element connections (ignored)."
               <<"\n     Please check the FE data file.\n"<< std::endl;

  if (ASMbase::modelSize == 1.0)
  {
    FaVec3 max, min;
    fem.getExtents(max,min);
    ASMbase::modelSize = (max-min).length();
  }
  IFEM::cout <<"Model extension (diameter):     "
             << ASMbase::modelSize << std::endl;

  PROFILE("Convert to IFEM");

  myMLGN.reserve(nnod);
  myMLGE.reserve(nel);
  myCoord.reserve(nnod);
  elmProp.reserve(nel);
  beamNodes.reserve(nBel);
  beamElms.reserve(nBel);
  beamMNPC.reserve(nBel);

  std::unordered_map<int,int> MEXINN; // External to internal node number map
  // Using an unordered_map instead of relying on ASMbase::getNodeIndex,
  // which is way to slow for very large models.

  // Extract the nodal points
  for (size_t inod = 1; inod <= nnod; inod++)
  {
    FFlNode* node = fem.getFENode(inod);
    int nid = node->getID();
    const FaVec3& X = node->getPos();
#if INT_DEBUG > 10
    std::cout << MLGN.size() <<" "<< nid <<": "<< X << std::endl;
#endif
    myMLGN.push_back(nid);
    myCoord.push_back(Vec3(X.x(),X.y(),X.z()));

    if (!MEXINN.emplace(nid,inod).second)
    {
      IFEM::cout <<"  ** Multiple instances of node number "<< nid << std::endl;
      iErr++;
    }

    if (node->isExternal()) // Create a node set for the supernodes
      this->addToNodeSet("ASET",MLGN.size());
    else if (node->isFixed())
    {
      // Create a node set for prescribed nodes,
      // one for each DOF constellation
      std::string cstat = std::to_string(-node->getStatus(-64));
      this->addToNodeSet("SPC"+cstat,MLGN.size());
    }
  }
  if (iErr > 0)
    std::cerr <<" *** A total of "<< iErr
              <<" multiple nodes where detected."<< std::endl;

  // Extract the element data
  for (ElementsCIter e = fem.elementsBegin(); e != fem.elementsEnd(); ++e)
  {
    int eid = (*e)->getID();

    // Nodal connectivities
    IntVec mnpc((*e)->getNodeCount());
    for (size_t i = 0; i < mnpc.size(); i++)
      mnpc[i] = MEXINN.at((*e)->getNodeID(1+i)) - 1;

    if ((*e)->getCathegory() == FFlTypeInfoSpec::BEAM_ELM && mnpc.size() == 2)
    {
#if INT_DEBUG > 1
      std::cout <<"Beam element "<< beamElms.size() <<" "<< eid <<":";
      for (int inod : mnpc) std::cout <<" "<< MLGN[inod];
#endif
      this->addBeamElement(*e,eid,mnpc,beamMNPC,beamElms,beamNodes,
                           nErr, ASM::useBeam == 1);
    }
    else if ((*e)->getCathegory() == FFlTypeInfoSpec::SHELL_ELM)
    {
      if (mnpc.size() == 4)
      {
        // Need to swap nodes 3 and 4 into IFEM order
        std::swap(mnpc[2],mnpc[3]);
        swapNode34 = true;
      }

#if INT_DEBUG > 10
      std::cout <<"Shell element "<< MLGE.size() <<" "<< eid <<":";
      for (int inod : mnpc) std::cout <<" "<< MLGN[inod];
#endif
      this->addShellElement(*e,eid,mnpc);
    }
    else if ((*e)->getTypeName() == "RGD" && mnpc.size() > 1)
    {
#if INT_DEBUG > 5
      std::cout <<"Rigid element "<< eid <<": master = "<< MLGN[mnpc.front()]
                <<" slaves =";
      for (size_t i = 1; i < mnpc.size(); i++) std::cout <<" "<< MLGN[mnpc[i]];
      std::cout << std::endl;
#endif
      spiders.push_back(mnpc);
      for (int& inod : mnpc) ++inod; // Need 1-based node indices
      this->addRigidCouplings((*e)->getNodeID(1),coord[mnpc.front()-1],
                              IntVec(mnpc.begin()+1,mnpc.end()));
    }
    else if ((*e)->getTypeName() == "WAVGM" && mnpc.size() > 1)
    {
      unsigned int refCount = (*e)->getNode(1)->getRefCount();
#ifdef FFL_REFCOUNT
      if (refCount < 2)
      {
        IFEM::cout <<"  ** Ignoring constraint element WAVGM "<< eid
                   <<" since it's reference node "<< MLGN[mnpc.front()]
                   <<" is not connected to any other element"<< std::endl;
	continue;
      }
#endif
#if INT_DEBUG > 1
      std::cout <<"Weighted average motion element "<< eid
                <<": reference node = "<< MLGN[mnpc.front()]
                <<" (ref.count = "<< refCount <<")";
#if INT_DEBUG > 5
      std::cout <<"\n\tmasters =";
      for (size_t i = 1; i < mnpc.size(); i++) std::cout <<" "<< MLGN[mnpc[i]];
#endif
      std::cout << std::endl;
#endif
      spiders.push_back(mnpc);
      if (ASM::replRBE3)
      {
        // Replace the flexible RBE3 element into an equivalent RBE2
        for (int& inod : mnpc) ++inod; // Need 1-based node indices
        this->addRigidCouplings((*e)->getNodeID(1),coord[mnpc.front()-1],
                                IntVec(mnpc.begin()+1,mnpc.end()));
      }
      else
        this->addFlexibleCouplings(*e,eid,mnpc);
    }
    else if ((*e)->getTypeName() == "CMASS" && !mnpc.empty())

      this->addMassElement(*e,eid,mnpc.front());

    else if ((*e)->getTypeName() == "BUSH" && mnpc.size() == 2)

      this->addSpringElement(*e,eid,mnpc,nErr);

    else
      IFEM::cout <<"  ** Ignored element "<< (*e)->getTypeName()
                 <<" "<< eid << std::endl;
  }
  if (nErr > 20)
    std::cerr <<" *** A total of "<< nErr
              <<" elements lack cross section properties."<< std::endl;

  // Extract the pressure loads, if any
  std::set<int> loadCases;
  if (!ASM::skipLoad)
    fem.getLoadCases(loadCases);
  for (int lc : loadCases)
  {
    LoadsVec loads;
    fem.getLoads(lc,loads);
    for (FFlLoadBase* load : loads)
    {
      std::vector<FaVec3> Pe;
      int iel, iface = 0;
      while ((iel = load->getLoad(Pe,iface)) > 0)
        if (iface >= 0)
        {
          Vec3Vec& elLoad = myLoads[iel];
          elLoad.clear();
          elLoad.reserve(Pe.size());
          for (const FaVec3& P : Pe)
            elLoad.push_back(P.getPt());
          if (elLoad.size() == 4)
            std::swap(elLoad[2],elLoad[3]);
#if INT_DEBUG > 5
          std::cout <<"Surface pressure on element "<< iel <<":";
          for (const Vec3& p : elLoad) std::cout <<"  "<< p;
          std::cout << std::endl;
#endif
        }
    }
  }

  if (nBel > 0 && nErr == 0)
  {
    // Create a separate patch for the beam elements
    Vec3Vec bXYZ;
    bXYZ.reserve(beamNodes.size());
    for (int node : beamNodes)
      bXYZ.push_back(coord[MEXINN.at(node)-1]);
    beamPatch = new ASMuBeam(bXYZ,beamMNPC,beamNodes,beamElms,myBprops,nsd,nf);
    IFEM::cout <<"\tCreated beam patch "<< beamPatch->idx+1
               <<" with "<< beamNodes.size() <<" nodes"<< std::endl;
  }

  // Extract the element sets, if any
  for (GroupCIter g = fem.groupsBegin(); g != fem.groupsEnd(); ++g)
  {
    std::string name = g->second->getName() + "_" + std::to_string(g->first);
    for (const GroupElemRef& elm : *g->second)
      if (utl::findIndex(beamElms,elm->getID()) < 0)
        this->addToElemSet(name,elm->getID(),true);
      else
        IFEM::cout <<"  ** Ignoring beam element "<< elm->getID()
                   <<" in element group "<< name << std::endl;
  }
#endif

  if (!nodeSets.empty())
  {
    IFEM::cout <<"Pre-defined node sets:          "<< nodeSets.size();
    for (const ASM::NodeSet& ns : nodeSets)
      IFEM::cout <<"\n\t\""<< ns.first <<"\"\t"<< ns.second.size() <<" nodes";
    IFEM::cout << std::endl;
  }

  if (!elemSets.empty())
  {
    IFEM::cout <<"Pre-defined element sets:       "<< elemSets.size();
    for (const ASM::NodeSet& eset : elemSets)
      IFEM::cout <<"\n\t\""<< eset.first <<"\"\t"<< eset.second.size()
                 <<" elements";
    IFEM::cout << std::endl;
  }

  if (this->empty() && nBel == 0)
  {
    std::cerr <<"\n *** No elements in this patch (ignored)."<< std::endl;
    return false;
  }

  return nErr == 0 && iErr == 0;
}


#ifdef HAS_FFLLIB
#define GET_ATTRIBUTE(el,att) dynamic_cast<FFl##att*>(el->getAttribute(#att))

void ASMu2DNastran::addBeamElement (FFlElementBase* elm, int eId,
                                    const IntVec& mnpc, IntMat& beamMNPC,
                                    IntVec& beamElms, IntVec& beamNodes,
                                    int& nErr, bool useEcc)
{
  beamElms.push_back(eId);
  beamMNPC.push_back(mnpc);
  for (int& inod : beamMNPC.back())
    if (int node = MLGN[inod]; (inod = utl::findIndex(beamNodes,node)) < 0)
    {
      beamNodes.push_back(node);
      inod = beamNodes.size() - 1;
    }

  BeamProps& bprop = myBprops[eId];
  if (FFlPMAT* mat = GET_ATTRIBUTE(elm,PMAT); mat)
  {
    bprop.Emod = mat->youngsModule.getValue();
    bprop.Gmod = mat->shearModule.getValue();
    bprop.Rho  = mat->materialDensity.getValue();
  }
  else
    IFEM::cout <<"  ** No material attached to element "<< eId
               <<", using default properties"<< std::endl;

  if (FFlPBEAMSECTION* bsec = GET_ATTRIBUTE(elm,PBEAMSECTION); bsec)
  {
    size_t i = 0;
    for (FFlFieldBase* field : *bsec)
      if (i < bprop.cs.size())
        bprop.cs[i++] = static_cast<FFlField<double>*>(field)->getValue();
    for (i = 4; i < 6; i++)  // Invert the shear reduction factors, since
      if (bprop.cs[i] > 0.0) // in FFlLib the factors As/A are assumed
        bprop.cs[i] = 1.0 / bprop.cs[i]; // but we need the A/As factors
  }
  else if (++nErr <= 20)
    std::cerr <<" *** No beam cross section attached to beam element "<< eId
              << std::endl;

  if (FFlPORIENT* bori = GET_ATTRIBUTE(elm,PORIENT); bori)
    bprop.Zaxis = Vec3(bori->directionVector.getValue().getPt());

  if (FFlPBEAMECCENT* bEcc = GET_ATTRIBUTE(elm,PBEAMECCENT); bEcc && useEcc)
    bprop.eccN = {
      -Vec3(bEcc->node1Offset.getValue().getPt()),
      -Vec3(bEcc->node2Offset.getValue().getPt())
    };

#if INT_DEBUG > 1
  std::cout <<" E="<< bprop.Emod <<" G="<< bprop.Gmod
            <<" rho="<< bprop.Rho <<"\nCS:";
  for (double cv : bprop.cs) std::cout <<" "<< cv;
  if (!bprop.Zaxis.isZero()) std::cout <<"\nZ-axis: "<< bprop.Zaxis;
  if (!bprop.eccN[0].isZero()) std::cout <<"\necc1: "<< bprop.eccN[0];
  if (!bprop.eccN[1].isZero()) std::cout <<"\necc2: "<< bprop.eccN[1];
  std::cout << std::endl;
#endif
}


void ASMu2DNastran::addShellElement (FFlElementBase* elm, int eId,
                                     const IntVec& mnpc)
{
  myMLGE.push_back(eId);
  myMNPC.push_back(mnpc);

  // Extract material properties
  ShellProps sprop;
  if (FFlPMAT* mat = GET_ATTRIBUTE(elm,PMAT); mat)
  {
    sprop.id1  = mat->getID();
    sprop.Emod = mat->youngsModule.getValue();
    sprop.Rny  = mat->poissonsRatio.getValue();
    sprop.Rho  = mat->materialDensity.getValue();
  }
  else
    IFEM::cout <<"  ** No material attached to element "<< eId
               <<", using default properties"<< std::endl;

  // Extract shell thickness
  if (FFlPTHICK* thk = GET_ATTRIBUTE(elm,PTHICK); thk)
  {
    sprop.id2   = thk->getID();
    sprop.Thick = thk->thickness.getValue();
  }
  else
    IFEM::cout <<"  ** No shell thickness attached to element "<< eId
               <<", using default value "<< sprop.Thick << std::endl;

  std::vector<ShellProps>::const_iterator it;
  if ((it = std::find(myProps.begin(),myProps.end(),sprop)) == myProps.end())
  {
    elmProp.push_back(myProps.size());
    myProps.emplace_back(sprop);
#if INT_DEBUG > 10
    std::cout <<" mid="<< sprop.id1 <<" tid="<< sprop.id2
              <<" E="<< sprop.Emod  <<" nu="<< sprop.Rny <<" rho="<< sprop.Rho
              <<" t="<< sprop.Thick << std::endl;
#endif
  }
  else
  {
    elmProp.push_back(it - myProps.begin());
#if INT_DEBUG > 10
    std::cout <<" mid="<< sprop.id1 <<" tid="<< sprop.id2 << std::endl;
#endif
  }
}


void ASMu2DNastran::addMassElement (FFlElementBase* elm, int eId, int inod)
{
  FFlPMASS* mass = GET_ATTRIBUTE(elm,PMASS);
  if (!mass)
  {
    IFEM::cout <<"  ** No mass property attached to mass element "<< eId
               <<" (ignored)"<< std::endl;
    return;
  }
#if INT_DEBUG > 5
  std::cout <<"Mass element "<< MLGE.size() <<" "<< eId
            <<": node = "<< MLGN[inod] << std::endl;
#endif

  const RealArray& Mvec = mass->M.getValue();
  RealArray::const_iterator m = Mvec.begin();

  if (myMass.empty() || *m > massMax)
    massMax = *m;

  elmProp.push_back(0);
  myMLGE.push_back(eId);
  myMNPC.push_back({inod});
  Matrix& M = myMass[eId];
  M.resize(6,6);
  for (int i = 1; i <= 6; i++)
    for (int j = 1; j <= i && m != Mvec.end(); j++)
    {
      M(i,j) = *(m++);
      if (j < i) M(j,i) = M(i,j);
    }
}


void ASMu2DNastran::addSpringElement (FFlElementBase* elm, int eId,
                                      const IntVec& mnpc, int& nErr)
{
  FFlPBUSHCOEFF* bush = GET_ATTRIBUTE(elm,PBUSHCOEFF);
  if (!bush)
  {
    IFEM::cout <<"  ** No property attached to bush element "<< eId
               <<" (ignored)"<< std::endl;
    return;
  }
#if INT_DEBUG > 5
  std::cout <<"\nBush element "<< MLGE.size() <<" "<< eId
            <<": nodes = "<< MLGN[mnpc.front()] <<" "<< MLGN[mnpc.back()];
#endif

  elmProp.push_back(0);
  myMLGE.push_back(eId);
  myMNPC.push_back(mnpc);
  Matrix& K = myStiff[eId];
  K.resize(12,12,true);
  for (int i = 1; i <= 6; i++)
  {
    K(i,i) = K(6+i,6+i) = bush->K[i-1].getValue();
    K(i,6+i) = K(6+i,i) = -K(i,i);
  }

  // Transform to global coordinate axes
  double X[3], Y[3], Z[3], T[9];
  if (elm->getNodalCoor(X,Y,Z) >= 0 && elm->getLocalSystem(T))
  {
    Vec3 ecc1(X[1]-X[0], Y[1]-Y[0], Z[1]-Z[0]);
    Vec3 ecc2(X[2]-X[0], Y[2]-Y[0], Z[2]-Z[0]);
    Matrix Tlg(3,3); Tlg.fill(T);
    // TODO: We need to use Tlg.transpose() here instead,
    // to get the same transformation matrix as in FEDEM.
    // But not sure FEDEM is right there.
    if (utl::transform(K,Tlg))
      ElasticBeam::transform(K,ecc1,ecc2);
#if INT_DEBUG > 5
    std::cout <<"\nTransformation matrix:"<< Tlg
              <<"Ecc1: "<< ecc1 <<"\nEcc2: "<< ecc2;
#endif
  }
  else if (++nErr <= 20)
    std::cerr <<" *** No local axes for bush element "<< eId << std::endl;

#if INT_DEBUG > 5
  std::cout <<"\nStiffness matrix:"<< K;
#endif
}


void ASMu2DNastran::addFlexibleCouplings (FFlElementBase* elm, int eId,
                                          const IntVec& mnpc)
{
  int indC[6] = { -1, -1, -1, 0, 0, 0 };
  std::set<int> refC;
  RealArray weights;
  if (FFlPWAVGM* wavgm = GET_ATTRIBUTE(elm,PWAVGM); wavgm)
  {
    int div = 100000;
    int dofIds = wavgm->refC.getValue();
    if (dofIds > 0)
      for (int i = 0; i < 6; i++)
      {
        int dof = dofIds/div;
        if (dof > 0 && dof < 7) refC.insert(dof);
        dofIds -= dof*div;
        div    /= 10;
      }
    for (int i = 0; i < 6; i++)
      indC[i] = wavgm->indC[i].getValue();
    weights = wavgm->weightMatrix.getValue();
  }
  else
    for (int lDof = 1; lDof <= 6; lDof++)
      refC.insert(lDof);

  size_t icol = 0;
  Matrix Xnod(nsd,mnpc.size());
  for (int inod : mnpc)
    Xnod.fillColumn(++icol,coord[inod].ptr());

  double* work = new double[10*mnpc.size()-7];
  for (int lDof : refC)
    this->addFlexibleCoupling(eId,lDof,indC,weights,work,mnpc,Xnod);
  delete[] work;
}
#endif


bool ASMu2DNastran::getMassProp (int eId, size_t igel,
                                 double& rho, double& t) const
{
  size_t iprop = igel < firstEl+elmProp.size() ? elmProp[igel-firstEl] : 0;
  if (iprop > myProps.size())
  {
    std::cerr <<" *** No properties for shell element "<< eId << std::endl;
    return false;
  }

  rho = myProps[iprop].Rho;
  t   = myProps[iprop].Thick;

  return true;
}


bool ASMu2DNastran::getStiffProp (int eId, size_t igel,
                                  double& E, double& nu) const
{
  size_t iprop = igel < firstEl+elmProp.size() ? elmProp[igel-firstEl] : 0;
  if (iprop > myProps.size())
  {
    std::cerr <<" *** No properties for shell element "<< eId << std::endl;
    return false;
  }

  E  = myProps[iprop].Emod;
  nu = myProps[iprop].Rny;

  return true;
}


bool ASMu2DNastran::getThickness (size_t iel, double& t) const
{
  int eId = this->getElmID(iel);
  if (eId < 1) return false;

  if (myStiff.find(eId) != myStiff.end())
    t = 0.0; // Silently ignore spring elements
  else if (myMass.find(eId) != myMass.end())
    t = 0.0; // Silently ignore mass elements
  else
  {
    size_t iprop = iel < elmProp.size() ? elmProp[iel] : 0;
    if (iprop > myProps.size())
    {
      std::cerr <<" *** No properties for shell element "<< eId << std::endl;
      return false;
    }
    t = myProps[iprop].Thick;
  }

  return true;
}


bool ASMu2DNastran::getStiffnessMatrix (int eId, Matrix& K) const
{
  std::map<int,Matrix>::const_iterator it = myStiff.find(eId);
  if (it == myStiff.end())
  {
    std::cerr <<" *** No stiffness matrix for element "<< eId << std::endl;
    return false;
  }

  K = it->second;

  return true;
}


bool ASMu2DNastran::getMassMatrix (int eId, Matrix& M) const
{
  std::map<int,Matrix>::const_iterator it = myMass.find(eId);
  if (it == myMass.end())
  {
    std::cerr <<" *** No mass matrix for element "<< eId << std::endl;
    return false;
  }

  M = it->second;

  return true;
}


bool ASMu2DNastran::getLoadVector (int eId, const Vec3& g, Vector& S) const
{
  std::map<int,Matrix>::const_iterator it = myMass.find(eId);
  if (it == myMass.end())
  {
    std::cerr <<" *** No mass matrix for element "<< eId << std::endl;
    return false;
  }

  for (size_t i = 1; i <= 3 && i < it->second.rows() && i <= S.size(); i++)
    S(i) = it->second(i,1)*g.x + it->second(i,2)*g.y + it->second(i,3)*g.z;

  return true;
}


bool ASMu2DNastran::addPressureAt (Vec3& p, int eId, const RealArray& N) const
{
  std::map<int,Vec3Vec>::const_iterator it = myLoads.find(eId);
  if (it == myLoads.end() || it->second.empty()) return false;

  size_t n = it->second.size();
  if (n == 1 || N.size() > n)
    p += it->second.front(); // constant pressure
  else if (N.size() == 1 && N.front() < 0.0)
  {
    // N[0] contains the nodal index
    size_t inod = -N.front();
    if (inod > 0 && inod <= n)
      p += it->second[inod-1];
  }
  else // interpolate the nodal values
    for (size_t i = 0; i < N.size(); i++)
      p += N[i]*it->second[i];

  return true;
}


bool ASMu2DNastran::initPressureCache ()
{
  bool ok = true;
  size_t nPressElms = 0;
  elmPres.resize(nel,0);
  for (int iSet = 1; true; iSet++)
    if (const IntVec& eSet = this->getElementSet(iSet); eSet.empty())
      break;
    else for (int iel : eSet)
      if (iel > 0 && iel <= static_cast<int>(nel))
      {
        if (!elmPres[iel-1])
        {
          ++nPressElms;
          elmPres[iel-1] = iSet;
        }
        else
        {
          ok = false;
          IFEM::cout <<"  ** Element "<< iel <<" ("<< MLGE[iel-1]
                     <<") is present in more than one element set.\n    "
                     <<" Can't use element index cache for pressure lookup."
                     << std::endl;
        }
      }

  if (!ok || !nPressElms)
    elmPres.clear();
  else
    IFEM::cout <<"\nTotal element set size: "<< nPressElms << std::endl;

  return ok;
}


/*!
  This method is used to speed up the check for whether an element is associated
  with the specified element set for pressure load or not. For very large models
  the ASMu2DLag::isInElementSet() method is too slow as it involves a search in
  a vector of element numbers. Instead, the index cache \ref elmPres is used,
  which for a given element index tells which element set (if any) it is
  associated with. The limitation is that each element can only be associated
  with one such set. If that does not hold, the slow search has to be used.
*/

bool ASMu2DNastran::checkPressSet (int iEl, size_t idx, int iSet) const
{
  if (idx >= firstEl && idx < firstEl+elmPres.size())
    return static_cast<int>(elmPres[idx-firstEl]) == iSet;

  // No element pressure index cache, fall-back for slow search
  return this->isInElementSet(iSet,-iEl);
}


bool ASMu2DNastran::getShellNormals (std::vector<Vec3Pair>& normals) const
{
  normals.reserve(nel);
  for (size_t iel = 0; iel < nel; iel++)
    if (const IntVec& mnpc = MNPC[iel]; mnpc.size() == 3 || mnpc.size() == 4)
    {
      const size_t n1 = mnpc.size() - 3;
      const size_t n3 = mnpc.size() == 4 ? 3 : 1;
      Vec3 V3(coord[mnpc[n3]]-coord[mnpc[0]], coord[mnpc[2]]-coord[mnpc[n1]]);
      V3.normalize();
      normals.emplace_back(this->getElementCenter(1+iel),V3);
    }

  return !normals.empty();
}


#ifdef HAS_ANDES
extern "C" void ifem_elaxes_(const int& nelnod, const double* Xnod,
                             double* eX, double* eY, double* eZ, int& ierr);
#endif

/*!
  This method calculates and prints out the globalized shell surface axes.
*/

void ASMu2DNastran::printElmInfo (int iel, const IntegrandBase* integr) const
{
  if (integr && integr->getNoFields(2) <= 2) return;

#ifdef HAS_ANDES
  if (Matrix Xnod; this->getElementCoordinates(Xnod,iel))
  {
    double tmp[9];
    int ierr, nelnod = Xnod.cols();
    ifem_elaxes_(nelnod,Xnod.ptr(),tmp,tmp+3,tmp+6,ierr);
    if (ierr == 0)
      IFEM::cout <<"\teX = "<< Vec3(tmp) <<"\teY = "<< Vec3(tmp+3) << std::endl;
  }
#else
  IFEM::cout <<"\tXc = "<< this->getElementCenter(iel) << std::endl;
#endif
}


#ifdef HAS_ANDES
extern "C" void wavgmconstreqn_(const int& iel, const int& lDof,
                                const int& nM, const int& nW, const int* indC,
                                const double* tenc, const double* weight,
                                const double& epsX, double* dX,
                                double* work, double* omega,
                                const int& ipsw, const int& lpu);
#endif


void ASMu2DNastran::addFlexibleCoupling (int eId, int lDof, const int* indC,
                                         const RealArray& weights, double* work,
                                         const IntVec& mnpc, const Matrix& Xnod)
{
#ifdef HAS_ANDES
  const int nM = Xnod.cols() - 1;
  const int nW = weights.size();
#if INT_DEBUG > 2
  const int ips = INT_DEBUG;
#else
  const int ips = 0;
#endif
  const int lpu = 6;
  const double epsX = 1.0e-4;
  const double Zero = 1.0e-8;
  double* dX = work + nM;
  double* omega = dX + 3*(1+nM);
  wavgmconstreqn_(eId,lDof,nM,nW,indC,Xnod.ptr(),weights.data(),epsX,
                  dX,work,omega,ips,lpu);

  MPC* cons = new MPC(MLGN[mnpc.front()],lDof);
  if (this->addMPC(cons) && cons)
    for (int iM = 1; iM <= nM; iM++)
      for (int mDof = 1; mDof <= 6; mDof++, omega++)
        if (fabs(*omega) > Zero)
          cons->addMaster(MLGN[mnpc[iM]],mDof,*omega);
#if INT_DEBUG > 1
        else
          std::cout <<"  ** Ignoring small coupling coefficient "<< *omega
                    <<" to local dof "<< mDof
                    <<" of master node "<< MLGN[mnpc[iM]]
                    <<" in RBE3 element "<< eId << std::endl;
#endif
#endif
}


void ASMu2DNastran::addBlock (int idx, ElementBlock* blk)
{
  myBlocks.emplace_back(idx,blk);
  nGnod[idx < 0 ? 1 : 0] += blk->getNoNodes();
}


/*!
  This method creates an element block visualizing the point mass elements
  as spheres with radii proportional to the mass magnintudes.
*/

ElementBlock* ASMu2DNastran::immersedGeometry (char* name) const
{
  ElementBlock* geo = nullptr;
  if (myMass.empty() || ASM::skipMass) return geo;

  // Let the largest point mass be visualized as a sphere
  // with diameter ~5% of the total length of the structure
  const double massScale = 0.025*modelSize/massMax;

  for (const std::pair<const int,Matrix>& mass : myMass)
  {
    int inod = MNPC[this->getElmIndex(mass.first)-1].front();
    double R = mass.second(1,1)*massScale;
    ElementBlock* newBlock = new SphereBlock(coord[inod],R,8,6);
    if (!geo)
      geo = new ElementBlock(*newBlock);
    else
      geo->merge(*newBlock,false);

    const_cast<ASMu2DNastran*>(this)->addBlock(inod,newBlock);
  }

  if (name)
    sprintf(name,"Point masses for Patch %zu",idx+1);

  return geo;
}


/*!
  This method creates an element block visualizing constraint elements as a set
  of line elements with one common node. The other node indices are stored as
  external element ID of each line element in the element block.
  This is then used to assign correct deformation values in extraSolution().
*/

ElementBlock* ASMu2DNastran::extraGeometry (char* name) const
{
  ElementBlock* geo = nullptr;
  if (spiders.empty()) return geo;

  geo = new ElementBlock(2);
  geo->unStructResize(0,spiders.size());

  size_t iref = 0;
  for (const IntVec& mnpc : spiders)
    geo->setCoor(iref++,coord[mnpc.front()]);

  iref = 0;
  for (const IntVec& mnpc : spiders)
  {
    const Vec3& X = coord[mnpc.front()];
    for (size_t i = 1; i < mnpc.size(); i++)
      if (!X.equal(coord[mnpc[i]],1.0e-4*modelSize))
        geo->addLine(iref,coord[mnpc[i]],mnpc[i]);
    iref++;
  }

  if (geo->getNoElms() == 0)
  {
    delete geo;
    return nullptr;
  }

  const_cast<ASMu2DNastran*>(this)->addBlock(-1, new ElementBlock(*geo));

  if (name)
    sprintf(name,"Couplings for Patch %zu",idx+1);

  return geo;
}


ElementBlock* ASMu2DNastran::sensorGeometry (int idx, bool nodal) const
{
  if (idx < 1 || idx > static_cast<int>(nodal ? nnod : nel))
    return nullptr;

  Vec3 XYZ = nodal ? coord[idx-1] : this->getElementCenter(idx);
  IFEM::cout <<"   * Sensor "<< idx <<" at "<< XYZ << std::endl;
  double L = ASMbase::modelSize/400.0;
  return new CubeBlock(XYZ,L);
}


bool ASMu2DNastran::evalSolution (Matrix& sField, const IntegrandBase& integr,
                                  const int*, char) const
{
  return this->evalSecSolution(sField,integr,true);
}


bool ASMu2DNastran::evalSolution (Matrix& sField, const IntegrandBase& integr,
                                  const RealArray*, bool atElmCenters) const
{
  return this->evalSecSolution(sField,integr,!atElmCenters);
}


bool ASMu2DNastran::evalSolution (Matrix& sField, const IntegrandBase& integr,
                                  const IntVec& elements) const
{
  return this->evalSecSolution(sField,integr,false,elements);
}


/*!
  This method evaluates the secondary solution at either the element nodes
  and then performs nodal averaging to obtain the unique nodal values,
  or perform direct evaluation at the element centers.
  It is assumed that all calculations are performed by the
  IntegrandBase::evalSol() call, therefore no basis function evaluations here.
*/

bool ASMu2DNastran::evalSecSolution (Matrix& sField,
                                     const IntegrandBase& integr, bool atNodes,
                                     const IntVec& elements) const
{
  sField.clear();

  FiniteElement fe;
  Vector        solPt;
  Vectors       globSolPt(atNodes ? nnod : 0);
  IntVec        checkPt(atNodes ? nnod : 0,0);

  // Number of evaluation points
  size_t npt = atNodes ? nnod : (elements.empty() ? nel : elements.size());
  size_t ipt = 0;

  // Evaluate the secondary solution field at each element node or center
  for (size_t iel = 1; iel <= nel; iel++)
    if ((fe.iel = MLGE[iel-1]) > 0) // ignore the zero-area elements
    {
      if (!atNodes)
      {
        if (elements.empty())
        {
          // Consider 3- and 4-noded elements only
          if (MNPC[iel-1].size() < 3)
            continue;
        }
        else if (utl::findIndex(elements,iel) < 0)
          continue;
      }

      if (!this->getElementCoordinates(fe.Xn,iel))
        return false;

      const IntVec& mnpc = MNPC[iel-1];
      const size_t nenod = atNodes ? mnpc.size() : 1;
      for (size_t loc = 0; loc < nenod; loc++)
      {
        if (nenod == 3)
          switch (1+loc) {
          case 1: fe.xi = 1.0; fe.eta = 0.0; break;
          case 2: fe.xi = 0.0; fe.eta = 1.0; break;
          case 3: fe.xi = 0.0; fe.eta = 0.0; break;
          }
        else if (nenod == 4)
        {
          fe.xi  = -1.0 + 2.0*(loc%2);
          fe.eta = -1.0 + 2.0*(loc/2);
        }
        else if (mnpc.size() == 3)
          fe.xi = fe.eta = 1.0/3.0;
        else
          fe.xi = fe.eta = 0.0;

        if (!integr.evalSol(solPt,fe,Vec3(),mnpc))
          return false;
        else if (solPt.empty())
          break; // a valid element with no secondary solution

        if (sField.empty())
          sField.resize(solPt.size(),npt,true);

        if (!atNodes)
          sField.fillColumn(++ipt,solPt);
        else if (++checkPt[mnpc[loc]] == 1)
          globSolPt[mnpc[loc]] = solPt;
        else
          globSolPt[mnpc[loc]] += solPt;
      }
    }

  // Nodal averaging
  if (!atNodes)
    sField.resize(sField.rows(),ipt);
  for (size_t i = 0; i < checkPt.size(); i++)
    if (checkPt[i])
      sField.fillColumn(1+i, globSolPt[i] /= static_cast<double>(checkPt[i]));

  return true;
}


bool ASMu2DNastran::immersedSolution (Matrix& field, const Vector& locSol) const
{
  field.resize(3,nGnod.front());
  if (nGnod.front() < 1) return false;

  size_t inod = 0;
  for (const ElementBlockID& geo : myBlocks)
    if (geo.first >= 0)
    {
      unsigned int ipnod = nf*geo.first;
      Vec3 X0(coord[geo.first]);
      Vec3 U0(locSol.ptr() + ipnod);
      Tensor T0(locSol[ipnod+3],locSol[ipnod+4],locSol[ipnod+5]);
      for (size_t i = 0; i < geo.second->getNoNodes(); i++)
      {
        const Vec3& X1 = geo.second->getCoord(i);
        Vec3 U1 = (X0+U0) + T0*(X1-X0) - X1;
        field.fillColumn(++inod, U1.ptr());
      }
    }

  return true;
}


bool ASMu2DNastran::extraSolution (Matrix& field, const Vector& locSol) const
{
  field.resize(3,nGnod.back());
  if (nGnod.back() < 1) return false;

  size_t inod = 0;
  for (const IntVec& mnpc : spiders)
    field.fillColumn(++inod, locSol.ptr() + nf*mnpc.front());

  for (const ElementBlockID& geo : myBlocks)
    if (geo.first < 0)
      for (size_t i = 1; i <= geo.second->getNoElms(); i++)
        field.fillColumn(++inod, locSol.ptr() + nf*geo.second->getElmId(i));

  return true;
}


ASMuBeam::ASMuBeam (const Vec3Vec& nodes, const IntMat& mmnpc,
                    const IntVec& mlgn, const IntVec& mlge,
                    const std::map<int,BeamProps>& props,
                    unsigned char n, unsigned char n_f)
  : ASMu1DLag(n,n_f), myProps(props)
{
  myCoord = nodes;
  myMNPC  = mmnpc;
  myMLGN  = mlgn;
  myMLGE  = mlge;
}


bool ASMuBeam::getProps (int eId, double& E, double& G,
                         double& rho, BeamProperty& bprop) const
{
  std::map<int,BeamProps>::const_iterator it = myProps.find(eId);
  if (it == myProps.end())
  {
    std::cerr <<" *** No properties for beam element "<< eId << std::endl;
    return false;
  }

  E   = it->second.Emod;
  G   = it->second.Gmod;
  rho = it->second.Rho;
  bprop.setConstant(RealArray(it->second.cs.begin(),it->second.cs.end()));
  bprop.setEccentric(it->second.eccN.front(),it->second.eccN.back());
  return true;
}


/*!
  This method overrides the parent class method to account for possible
  element-level \a Zaxis properties, and a rotation angle between the
  principal axes of inertia for the beam cross section relative to
  the local element axes.
*/

bool ASMuBeam::initLocalElementAxes (const Vec3& Zaxis)
{
  size_t iel = 0;
  for (Tensor& Tlg : myCS)
  {
    int eId = MLGE[iel++];
    std::map<int,BeamProps>::const_iterator it = myProps.find(eId);
    if (it == myProps.end())
    {
      std::cerr <<" *** No properties for beam element "<< eId << std::endl;
      return false;
    }

    // Set up the global-to-local transformation matrix
    int n1 = MNPC[iel-1].front();
    int n2 = MNPC[iel-1].back();
    const Vec3& X1 = coord[n1];
    const Vec3& X2 = coord[n2];
    if (!it->second.Zaxis.isZero())
      Tlg = Tensor(X2-X1,it->second.Zaxis,false,true);
    else if (!Zaxis.isZero())
      Tlg = Tensor(X2-X1,Zaxis,false,true);
    else
      Tlg = Tensor(X2-X1,true);

    // Rotate from local element axes to principal axes
    const double phi = it->second.cs.back();
    if (fabs(phi) > 1.0e-6)
      Tlg.rotate(phi*M_PI/180.0,1);

#ifdef INT_DEBUG
    std::cout <<"\nLocal axes for beam element "<< eId
              <<",\nfrom  N"<< MLGN[n1] <<" = "<< X1
              <<"  to  N"<< MLGN[n2] <<" = "<< X2;
    if (!it->second.Zaxis.isZero())
      std::cout <<"  with Z-axis = "<< it->second.Zaxis;
    if (fabs(phi) > 1.0e-6)
      std::cout <<"  and phi = "<< phi;
    std::cout <<":\n";
    Tlg.print(std::cout,15);
#endif
  }

  return true;
}
