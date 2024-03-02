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
#ifdef HAS_FFLLIB
#include "FFlLinkHandler.H"
#include "FFlNastranReader.H"
#include "FFlElementBase.H"
#include "FFlFEParts/FFlNode.H"
#include "FFlFEParts/FFlPMAT.H"
#include "FFlFEParts/FFlPTHICK.H"


class MyNastranReader : public FFlNastranReader
{
public:
  MyNastranReader(FFlLinkHandler& fePart, int lCount)
    : FFlNastranReader(&fePart,lCount) {}
  virtual ~MyNastranReader() {};

  bool readAndResolve(std::istream& is)
  {
    if (!this->resolve(this->read(is)))
      myLink->deleteGeometry(); // Parsing failure, delete all FE data
    else if (nWarnings+nNotes > 0)
      std::cout <<"\n  ** Parsing FE data succeeded."
                <<"\n     However, "<< nWarnings
                <<" warning(s) and "<< nNotes <<" note(s) were reported.\n"
                <<"     Review the messages and check the FE data file.\n";
    return myLink->hasGeometry();
  }
};
#endif


bool ASMu2DNastran::read (std::istream& is)
{
  myMLGE.clear();
  myMLGN.clear();
  myCoord.clear();

  // Fast-forward until "BEGIN BULK"
  int lCount = 0;
  char cline[256];
  while (is.getline(cline,255))
    if (strncmp(cline,"BEGIN BULK",10))
      ++lCount;
    else
      break;

  if (!is) return false; // No bulk data file

#ifdef HAS_FFLLIB
#define GET_ATTRIBUTE(el,att) dynamic_cast<FFl##att*>((*el)->getAttribute(#att))

  FFlLinkHandler  fem;
  MyNastranReader reader(fem,lCount);
  if (!reader.readAndResolve(is) || !fem.resolve())
  {
    std::cout <<"\n *** Parsing/resolving FE data failed.\n"
              <<"     The FE model is probably not consistent and has not been"
              <<" resolved completely.\n";
    return false;
  }

  nnod = fem.getNodeCount();
  nel  = fem.getElementCount(FFlTypeInfoSpec::SHELL_ELM);
#ifdef INT_DEBUG
  fem.dump();
#else
  std::cout <<"\nTotal number of nodes:          "<< nnod
            <<"\nNumber of shell elements:       "<< nel
            <<"\nNumber of constraint elements:  "
            << fem.getElementCount(FFlTypeInfoSpec::CONSTRAINT_ELM)
            <<"\nNumber of other elements:       "
            << fem.getElementCount(FFlTypeInfoSpec::OTHER_ELM)
            << std::endl;
#endif

  myMLGN.reserve(nnod);
  myMLGE.reserve(nel);
  myCoord.reserve(nnod);

  // Extract the nodal points
  for (NodesCIter n = fem.nodesBegin(); n != fem.nodesEnd(); ++n)
  {
    int nid = (*n)->getID();
    const FaVec3& X = (*n)->getPos();
#if INT_DEBUG > 1
    std::cout << myMLGN.size() <<" "<< nid <<": "<< X << std::endl;
#endif
    myMLGN.push_back(nid);
    myCoord.push_back(Vec3(X.x(),X.y(),X.z()));

    if ((*n)->isExternal()) // Create a node set for the supernodes
      this->getNodeSet("ASET",lCount).push_back(myMLGN.size());
  }

  // Extract the element data
  for (ElementsCIter e = fem.elementsBegin(); e != fem.elementsEnd(); ++e)
  {
    int eid = (*e)->getID();

    // Nodal connectivities
    IntVec mnpc((*e)->getNodeCount());
    for (size_t i = 0; i < mnpc.size(); i++)
      mnpc[i] = this->getNodeIndex((*e)->getNodeID(1+i)) - 1;

    if ((*e)->getCathegory() == FFlTypeInfoSpec::SHELL_ELM)
    {
#if INT_DEBUG > 1
      std::cout <<"Shell element "<< myMLGE.size() <<" "<< eid <<":";
      for (int node : mnpc) std::cout <<" "<< MLGN[node];
#endif
      myMLGE.push_back(eid);
      myMNPC.push_back(mnpc);

      // Extract material properties
      ShellProps& sprop = myProps[eid];
      FFlPMAT* mat = GET_ATTRIBUTE(e,PMAT);
      if (mat)
      {
        sprop.Emod = mat->youngsModule.getValue();
        sprop.Rny  = mat->poissonsRatio.getValue();
        sprop.Rho  = mat->materialDensity.getValue();
      }
      else
        std::cout <<"  ** No material attached to element "<< eid
                  <<", using default properties"<< std::endl;

      // Extract shell thickness
      FFlPTHICK* thk = GET_ATTRIBUTE(e,PTHICK);
      if (thk)
        sprop.Thick = thk->thickness.getValue();
      else
        std::cout <<"  ** No shell thickness attached to element "<< eid
                  <<", using default value "<< sprop.Thick << std::endl;

#if INT_DEBUG > 1
      std::cout <<" t="<< sprop.Thick <<" E="<< sprop.Emod <<" nu="<< sprop.Rny
                <<" rho="<< sprop.Rho << std::endl;
#endif
    }
    else if ((*e)->getTypeName() == "RGD" && mnpc.size() > 1)
    {
#if INT_DEBUG > 1
      std::cout <<"Rigid element "<< eid <<": master = "<< MLGN[mnpc.front()]
                <<" slaves =";
      for (size_t i = 1; i < mnpc.size(); i++) std::cout <<" "<< MLGN[mnpc[i]];
      std::cout << std::endl;
#endif
      for (int& n : mnpc) ++n; // Need 1-based node indices
      this->addRigidCouplings((*e)->getNodeID(1),myCoord[mnpc.front()-1],
                              IntVec(mnpc.begin()+1,mnpc.end()));
    }
    else
      std::cout <<"  ** Ignored element "<< (*e)->getTypeName()
                <<" "<< eid << std::endl;
  }
#endif

  return true;
}


bool ASMu2DNastran::getProps (int eId, double& E, double& nu,
                              double& rho, double& t) const
{
  std::map<int,ShellProps>::const_iterator it = myProps.find(eId);
  if (it == myProps.end())
  {
    std::cerr <<" *** No properties for shell element "<< eId << std::endl;
    return false;
  }

  E   = it->second.Emod;
  nu  = it->second.Rny;
  rho = it->second.Rho;
  t   = it->second.Thick;

  return true;
}
