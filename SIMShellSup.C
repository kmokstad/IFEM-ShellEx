// $Id$
//==============================================================================
//!
//! \file SIMShellSup.C
//!
//! \date May 24 2025
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for linear elastic superelement FEM analysis.
//!
//==============================================================================

#include "SIMShellSup.h"
#include "SIMAndesShell.h"
#include "SAM.h"
#include "ASMsupel.h"
#include "HasGravityBase.h"
#include "ElementBlock.h"
#include "Utilities.h"
#include "Tensor.h"
#include "VTF.h"
#include "IFEM.h"
#include "tinyxml2.h"


SIMShellSup::SIMShellSup (const char* hd, bool fd) : SIMsupel(hd,-6), fixDup(fd)
{
  // Integrand class for superelement analysis with gravity loads
  class MyProblem : public HasGravityBase
  {
  public:
    MyProblem() : HasGravityBase(3) { npv = 6; }
    void printLog() const { IFEM::cout <<"Linear superelement"<< std::endl; }
  };

  myProblem = new MyProblem();
}


SIMShellSup::~SIMShellSup()
{
  for (std::pair<const std::string,FEmodel>& sub : mySubSim)
  {
    delete sub.second.sim;
    delete sub.second.blk;
  }
}


bool SIMShellSup::parse (const tinyxml2::XMLElement* elem)
{
  if (strcasecmp(elem->Value(),"elasticity"))
    return this->SIMsupel::parse(elem);

  std::string supId;
  SIMgeneric* supEl = nullptr;
  if (utl::getAttribute(elem,"supId",supId))
  {
    // This <elasticity> block defines the FE model of a superelement.
    // Create a separate SIM-object for it here for the recovery operation, but
    // with no integrand since recovery is based on an externally reduced model.
    IFEM::cout <<"  Parsing FE model for superelement \""<< supId <<"\"\n";
    supEl = new SIMAndesShell(-1);
  }

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(),"supernodes") && supEl)
    {
      std::string setName;
      IFEM::cout <<"  Parsing <supernodes> for \""<< supId <<"\".\n";
      std::vector<int>& superNodes = mySubSim[supId].superNodes;
      if (utl::getAttribute(child,"set",setName))
	superNodes = supEl->getNodeSet(setName);
      if (superNodes.empty())
        IFEM::cout <<"  ** Empty or non-existing node set \""<< setName <<"\".";
      else
      {
        IFEM::cout <<"\t"<< superNodes.front();
        for (size_t i = 1; i < 10 && i < superNodes.size(); i++)
          IFEM::cout <<" "<< superNodes[i];
        if (superNodes.size() > 10)
          IFEM::cout <<" ... ("<< superNodes.size() <<" nodes).";
      }
      IFEM::cout << std::endl;
    }
    else if (!strcasecmp(child->Value(),"gravity") && !myModel.empty())
      IFEM::cout <<"  ** The superelement <gravity> definition will be ignored"
                 <<"     unless you specify <elasticity> before the"
                 <<"     <geometry> block of the superelement."<< std::endl;

    else if (myProblem->parse(child))
      continue; // parsed <gravity>

    else if (!strcasecmp(child->Value(),"nodeload") && child->FirstChild())
    {
      IFEM::cout <<"  Parsing <nodeload>"<< std::endl;

      this->parseLoad(child);
    }
    else if (supEl && !supEl->parse(child)) // parse the superelement <geometry>
    {
      delete supEl;
      std::cerr <<" *** Failure."<< std::endl;
      return false;
    }

  if (supEl)
  {
    mySubSim[supId].sim = supEl;
    IFEM::cout <<"  FE model \""<< supId <<"\" loaded."<< std::endl;
  }

  return true;
}


bool SIMShellSup::preprocessB ()
{
  // Set up links to the underlying FE model for each superelement
  for (SuperElm& sup : mySups)
  {
    std::map<std::string,FEmodel>::const_iterator sit = mySubSim.find(sup.id);
    if (sit != mySubSim.end()) sup.sim = sit->second.sim;
  }

  // Preprocess the superelement FE models for the recovery process.
  // We need to maintain the original nodal ordering (only eliminating holes
  // in the sequence), otherwise the recovery will fail if the patch contains
  // beam elements which are placed in a separate patch.
  bool preserved = SIMbase::preserveNOrder;
  bool ok = SIMbase::preserveNOrder = true;
  for (std::pair<const std::string,FEmodel>& sub : mySubSim)
    if (sub.second.sim)
    {
      IFEM::cout <<"\nSetting up for recovery of superelement \""
                 << sub.first <<"\""<< std::endl;
      ok &= (sub.second.sim->preprocess({},fixDup) &&
             sub.second.sim->setMode(SIM::RECOVERY));
    }

  IFEM::cout <<"\nSuperelement recovery setup done."<< std::endl;
  SIMbase::preserveNOrder = preserved;
  return ok;
}


bool SIMShellSup::renumberNodes (const std::map<int,int>& nodeMap)
{
  return (this->SIMsupel::renumberNodes(nodeMap) &&
          this->renumberLoadedNodes(nodeMap));
}


bool SIMShellSup::assembleDiscreteTerms (const IntegrandBase* itg,
                                         const TimeDomain& time)
{
  return this->assembleLoads(*myEqSys,time.t);
}


bool SIMShellSup::recoverInternalDOFs (const ASMbase* pch, SuperElm& sup,
                                       const Vector& supSol) const
{
  const ASMsupel* supch = dynamic_cast<const ASMsupel*>(pch);
  if (supch && supch->hasRecovery())
   // Recover the internal displacement state for this superelement
   return supch->recoverInternals(supSol,sup.sol);

  return this->SIMsupel::recoverInternalDOFs(pch,sup,supSol);
}


ElementBlock* SIMShellSup::tesselatePatch (size_t pidx) const
{
  if (pidx >= mySups.size() || !mySups[pidx].sim)
    return this->SIMsupel::tesselatePatch(pidx); // use simplified tesselation

  FEmodel& sub = const_cast<SIMShellSup*>(this)->mySubSim[mySups[pidx].id];

  // Create an ElementBlock containing the shell patches in this superelement.
  // Normally, there should only be one such patch, but...
  if (sub.sim && !sub.blk)
    for (ASMbase* pch : sub.sim->getFEModel())
      if (pch && !pch->empty() && pch->getNoParamDim() > 1) // skip beam patches
      {
        bool ok = true;
        if (!sub.blk)
        {
          sub.blk = new ElementBlock(4);
          ok = pch->tesselate(*sub.blk,nullptr);
        }
        else
        {
          ElementBlock newblk(4);
          if ((ok = pch->tesselate(newblk,nullptr)))
            sub.blk->merge(newblk,false);
        }
        if (!ok)
        {
          delete sub.blk;
          return sub.blk = nullptr;
        }
      }

  if (!sub.blk)
    return nullptr;

  // Make a copy of sub.blk and apply the superelement transformation to it
  ElementBlock* supblk = new ElementBlock(*sub.blk);
  supblk->transform(mySups[pidx].MVP);
  return supblk;
}


/*!
  \brief Static helper to write out scalar fields to VTF-file.
*/

static bool writeFields (const Matrix& field, int geomID,
                         int& nBlock, std::vector<IntVec>& sID, VTF* vtf)
{
  for (size_t j = 1; j <= field.rows(); j++)
    if (!vtf->writeNres(field.getRow(j), ++nBlock, geomID))
      return false;
    else if (j <= sID.size())
      sID[j-1].push_back(nBlock);
    else
      sID.push_back({nBlock});

  return true;
}


int SIMShellSup::writeGlvS1 (const Vector& psol, int iStep, int& nBlock,
                             double, const char*, int idBlock, int, bool)
{
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return 0;
  else if (psol.empty())
    return 0;

  VTF* vtf = this->getVTF();
  if (!vtf) return -99;

  // Recover internal displacements for all superelements
  if (!this->recoverInternalDOFs(psol))
    return -9;

  IntVec vID;
  std::vector<IntVec> sID;
  sID.reserve(this->getNoFields());

  int geomID = this->getStartGeo();
  for (size_t pidx = 0; pidx < myModel.size(); pidx++)
  {
    Matrix field, subfield;

    if (mySups[pidx].sim)
    {
      for (ASMbase* pch : mySups[pidx].sim->getFEModel())
        if (pch && !pch->empty() && pch->getNoParamDim() > 1)
        {
          // Extract displacement vector for this sub-patch
          Vector pchvec;
          pch->extractNodalVec(mySups[pidx].sol, pchvec,
                               mySups[pidx].sim->getSAM()->getMADOF());
#if INT_DEBUG > 2
          std::cout <<"\nSolution vector for sub-patch "<< pch->idx+1 << pchvec;
#endif

          // Evaluate the internal displacement field on this sub-patch
          if (!pch->evalSolution(subfield,pchvec,(const int*)nullptr))
            return -5;

          if (field.empty())
            field = subfield;
          else
            field.augmentCols(subfield);
        }
    }
    else // Evaluate displacement field on supernodes only
      if (!myModel[pidx]->evalSolution(field, mySups[pidx].sol, opt.nViz))
        return -1;

    if (msgLevel > 1)
      IFEM::cout <<"Writing primary solution for patch "
                 << myModel[pidx]->idx+1 <<" ("<< field.rows()
                 <<","<< field.cols() <<")"<< std::endl;

    // Output as vector field
    if (!vtf->writeVres(field, ++nBlock, ++geomID, this->getNoSpaceDim()))
      return -2;
    else
      vID.push_back(nBlock);

    // Output as scalar fields
    if (!writeFields(field, geomID, nBlock, sID, vtf))
      return -3;
  }

  // Write result block identifications

  const char* fieldNames[6] = { "u_x", "u_y", "u_z", "r_x", "r_y", "r_z" };
  bool ok = vID.empty() || vtf->writeDblk(vID,"Displacement",idBlock,iStep);
  for (size_t i = 0; i < sID.size() && !sID[i].empty() && ok; i++)
    ok = vtf->writeSblk(sID[i], fieldNames[i], idBlock++, iStep);

  return ok ? idBlock : -4;
}
