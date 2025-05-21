// $Id$
//==============================================================================
//!
//! \file SIMAndesShell.C
//!
//! \date Feb 20 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for FE analysis using the ANDES shell elements.
//!
//==============================================================================

#include "SIMAndesShell.h"
#include "ASMu2DNastran.h"
#include "AndesShell.h"
#include "AlgEqSystem.h"
#include "SAM.h"
#include "DataExporter.h"
#include "HDF5Writer.h"
#include "Functions.h"
#include "Utilities.h"
#include "ElementBlock.h"
#include "VTF.h"
#include "IFEM.h"
#include "tinyxml2.h"
#include <fstream>


char SIMAndesShell::useBeams = 1;
bool SIMAndesShell::readSets = true;


SIMAndesShell::SIMAndesShell (unsigned short int n, bool m) : nss(n), modal(m)
{
  nsd = 3;
  nf.front() = 6;
}


SIMAndesShell::~SIMAndesShell ()
{
  for (PointLoad& load : myLoads)
    delete load.p;
}


ElasticBase* SIMAndesShell::getIntegrand ()
{
  if (!myProblem)
    myProblem = new AndesShell(nss,modal,useBeams);

  return dynamic_cast<ElasticBase*>(myProblem);
}


bool SIMAndesShell::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),"geometry"))
  {
    const tinyxml2::XMLElement* child = elem->FirstChildElement("fixRBE3");
    if (child && child->FirstChild())
      utl::parseIntegers(ASMu2DNastran::fixRBE3,child->FirstChild()->Value());
  }

  if (!this->SIMElasticity<SIM2D>::parse(elem))
    return false;

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(elem->Value(),"postprocessing"))
    {
      if (myProblem && !strcasecmp(child->Value(),"vonMises_only"))
        static_cast<AndesShell*>(myProblem)->vonMisesOnly();
    }
    else if (strcasecmp(elem->Value(),"elasticity"))
      continue; // The remaining should be within the elasticity context

    else if (!strcasecmp(child->Value(),"pressure") && child->FirstChild())
    {
      IFEM::cout <<"  Parsing <pressure>"<< std::endl;

      std::string set, type;
      utl::getAttribute(child,"set",set);
      int code = this->getUniquePropertyCode(set,1);
      if (code == 0) utl::getAttribute(child,"code",code);
      if (code > 0)
      {
        utl::getAttribute(child,"type",type,true);
        IFEM::cout <<"\tPressure code "<< code;
        if (!type.empty()) IFEM::cout <<" ("<< type <<")";
        myScalars[code] = utl::parseRealFunc(child->FirstChild()->Value(),type);
        this->setPropertyType(code,Property::BODYLOAD);
        IFEM::cout << std::endl;
        AndesShell* shellp = static_cast<AndesShell*>(this->getIntegrand());
        if (set.empty()) // Applies to all elements in the model
          shellp->setPressure(myScalars[code],code);
        else for (const ASMbase* pch : myModel)
          if (shellp->setPressure(myScalars[code],code,set,pch))
            break; // Note: This assumes a set has elements from one patch only
      }
    }
    else if (!strcasecmp(child->Value(),"nodeload") && child->FirstChild())
    {
      IFEM::cout <<"  Parsing <nodeload>"<< std::endl;

      PointLoad load;
      utl::getAttribute(child,"node",load.inod);
      utl::getAttribute(child,"dof",load.ldof);

      if (load.inod > 0 && load.ldof > 0 && load.ldof <= 6)
      {
        std::string type("constant");
        utl::getAttribute(child,"type",type);

        IFEM::cout <<"\tNode "<< load.inod <<" dof "<< load.ldof <<" Load: ";
        if (type == "constant")
        {
          load.p = new ConstantFunc(atof(child->FirstChild()->Value()));
          IFEM::cout << (*load.p)(0.0) << std::endl;
        }
        else
          load.p = utl::parseTimeFunc(child->FirstChild()->Value(),type);

        myLoads.push_back(load);
      }
    }
    else if (!strcasecmp(child->Value(),"material"))
    {
      IFEM::cout <<"  Parsing <material>"<< std::endl;

      this->getIntegrand()->parseMatProp(child,true);
    }

  return true;
}


ASMbase* SIMAndesShell::readPatch (std::istream& isp, int, const CharVec&,
                                   const char*) const
{
  ASMbase* pch = NULL;
  ASMu2DNastran* shell = NULL;
  if (nf.size() == 2 && nf[1] == 'n') // Nastran bulk data file
    pch = shell = new ASMu2DNastran(nsd,nf.front(),readSets,useBeams);
  else if (!(pch = ASM2D::create(opt.discretization,nsd,nf)))
    return pch;

  if (!pch->read(isp))
  {
    delete pch;
    return nullptr;
  }

  if (shell)
  {
    // Check if we also have beam elements in the model.
    // They will be kept in a separate patch of 1D elements.
    ASMbase* bpch = shell->haveBeams();
    if (bpch)
    {
      bpch->idx = myModel.size();
      const_cast<SIMAndesShell*>(this)->myModel.push_back(bpch);
    }
    else
      useBeams = 0;
  }

  pch->idx = myModel.size();
  if (shell && useBeams && pch->getNoElms() > 0)
    IFEM::cout <<"\tCreated shell patch "<< pch->idx+1
               <<" with "<< pch->getNoElms() <<" elements"<< std::endl;
  else
    IFEM::cout <<"\tReading patch "<< pch->idx+1 << std::endl;

  if (shell)
  {
    // Add topology items for the predefined node- and element sets
    std::string name;
    SIMinput* sim = const_cast<SIMAndesShell*>(this);
    for (int iset = 1; shell->getNodeSet(iset,name); iset++)
      sim->topology(name).emplace(pch->idx+1,iset,4);
    for (int iset = 1; shell->getElementSet(iset,name); iset++)
      sim->topology(name).emplace(pch->idx+1,iset,5);
  }

  return pch;
}


void SIMAndesShell::getShellThicknesses (RealArray& elmThick) const
{
  // Include also the collapsed and non-shell elements,
  // because that is what the VTF-writer expects
  elmThick.clear();
  elmThick.reserve(this->getNoElms(false,true));

  int missing = 0;
  const ASMu2DNastran* shell;
  for (const ASMbase* pch : myModel)
    if ((shell = dynamic_cast<const ASMu2DNastran*>(pch)))
    {
      for (size_t iel = 1; iel <= pch->getNoElms(true); iel++)
        if (pch->getElementNodes(iel).size() > 1) // skip 1-noded mass elements
        {
          elmThick.push_back(0.0);
          int ielNo = pch->getElmID(iel);
          if (ielNo > 0 && !shell->getThickness(ielNo,elmThick.back()))
            ++missing;
        }
    }
    else
      elmThick.insert(elmThick.end(),pch->getNoElms(true),0);

  if (missing > 1)
    IFEM::cout <<" *** A total of "<< missing <<" elements lack thickness.\n"
               <<"     Please check your model for inconsistency."<< std::endl;
}


bool SIMAndesShell::getElementGroup (int iset, std::string& name,
                                     RealArray& elGroup) const
{
  elGroup.clear();
  elGroup.reserve(this->getNoElms(false,true));

  const ASMu2DLag* shell;
  for (const ASMbase* pch : myModel)
    if ((shell = dynamic_cast<const ASMu2DLag*>(pch)) &&
        shell->getElementSet(iset,name))
    {
      for (size_t iel = 1; iel <= pch->getNoElms(true); iel++)
        if (pch->getElementNodes(iel).size() > 1) // skip 1-noded mass elements
          elGroup.push_back(pch->getElmID(iel) > 0 &&
                            shell->isInElementSet(iset,iel));
    }
    else
      elGroup.insert(elGroup.end(),pch->getNoElms(true),0);

  return !elGroup.empty();
}


bool SIMAndesShell::renumberNodes (const std::map<int,int>& nodeMap)
{
  bool ok = this->SIMElasticity<SIM2D>::renumberNodes(nodeMap);

  for (PointLoad& load : myLoads)
    if (load.inod > 0)
      ok &= utl::renumber(load.inod,nodeMap,true);

  return ok;
}


/*!
  This method is invoked after the assembly loop over the elements is finished.
  It will therefore check if there were any failed elements during the assembly,
  and return \e false in that case such that the simulation will be aborted.
*/

bool SIMAndesShell::assembleDiscreteTerms (const IntegrandBase* itg,
                                           const TimeDomain& time)
{
  bool ok = static_cast<AndesShell*>(myProblem)->allElementsOK();

  if (itg != myProblem || !myEqSys || !myEqSys->getNoRHS())
    return ok;

  SystemVector* R = myEqSys->getVector(myEqSys->getNoRHS()-1);
  if (R) // Assemble external nodal point loads
    for (const PointLoad& load : myLoads)
    {
      double P = (*load.p)(time.t);
      int ldof = load.ldof;
      myEqSys->addScalar(P,ldof-1);
      ok &= mySam->assembleSystem(*R,P,{load.inod,ldof});
    }

  return ok;
}


DataExporter* SIMAndesShell::getHDF5writer (const Vector& psol,
                                            double dumpNodeMap) const
{
  IFEM::cout <<"\nWriting HDF5 file "<< opt.hdf5 <<".hdf5"<< std::endl;

  DataExporter* writer = new DataExporter(true,opt.saveInc);
  writer->registerWriter(new HDF5Writer(opt.hdf5,adm));

  int result = DataExporter::PRIMARY | DataExporter::DISPLACEMENT;
  if (!opt.pSolOnly) result |= DataExporter::SECONDARY;
  if (dumpNodeMap)   result |= DataExporter::L2G_NODE;
  writer->registerField("u","solution",DataExporter::SIM,result);
  writer->setFieldValue("u",this,&psol);

  return writer;
}


bool SIMAndesShell::writeGlvLoc (std::vector<std::string>& locfiles,
                                 bool nodalR, int& nBlock) const
{
  ElementBlock* sensor;
  const ASMu2DNastran* shell;
  for (const std::string& fName : locfiles)
  {
    ElementBlock* sensorBlock = new ElementBlock(8);
    std::ifstream locs(fName);
    while (locs.good())
    {
      Vec3 XYZloc;
      int idx = 0;
      locs >> idx;
      for (const ASMbase* pch : myModel)
        if ((shell = dynamic_cast<const ASMu2DNastran*>(pch)) &&
            (sensor = shell->sensorGeometry(idx,nodalR)))
        {
          sensorBlock->merge(*sensor,false);
          break;
        }
    }

    if (sensorBlock->getNoElms() < 1)
      delete sensorBlock;
    else if (!this->getVTF()->writeGrid(sensorBlock,fName.c_str(),++nBlock))
      return false;
  }

  return true;
}
