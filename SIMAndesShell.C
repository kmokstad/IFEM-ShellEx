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
#include "IFEM.h"
#include "tinyxml2.h"


static bool withBeams = false; //!< If \e true, the model contains beam elements


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
    myProblem = new AndesShell(nss,modal,withBeams);

  return dynamic_cast<ElasticBase*>(myProblem);
}


bool SIMAndesShell::parse (const tinyxml2::XMLElement* elem)
{
  if (!this->SIMElasticity<SIM2D>::parse(elem))
    return false;

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(),"pressure") && child->FirstChild())
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

  if (!strcasecmp(elem->Value(),"postprocessing") && myProblem)
    for (child = elem->FirstChildElement(); child;
         child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"vonMises_only"))
        static_cast<AndesShell*>(myProblem)->vonMisesOnly();

  return true;
}


ASMbase* SIMAndesShell::readPatch (std::istream& isp, int pchInd,
                                   const CharVec&, const char* whiteSpace) const
{
  ASMbase* pch = NULL;
  bool nastran = nf.size() == 2 && nf[1] == 'n';
  if (nastran) // Nastran bulk data file
    pch = new ASMu2DNastran(nsd,nf.front());
  else if (!(pch = ASM2D::create(opt.discretization,nsd,nf)))
    return pch;

  if (!pch->read(isp) || this->getLocalPatchIndex(pchInd+1) < 1)
  {
    delete pch;
    return nullptr;
  }

  if (whiteSpace)
    IFEM::cout << whiteSpace <<"Reading patch "<< pchInd+1 << std::endl;

  ASMbase* bpch = nullptr;
  if (nastran) // Check if we also have beam elements in the model
    if ((bpch = static_cast<ASMu2DNastran*>(pch)->haveBeams()))
    {
      withBeams = true;
      bpch->idx = myModel.size();
      const_cast<SIMAndesShell*>(this)->myModel.push_back(bpch);
    }

  pch->idx = myModel.size();
  return pch;
}


void SIMAndesShell::getShellThicknesses (RealArray& elmThick) const
{
  // Include also the collapsed and non-shell elements,
  // because that is what the VTF-writer expects
  elmThick.resize(this->getNoElms(false,true),0.0);

  int iel = 0, missing = 0;
  for (const ASMbase* pch : myModel)
  {
    const ASMu2DNastran* shell = dynamic_cast<const ASMu2DNastran*>(pch);
    if (shell)
      for (size_t jel = 1; jel <= pch->getNoElms(true); iel++, jel++)
      {
        int ielNo = pch->getElmID(jel);
        if (ielNo > 0 && !shell->getThickness(ielNo,elmThick[iel]))
          ++missing;
      }
    else
      iel += pch->getNoElms(true);
  }

  if (missing > 1)
    IFEM::cout <<" *** A total of "<< missing <<" elements lack thickness.\n"
               <<"     Please check your model for inconsistency."<< std::endl;
}


bool SIMAndesShell::getElementGroup (int iset, std::string& name,
                                     RealArray& elGroup) const
{
  elGroup.clear();

  int iel = 0;
  for (const ASMbase* pch : myModel)
  {
    const ASMu2DLag* shell = dynamic_cast<const ASMu2DLag*>(pch);
    if (shell && shell->getElementSet(iset,name))
    {
      if (elGroup.empty())
        elGroup.resize(this->getNoElms(false,true),0.0);
      for (size_t je = 1; je <= pch->getNoElms(true); iel++, je++)
        elGroup[iel] = pch->getElmID(je) > 0 && shell->isInElementSet(iset,je);
    }
    else
      iel += pch->getNoElms(true);
  }

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
