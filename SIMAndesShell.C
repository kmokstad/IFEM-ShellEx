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
#include "Functions.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml2.h"


SIMAndesShell::SIMAndesShell (unsigned char n, bool m) : nsv(n), modal(m)
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
    myProblem = new AndesShell(nsv,modal);

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

  return true;
}


ASMbase* SIMAndesShell::readPatch (std::istream& isp, int pchInd,
                                   const CharVec&, const char* whiteSpace) const
{
  ASMbase* pch = NULL;
  if (nf.size() == 2 && nf[1] == 'n') // Nastran bulk data file
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

  pch->idx = myModel.size();
  return pch;
}


void SIMAndesShell::getShellThicknesses (RealArray& elmThick) const
{
  int iel = 0, missing = 0;
  for (const ASMbase* pch : myModel)
  {
    const ASMu2DNastran* shell = dynamic_cast<const ASMu2DNastran*>(pch);
    if (shell)
      for (double& t : elmThick)
        if (!shell->getThickness(pch->getElmID(++iel),t))
          ++missing;
  }
  if (missing > 1)
    IFEM::cout <<" *** A total of "<< missing <<" elements lack thickness.\n"
               <<"     Please check your model for inconsistency."<< std::endl;
}


bool SIMAndesShell::renumberNodes (const std::map<int,int>& nodeMap)
{
  bool ok = this->SIMElasticity<SIM2D>::renumberNodes(nodeMap);

  for (PointLoad& load : myLoads)
    if (load.inod > 0)
      ok &= utl::renumber(load.inod,nodeMap,true);

  return ok;
}


bool SIMAndesShell::assembleDiscreteTerms (const IntegrandBase* itg,
                                           const TimeDomain& time)
{
  if (itg != myProblem || !myEqSys)
    return true;

  bool ok = true;
  SystemVector* R = myEqSys->getVector(0);
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
