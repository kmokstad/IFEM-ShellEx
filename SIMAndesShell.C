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
#include "Functions.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml2.h"


SIMAndesShell::SIMAndesShell (unsigned char n, bool m) : nsv(n), modal(m)
{
  nsd = 3;
  nf.front() = 6;
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
