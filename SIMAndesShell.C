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
#include "IFEM.h"


SIMAndesShell::SIMAndesShell (unsigned char n) : nsv(n)
{
  nsd = 3;
  nf.front() = 6;
}


ElasticBase* SIMAndesShell::getIntegrand ()
{
  if (!myProblem)
    myProblem = new AndesShell(nsv);

  return dynamic_cast<ElasticBase*>(myProblem);
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
