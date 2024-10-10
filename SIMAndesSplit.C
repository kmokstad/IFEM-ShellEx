// $Id$
//==============================================================================
//!
//! \file SIMAndesSplit.C
//!
//! \date Oct 10 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for FE analysis using the ANDES shell elements.
//!
//==============================================================================

#include "SIMAndesSplit.h"
#include "AndesShell.h"
#include "ASMbase.h"
#include "AlgEqSystem.h"
#include "SparseMatrix.h"
#include "SAM.h"
#include "IFEM.h"


SIMAndesSplit::~SIMAndesSplit ()
{
  for (Region& material : myRegions)
    delete material.myEqSys;
}


bool SIMAndesSplit::preprocessB ()
{
  AndesShell* shellp = dynamic_cast<AndesShell*>(myProblem);
  if (!shellp)
  {
    std::cerr<<" *** SIMAndesSplit::preprocessB: No shell problem yet"
             << std::endl;
    return false;
  }

  for (const ASMbase* pch : myModel)
    for (size_t iel = 1; iel <= pch->getNoElms(true); iel++)
      if (shellp->isInLossArea(pch->getElementCenter(iel)))
        myRegions[1].myElements.push_back(pch->getElmID(iel));
      else
        myRegions[0].myElements.push_back(pch->getElmID(iel));

  IFEM::cout <<"\nSplitting the domain into two material regions: "
             << myRegions.front().myElements.size() <<" "
             << myRegions.back().myElements.size();
#ifdef INT_DEBUG
  int count = 0;
  const char* newline = "\n                           ";
  IFEM::cout <<"\nElements in second region:";
  for (int e :  myRegions.back().myElements)
    std::cout << (++count == 1 || (count-1)%10 ? " " : newline) << e;
#endif
  IFEM::cout << std::endl;

  return true;
}


bool SIMAndesSplit::assembleSystem (const TimeDomain& time,
                                    const Vectors& prevSol,
                                    bool newLHS, bool poorConvg)
{
  if (myRegions[1].myElements.empty()) // We have only one material region
    return this->SIMAndesShell::assembleSystem(time,prevSol,newLHS,poorConvg);

  // Assemble the two material regions separately
  AlgEqSystem* fullEqSys = myEqSys;
  for (Region& material : myRegions)
  {
    // Swap equation system to be assembled
    if (!material.myEqSys)
      material.myEqSys = new AlgEqSystem(*fullEqSys);
    myEqSys = material.myEqSys;

    // Activate elements in this material region
    for (ASMbase* pch : myModel) pch->setActiveElements(&material.myElements);

    // Assemble this material region
    if (!this->SIMAndesShell::assembleSystem(time,prevSol,newLHS,poorConvg))
      return false;

    // Output separated matrices for the ROM analysis
    this->dumpEqSys(msgLevel > 1);
  }

  // Combine the two regions into one
  myEqSys = fullEqSys;
  myEqSys->copy(*myRegions[0].myEqSys).add(*myRegions[1].myEqSys);

  // Restore activation of all elements
  for (ASMbase* pch : myModel) pch->setActiveElements(nullptr);

  return true;
}


bool SIMAndesSplit::assembleDiscreteTerms (const IntegrandBase* itg,
                                           const TimeDomain& time)
{
  // This assumes we have point load(s) only on the first material region
  if (myEqSys == myRegions.front().myEqSys || !myRegions.front().myEqSys)
    return this->SIMAndesShell::assembleDiscreteTerms(itg,time);

  return true;
}


bool SIMAndesSplit::solveSystem (Vector& solution, int printSol, double* rCond,
                                 const char* compName, size_t idxRHS)
{
  std::vector<int> singularDofs;
  SparseMatrix* Kptr = dynamic_cast<SparseMatrix*>(myEqSys->getMatrix());
  if (Kptr && !Kptr->isFactored())
  {
    // Check if the stiffness matrix is singular (due to absent material)
    // And fill in with non-zero diagonal elements where needed
    SparseMatrix& Kmat = *Kptr;
    double Kstiff = Kmat.Linfnorm();
    for (size_t ieq = 1; ieq <= Kmat.dim(); ieq++)
      if (fabs(Kmat(ieq,ieq)) < 1.0e-12)
      {
#if INT_DEBUG > 1
        std::cout <<"Resetting pivot element for equation "<< ieq
                  <<" to "<< Kstiff << std::endl;
#endif
        singularDofs.push_back(ieq);
        Kmat(ieq,ieq) = Kstiff;
      }
  }

  bool ok = this->solveEqSystem(solution,idxRHS,rCond,printSol,true,compName);

  if (!singularDofs.empty())
  {
    // Reset solution at the singular DOFs to zero
    const int* meqn = this->getSAM()->getMEQN();
    const int nDofs = this->getSAM()->getNoDOFs();
    for (int ieq : singularDofs)
    {
      int iDof = std::find(meqn,meqn+nDofs,ieq) - meqn;
      if (iDof >= 0 && iDof < nDofs)
        solution[iDof] = 0.0;
    }
  }

  return ok;
}
