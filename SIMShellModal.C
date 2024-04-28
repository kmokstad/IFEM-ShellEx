// $Id$
//==============================================================================
//!
//! \file SIMShellModal.C
//!
//! \date Apr 28 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for modal linear shell FEM analysis.
//!
//==============================================================================

#include "SIMShellModal.h"
#include "IntegrandBase.h"


bool SIMShellModal::assembleSystem (const TimeDomain& time,
                                    const Vectors& mSol, bool, bool)
{
  // Assemble the eigenvalue system
  if (myProblem->getMode() == SIM::VIBRATION)
    return this->SIM2D::assembleSystem(time,Vectors());

  if (time.it > 0)
    // Swap back to the full equation system for assembly of load vector
    this->swapSystem(myEqSys,mySam);

  else
  {
    // Assemble the load vector of this time step.
    // We need to do this in the first iteration only, as for linear systems
    // the load vector is not supposed to change during the iterations.
    if (!this->SIM2D::assembleSystem(time,sol,false))
      return false;

    // Extract the load vector in DOF-order
    if (!this->extractLoadVec(Rhs))
      return false;
  }

  // Assemble the modal equation system
  if (!this->assembleModalSystem(time,mSol,
                                 myProblem->getIntegrationPrm(2),
                                 myProblem->getIntegrationPrm(3)))
    return false;

  // Swap the equation systems such that the dynamic simulation driver
  // operates on the modal system
  return this->swapSystem(myEqSys,mySam);
}


const Vectors& SIMShellModal::expandSolution (const Vectors& mSol, bool swapBck)
{
  // Swap back to the full equation system data for postprocessing
  // and assembly of load vector for the next time step
  if (swapBck)
    this->swapSystem(myEqSys,mySam);

  return this->expandSolution(mSol);
}


bool SIMShellModal::parse (const tinyxml2::XMLElement* elem)
{
  return this->parseParams(elem) || this->SIMAndesShell::parse(elem);
}


bool SIMShellModal::preprocessB ()
{
  parsed = true;
  this->setIntegrationPrm(0,alpha1);
  this->setIntegrationPrm(1,alpha2);
  return this->SIMAndesShell::preprocessB();
}
