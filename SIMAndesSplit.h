// $Id$
//==============================================================================
//!
//! \file SIMAndesSplit.h
//!
//! \date Oct 10 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for FE analysis using the ANDES shell elements.
//!
//==============================================================================

#ifndef _SIM_ANDES_SPLIT_H
#define _SIM_ANDES_SPLIT_H

#include "SIMAndesShell.h"
#include <array>


/*!
  \brief Driver class for FE analysis using the ANDES shell elements.
*/

class SIMAndesSplit : public SIMAndesShell
{
public:
  //! \brief Default constructor.
  explicit SIMAndesSplit(unsigned short int n = 0) : SIMAndesShell(n) {}
  //! \brief The destructor deallocates the per region equation systems.
  virtual ~SIMAndesSplit();

  //! \brief Administers assembly of the linear equation system.
  virtual bool assembleSystem(const TimeDomain& time, const Vectors& pSol,
                              bool newLHS, bool poorConvg);

  //! \brief Solves the assembled linear system of equations for a given load.
  virtual bool solveSystem(Vector& solution, int printSol, double* rCond,
                           const char* compName, size_t idxRHS);

protected:
  //! \brief Assembles the nodal point loads, if any.
  virtual bool assembleDiscreteTerms(const IntegrandBase* itg,
                                     const TimeDomain& time);

  //! \brief Preprocessing performed after the system assembly initialization.
  virtual bool preprocessB();

private:
  //! \brief Struct defining a material region.
  struct Region
  {
    std::vector<int> myElements;        //!< The elements defining this region
    AlgEqSystem*     myEqSys = nullptr; //!< System matrices for this region
  };

  std::array<Region,2> myRegions; //!< Separate material regions
};

#endif
