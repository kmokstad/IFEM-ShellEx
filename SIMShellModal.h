// $Id$
//==============================================================================
//!
//! \file SIMShellModal.h
//!
//! \date Apr 28 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for modal linear shell FEM analysis.
//!
//==============================================================================

#ifndef _SIM_SHELL_MODAL_H
#define _SIM_SHELL_MODAL_H

#include "SIMAndesShell.h"
#include "SIMmodal.h"


/*!
  \brief Modal driver for FEM analysis of shell problems.
  \details The main feature of this class is that it overrides the
  SIMbase::assembleSystem() method, to deal with both assembly of the eigenvalue
  problem, and then also the modal equation system during the time integration.
*/

class SIMShellModal : public SIMAndesShell, public SIMmodal
{
public:
  //! \brief The constructor forwards to the parent class constructors.
  //! \param[in] modes Array of eigenmodes for the elasticity problem
  explicit SIMShellModal(std::vector<Mode>& modes)
    : SIMAndesShell(0,true), SIMmodal(modes) {}
  //! \brief Empty destructor.
  virtual ~SIMShellModal() {}

  using SIMAndesShell::assembleSystem;
  //! \brief Administers assembly of the linear equation system.
  //! \param[in] time Parameters for time-dependent simulations
  //! \param[in] mSol Previous modal solution vectors
  //!
  //! \details This method assembles the eigenvalue problem in the first call.
  //! Then it is used to build up the modal dynamic equation system during the
  //! time integration. The modal equation system is kept in a separate
  //! AlgEqSystem object (and an associated SAM object), and the pointers to
  //! those modal objects are swapped with the original ones depending on which
  //! system we are currently working on.
  virtual bool assembleSystem(const TimeDomain& time,
                              const Vectors& mSol, bool, bool);

  using SIMmodal::expandSolution;
  //! \brief Expands and returns the current dynamic solution.
  //! \param[in] mSol Current modal solution
  //! \param[in] swapBck If \e true, the equation systems are swapped
  virtual const Vectors& expandSolution(const Vectors& mSol, bool swapBck);

protected:
  using SIMAndesShell::parse;
  //! \brief Parses a data section from an XML element.
  //! \details Overrides the parent class method to do nothing when invoked
  //! during the second time parsing for the time integration setup only.
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details In addition to invoking the inherited method,
  //! this method sets the \a parsed flag of the parent class SIMmodal,
  //! such that the model parsing is skipped when the input file is parsed
  //! for the second time while doing the time integration setup.
  virtual bool preprocessB();
};

#endif
