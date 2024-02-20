// $Id$
//==============================================================================
//!
//! \file AndesShell.h
//!
//! \brief Class representing a linear elastic shell with rotational DOFs.
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \date Feb 15 2024
//!
//==============================================================================

#ifndef _ANDES_SHELL_H
#define _ANDES_SHELL_H

#include "ElasticBase.h"


/*!
  \brief Class representing the integrand of an ANDES shell.
*/

class AndesShell : public ElasticBase
{
public:
  //! \brief Default constructor.
  //! \param[in] n Number of consequtive solution vectors to reside in core
  explicit AndesShell(unsigned short int n = 1);
  //! \brief Empty destructor.
  virtual ~AndesShell() {}

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  using ElasticBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t, bool) const;

  using ElasticBase::finalizeElement;
  //! \brief Finalizes the element matrices after the numerical integration.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Nodal and integration point data for current element
  virtual bool finalizeElement(LocalIntegral& elmInt, const FiniteElement& fe,
                               const TimeDomain&, size_t);

private:
  double Thick; //!< Current shell thickness
  double Emod;  //!< Current Young's modules
  double Rny;   //!< Current Poisson's ratio
  double Rho;   //!< Current mass density
};

#endif
