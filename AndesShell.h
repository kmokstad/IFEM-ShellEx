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
#include <set>

class ASMu2DNastran;
class RealFunc;


/*!
  \brief Class representing the integrand of an ANDES shell.
*/

class AndesShell : public ElasticBase
{
public:
  //! \brief Default constructor.
  //! \param[in] n Number of consequtive solution vectors to reside in core
  explicit AndesShell(unsigned short int n = 1);
  //! \brief The destructor writes out the ignored bad elements.
  virtual ~AndesShell();

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Initialization of integrand with patch-specific data.
  virtual void initForPatch(const ASMbase* pch);

  //! \brief Defines the pressure field.
  void setPressure(RealFunc* pf = nullptr);

  using ElasticBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t, bool) const;

  using ElasticBase::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] fe Nodal and integration point data for current element
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC,
                           const FiniteElement& fe, const Vec3&, size_t,
                           LocalIntegral& elmInt);

  using ElasticBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

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

  const ASMu2DNastran* currentPatch; //!< Pointer to underlying FE model

  std::vector<RealFunc*> presFld; //!< Pointers to pressure field functions

  std::set<int> degenerated;  //!< List of detected degenerated elements
  std::set<int> straightline; //!< List of detected straight line elements
};

#endif
