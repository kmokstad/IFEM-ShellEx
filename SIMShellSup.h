// $Id$
//==============================================================================
//!
//! \file SIMShellSup.h
//!
//! \date May 24 2025
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for linear elastic superelement FEM analysis.
//!
//==============================================================================

#ifndef _SIM_SHELL_SUP_H
#define _SIM_SHELL_SUP_H

#include "SIMsupel.h"


/*!
  \brief Solution driver for linear elastic superelement FEM analysis.
*/

class SIMShellSup : public SIMsupel
{
public:
  //! \brief The constructor creates a dummy integrand with gravity only.
  //! \param[in] hd Sub-simulator heading
  //! \param[in] fd If \e true, merge duplicated FE nodes on patch interfaces
  SIMShellSup(const char* hd, bool fd);
  //! \brief The destructor deletes the FE substructure data.
  virtual ~SIMShellSup();

protected:
  using SIMsupel::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is overridden to preprocess the FE substructures.
  virtual bool preprocessB();

  using SIMsupel::recoverInternalDOFs;
  //! \brief Recovers the internal DOFs from static condensation.
  //! \param[in] pch The patch associated with the superelement \a sup
  //! \param sup The superelement to recover internal DOFs for
  //! \param[in] supSol Supernode solution values
  virtual bool recoverInternalDOFs(const ASMbase* pch, SuperElm& sup,
                                   const Vector& supSol) const;

  //! \brief Tesselates the superelement associated with specified patch.
  virtual ElementBlock* tesselatePatch(size_t pidx) const;

  //! \brief Returns a pointer to the problem-specific data object.
  IntegrandBase* getMyProblem() const;

  //! \brief Writes primary solution for a given time step to the VTF-file.
  //! \param[in] psol Primary solution vector
  //! \param[in] iStep Time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] idBlock Starting value of result block numbering
  virtual int writeGlvS1(const Vector& psol, int iStep, int& nBlock,
                         double, const char*, int idBlock, int, bool);

private:
  //! \brief Struct representing a FE substructure.
  struct FEmodel
  {
    SIMgeneric*   sim = nullptr; //!< Underlying FE model of the substructure
    ElementBlock* blk = nullptr; //!< Tesselated FE model for visualization
    std::vector<int> superNodes; //!< List of external supernode numbers
  };

  std::map<std::string,FEmodel> mySubSim; //!< FE substructure container

  bool fixDup; //!< If \e true, merge duplicated FE nodes on patch interfaces
};

#endif
