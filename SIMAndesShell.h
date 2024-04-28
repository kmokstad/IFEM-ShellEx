// $Id$
//==============================================================================
//!
//! \file SIMAndesShell.h
//!
//! \date Feb 20 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for FE analysis using the ANDES shell elements.
//!
//==============================================================================

#ifndef _SIM_ANDES_SHELL_H
#define _SIM_ANDES_SHELL_H

#include "SIMElasticity.h"
#include "SIM2D.h"


/*!
  \brief Driver class for FE analysis using the ANDES shell elements.
*/

class SIMAndesShell : public SIMElasticity<SIM2D>
{
public:
  //! \brief Default constructor.
  //! \param[in] n Number of consequtive solution vectors in core
  //! \param[in] m If \e true, a modal linear dynamics simulation is performed
  explicit SIMAndesShell(unsigned char n = 1, bool m = false);
  //! \brief Empty destructor.
  virtual ~SIMAndesShell() {}

protected:
  using SIMElasticity<SIM2D>::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Reads a patch from given input stream.
  //! \param[in] isp The input stream to read from
  //! \param[in] pchInd 0-based index of the patch to read
  //! \param[in] whiteSpace For message formatting
  virtual ASMbase* readPatch(std::istream& isp, int pchInd, const CharVec&,
                             const char* whiteSpace) const;

  //! \brief Returns the actual integrand.
  virtual ElasticBase* getIntegrand();

  //! \brief Dummy override, does nothing.
  virtual bool initBodyLoad(size_t) { return true; }

private:
  unsigned char nsv; //!< Number of consequtive solution vectors in core

  bool modal; //!< Modal dynamics simulation flag
};

#endif
