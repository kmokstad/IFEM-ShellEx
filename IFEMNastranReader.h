// $Id$
//==============================================================================
//!
//! \file IFEMNastranReader.h
//!
//! \date Feb 27 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief A Nastran Bulk Data File (BDF) parser for IFEM.
//!
//==============================================================================

#ifndef _IFEM_NASTRAN_READER_H
#define _IFEM_NASTRAN_READER_H

#include "FFlNastranReader.H"


/*!
  \brief Class for reading Nastran bulk data into a FFlLinkHandler object.
*/

class IFEMNastranReader : public FFlNastranReader
{
  int nPreBulk; //!< Number of lines before BEGIN BULK

public:
  //! \brief The constructor forwards to the parent class constructor.
  IFEMNastranReader(FFlLinkHandler& fePart, int lCount)
    : FFlNastranReader(&fePart,lCount), nPreBulk(lCount) {}

  //! \brief Reads the FE model from the Nastran file stream.
  bool readFE(std::istream& is, std::istream& iset, bool skipBeams = false);

  //! \brief Reads through a Nastran file to extract pre-bulk SET definitions.
  static int parseSets(std::istream& is, std::ostream& os,
                       const char* fileName = nullptr);
};

#endif
