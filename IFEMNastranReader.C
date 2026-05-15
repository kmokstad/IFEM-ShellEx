// $Id$
//==============================================================================
//!
//! \file IFEMNastranReader.C
//!
//! \date Feb 27 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief A Nastran Bulk Data File (BDF) parser for IFEM.
//!
//==============================================================================

#include <cstring>

#include "IFEMNastranReader.h"
#include "FFlLinkHandler.H"
#include "FFlElementBase.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"
#ifdef HAS_PROFILER
#include "Profiler.h"
#else
//! \cond DO_NOT_DOCUMENT
#define PROFILE(s)
//! \endcond
#endif

namespace ASM {
  char skipMass = 0; //!< 1: ignore one-noded mass elements, 2: no VTF output
  bool skipRBE2 = false; //!< If \e true, ignore all RBE2 elements
  bool readSets = true; //!< If \e true, read pre-bulk Nastran SET definitions
  std::vector<int> ignoredElms; //!< List of elements to ignore
}


int IFEMNastranReader::parseSets (std::istream& is, std::ostream& os,
                                  const char* fileName)
{
  // Fast-forward until "BEGIN BULK"
  int lCount = 0;
  char cline[256];
  while (is.getline(cline,255))
    if (!strncmp(cline,"BEGIN BULK",10))
      break;
    else
    {
      ++lCount;
      // Copy element SET definitions to a second stream,
      // since they have to be parsed after the FE model is loaded
      if (ASM::readSets && (!strncmp(cline,"SET ",4) || os.tellp() > 0))
        os << cline << '\n';
    }

  if (!is)
  {
    if (fileName)
      std::cerr <<"\n *** Parsing FE data file "<< fileName
                <<" failed."<< std::endl;
    return -1; // No bulk data file, or second pass when parsing for more
  }
  else if (fileName)
    ListUI <<"Parsing FE data file "<< fileName <<"\n";

  if (lCount > 0)
    ListUI <<"\tNastran bulk data starting at line "<< lCount+1 <<"\n";

  return lCount;
}


bool IFEMNastranReader::readFE (std::istream& is, std::istream& iset,
                                bool skipBeams)
{
  PROFILE("Nastran file parser");

  if (!this->resolve(this->read(is)))
    myLink->deleteGeometry(); // Parsing failure, delete all FE data
  else if (nWarnings+nNotes > 0)
    ListUI <<"\n  ** Parsing FE data succeeded."
           <<"\n     However, "<< nWarnings
           <<" warning(s) and "<< nNotes <<" note(s) were reported.\n"
           <<"     Review the messages and check the FE data file.\n\n";
  if (!myLink->hasGeometry())
    return false;

  // Remove all solid elements (not yet supported),
  // and (optionally) all the mass and beam elements, and ...
  ElementsVec toBeErased;
  ElementsCIter it;
  for (it = myLink->elementsBegin(); it != myLink->elementsEnd(); ++it)
    if ((*it)->getCathegory() == FFlTypeInfoSpec::SOLID_ELM ||
        std::find(ASM::ignoredElms.begin(),ASM::ignoredElms.end(),
                  (*it)->getID()) != ASM::ignoredElms.end() ||
        (ASM::skipRBE2 && (*it)->getTypeName() == "RGD") ||
        (ASM::skipMass == 1 && (*it)->getTypeName() == "CMASS") ||
        (skipBeams && (*it)->getCathegory() == FFlTypeInfoSpec::BEAM_ELM))
        toBeErased.push_back(*it);
  if (!toBeErased.empty())
  {
    ListUI <<"  ** Erasing "<< toBeErased.size()
           <<" elements from the model.\n";
    myLink->removeElements(toBeErased);
  }

  // Now parse the element set definitions, if any
  lastComment = { 0, "" };
  return iset ? this->processAllSets(iset,nPreBulk) : true;
}
