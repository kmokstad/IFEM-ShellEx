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

class DataExporter;


/*!
  \brief Driver class for FE analysis using the ANDES shell elements.
*/

class SIMAndesShell : public SIMElasticity<SIM2D>
{
public:
  //! \brief Default constructor.
  //! \param[in] n Number of consecutive solutions in core (0 = linear analysis)
  //! \param[in] m If \e true, a modal linear dynamics simulation is performed
  explicit SIMAndesShell(unsigned short int n = 0, bool m = false);
  //! \brief The destructor deletes the nodal point load functions.
  virtual ~SIMAndesShell();

  //! \brief Creates a HDF5 data exporter for the simulator.
  //! \param[in] psol Primary solution vector
  //! \param[in] dumpNodeMap If \e true, write node mapping to the HDF5 as well
  DataExporter* getHDF5writer(const Vector& psol, double dumpNodeMap) const;

  //! \brief Retrieves the shell thickness of all elements in the model.
  void getShellThicknesses(RealArray& elmThick) const;
  //! \brief Retrieves the specified element group \a iset as a scalar field.
  bool getElementGroup(int iset, std::string& name, RealArray& elGroup) const;

  //! \brief Writes additional geometries illustrating sensor locations.
  //! \param[in] locfiles Files with sensor locations (node or element IDs)
  //! \param[in] nodal If \e true, nodal sensors. Othwerwise element sensors.
  //! \param nBlock Running geometry block counter
  bool writeGlvLoc(std::vector<std::string>& locfiles,
                   bool nodal, int& nBlock) const;

protected:
  using SIMElasticity<SIM2D>::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Reads a patch from given input stream.
  virtual ASMbase* readPatch(std::istream& isp,
                             int, const CharVec&, const char*) const;

  //! \brief Returns the actual integrand.
  virtual ElasticBase* getIntegrand();

  //! \brief Dummy override, does nothing.
  virtual bool initBodyLoad(size_t) { return true; }

  //! \brief Renumbers the global node numbers of the nodal point loads.
  virtual bool renumberNodes(const std::map<int,int>& nodeMap);

  //! \brief Assembles the nodal point loads, if any.
  virtual bool assembleDiscreteTerms(const IntegrandBase* itg,
                                     const TimeDomain& time);

public:
  static char useBeams; //!< If non-zero, include beam elements, if any
  static bool readSets; //!< If \e true, also read Nastran SET definitions

private:
  //! \brief Struct defining a nodal point load.
  struct PointLoad
  {
    int         inod; //!< Node or patch index
    int         ldof; //!< Local DOF number
    ScalarFunc* p;    //!< Load magnitude
    //! \brief Default constructor.
    PointLoad(int n = 0) : inod(n), ldof(0), p(nullptr) {}
  };

  std::vector<PointLoad> myLoads; //!< Nodal point loads

  unsigned short int nss; //!< Number of consequtive solution states in core

  bool modal; //!< Modal dynamics simulation flag
};

#endif
