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
  //! \param[in] n Number of consequtive solution vectors in core
  //! \param[in] m If \e true, a modal linear dynamics simulation is performed
  explicit SIMAndesShell(unsigned char n = 1, bool m = false);
  //! \brief The destructor deletes the nodal point load functions.
  virtual ~SIMAndesShell();

  //! \brief Creates a HDF5 data exporter for the simulator.
  //! \param[in] psol Primary solution vector
  //! \param[in] dumpNodeMap If \e true, write node mapping to the HDF5 as well
  DataExporter* getHDF5writer(const Vector& psol, double dumpNodeMap) const;

  //! \brief Retrieves the shell thickness of all elements in the model.
  void getShellThicknesses(RealArray& elmThick) const;
  //! \brief Retrieves the specified element group \a idx as a scalar field.
  bool getElementGroup(int idx, RealArray& elGroup) const;

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

  //! \brief Renumbers all global nodes number if the model.
  //! \param[in] nodeMap Mapping from old to new node number
  virtual bool renumberNodes(const std::map<int,int>& nodeMap);

  //! \brief Assembles the nodal point loads, if any.
  virtual bool assembleDiscreteTerms(const IntegrandBase* itg,
                                     const TimeDomain& time);

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

  unsigned char nsv; //!< Number of consequtive solution vectors in core

  bool modal; //!< Modal dynamics simulation flag
};

#endif
