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
  explicit SIMAndesShell(short int n = 0, bool m = false);
  //! \brief The destructor deletes the nodal point load functions.
  virtual ~SIMAndesShell();

  //! \brief Enables the element matrix cache for reaction force calculation.
  void initForSingleStep() { if (!myRFset.empty()) this->initLHSbuffers(); }
  //! \brief Switches off separate reaction force calculations when multi-steps.
  virtual void initForMultiStep() { myRFset.clear(); }

  using SIMElasticity<SIM2D>::solveSystem;
  //! \brief Solves the assembled linear system of equations for a given load.
  //! \param[out] solution Global primary solution vector
  //! \param[in] printSol Print solution if its size is less than \a printSol
  //! \param[out] rCond Reciprocal condition number
  //!
  //! \details This method is overridden to also compute nodal reaction forces.
  //! This requires an additional assembly loop calculating the internal forces,
  //! since we only are doing a linear solve here.
  virtual bool solveSystem(Vector& solution, int printSol, double* rCond,
                           const char*, size_t);

  //! \brief Returns current reaction force container.
  virtual const RealArray* getReactionForces() const;

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

  //! \brief Renumbers the global node numbers of the springs and point loads.
  virtual bool renumberNodes(const std::map<int,int>& nodeMap);

  //! \brief Assembles the DOF springs and nodal point loads, if any.
  virtual bool assembleDiscreteTerms(const IntegrandBase* itg,
                                     const TimeDomain& time);

public:
  static char useBeams; //!< If non-zero, include beam elements, if any
  static bool readSets; //!< If \e true, also read Nastran SET definitions
  static bool replRBE3; //!< If \e true, replace all RBE3 elements by RBE2s

private:
  //! \brief Struct defining a DOF spring.
  struct DOFspring
  {
    int    inod;  //!< Node index
    int    ldof;  //!< Local DOF number
    double coeff; //!< Stiffness coefficient
    //! \brief Default constructor.
    explicit DOFspring(int n = 0, int d = 0, double c = 0.0)
      : inod(n), ldof(d), coeff(c) {}
  };

  //! \brief Struct defining a nodal point load.
  struct PointLoad
  {
    int         inod; //!< Node index
    int         ldof; //!< Local DOF number
    ScalarFunc* p;    //!< Load magnitude
    //! \brief Default constructor.
    explicit PointLoad(int n = 0, int d = 0, ScalarFunc* f = nullptr)
      : inod(n), ldof(d), p(f) {}
  };

  std::vector<DOFspring> mySprings; //!< Global DOF springs
  std::vector<PointLoad> myLoads;   //!< Nodal point loads

  RealArray   myReact; //!< Nodal reaction forces
  std::string myRFset; //!< Node set for calculation of reaction forces
  std::string myPath;  //!< Relative path of the patch file

  short int nss; //!< Number of consequtive solution states in core

  bool modal; //!< Modal dynamics simulation flag
};

#endif
