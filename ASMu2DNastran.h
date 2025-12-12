// $Id$
//==============================================================================
//!
//! \file ASMu2DNastran.h
//!
//! \date Feb 27 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of 2D %Lagrange FE models from Nastran Bulk Data File (BDF).
//!
//==============================================================================

#ifndef _ASM_U2D_NASTRAN_H
#define _ASM_U2D_NASTRAN_H

#include "ASMu2DLag.h"
#include "ASMu1DLag.h"

class ElementBlock;
class BeamProperty;
class FFlElementBase;


/*!
  \brief Driver for assembly of unstructured 2D %Lagrange FE models.
  \details This class overrides the read() method of its parent class
  assuming the mesh data is stored in the Nastran Bulk Data File format.
*/

class ASMu2DNastran : public ASMu2DLag
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  ASMu2DNastran(unsigned char n, unsigned char n_f, const std::string& path);
  //! \brief Disable default copy constructor.
  ASMu2DNastran(const ASMu2DNastran&) = delete;
  //! \brief The destructor deletes the immersed/extra element blocks, if any.
  virtual ~ASMu2DNastran();

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream& is);

  //! \brief Retrieves the shell thickness for element \a iel (1-based index).
  bool getThickness(size_t iel, double& t) const;
  //! \brief Retrieves mass properties for element with external ID \a eId.
  bool getMassProp(int eId, size_t igel, double& rho, double& t) const;
  //! \brief Retrieves stiffness properties for element with external ID \a eId.
  bool getStiffProp(int eId, size_t igel, double& E, double& nu) const;
  //! \brief Retrieves the stiffness matrix for element with external ID \a eId.
  bool getStiffnessMatrix(int eId, Matrix& K) const;
  //! \brief Retrieves the mass matrix for element with external ID \a eId.
  bool getMassMatrix(int eId, Matrix& M) const;
  //! \brief Retrieves the load vector for element with external ID \a eId.
  bool getLoadVector(int eId, const Vec3& g, Vector& S) const;

  //! \brief Evaluates the surface pressure at current integration point.
  bool addPressureAt(Vec3& p, int eId, const RealArray& N) const;

  //! \brief Initializes the \ref elmPres member.
  bool initPressureCache();

  //! \brief Checks if an element is associated with a given element set or not.
  //! \param[in] iEl Global external element number (1-based)
  //! \param[in] idx GLobal element index (0-based)
  //! \param[in] iSet Element set index (1-based) for a pressure load
  bool checkPressSet(int iEl, size_t idx, int iSet) const;

  //! \brief Calculates the shell normal vectors and the element centers
  //! \param[out] normals Vector of element-center normal-vector pairs
  bool getShellNormals(std::vector<Vec3Pair>& normals) const;

  //! \brief Prints out additional element information.
  virtual void printElmInfo(int iel, const IntegrandBase* integr) const;

  //! \brief Returns an additional geometry to visualize (point masses, etc.).
  virtual ElementBlock* immersedGeometry(char* name) const;
  //! \brief Returns an additional geometry to visualize (constraints, etc.).
  virtual ElementBlock* extraGeometry(char* name) const;
  //! \brief Returns an additional geometry to visualize a sensor location.
  ElementBlock* sensorGeometry(int idx, bool nodal) const;

  using ASMu2DLag::evalSolution;
  //! \brief Evaluates the secondary solution field at all nodal points.
  //! \param[out] sField Solution field
  //! \param[in] integr Object with problem-specific data and methods
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integr,
                            const int*, char) const;
  //! \brief Evaluates the secondary solution field at the element centers.
  //! \param[out] sField Solution field
  //! \param[in] integr Object with problem-specific data and methods
  //! \param[in] atElmCenters If \e false, evaluate at all nodal points instead
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integr,
                            const RealArray*, bool atElmCenters) const;
  //! \brief Evaluates the secondary solution field at the element centers.
  //! \param[out] sField Solution field
  //! \param[in] integr Object with problem-specific data and methods
  //! \param[in] elements List of elements to evaluate at (all if empty)
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integr,
                            const IntVec& elements) const;

  //! \brief Evaluates the primary solution at immersed geometry points.
  //! \param[out] field Solution field values at immersed geometry points
  //! \param[in] locSol Solution vector local to current patch
  virtual bool immersedSolution(Matrix& field, const Vector& locSol) const;
  //! \brief Evaluates the primary solution at extra geometry points.
  //! \param[out] field Solution field values at extra geometry points
  //! \param[in] locSol Solution vector local to current patch
  virtual bool extraSolution(Matrix& field, const Vector& locSol) const;

  //! \brief Returns \e true if this patch has element-wise surface pressures.
  bool haveLoads() const { return !myLoads.empty(); }

  //! \brief Returns the sub-patch with beam elements, if any.
  ASMbase* haveBeams() const { return beamPatch; }

protected:
  //! \brief Adds an element block with additional geometry.
  void addBlock(int idx, ElementBlock* blk);

  //! \brief Adds a beam element to this patch.
  void addBeamElement(FFlElementBase* elm, int eId, const IntVec& mnpc,
                      IntMat& beamMNPC, IntVec& beamElms, IntVec& beamNodes,
                      int& nErr, bool useEcc = true);
  //! \brief Adds a shell element to this patch.
  bool addShellElement(FFlElementBase* elm, int eId, const IntVec& mnpc);
  //! \brief Adds a mass element to this patch.
  void addMassElement(FFlElementBase* elm, int eId, int inod);
  //! \brief Adds a spring element to this patch.
  void addSpringElement(FFlElementBase* elm, int eId, const IntVec& mnpc,
                        int& nErr);

  //! \brief Adds MPCs representing a set of flexible couplings to this patch.
  void addFlexibleCouplings(FFlElementBase* elm, int eId, const IntVec& mnpc);
  //! \brief Adds one MPC representing a flexible coupling to this patch.
  void addFlexibleCoupling(int eId, int lDof, const int* indC,
                           const RealArray& weights, double* work,
                           const IntVec& mnpc, const Matrix& Xnod);

  //! \brief Returns the index into \ref myProps for element \a eId.
  size_t getPropIndex(int eId, size_t igel, char label) const;


  //! \brief Evaluates the secondary solution field at element centers.
  //! \param[out] sField Solution field
  //! \param[in] integr Object with problem-specific data and methods
  //! \param[in] atNodes If \e true, evaluate at all nodal points instead
  //! \param[in] elements List of elements to evaluate at (all if empty)
  bool evalSecSolution(Matrix& sField, const IntegrandBase& integr,
                       bool atNodes, const IntVec& elements = {}) const;

public:
  //! \brief Data type for shell element properties.
  struct ShellProps
  {
    int    id1   = -1;     //!< Material property ID
    int    id2   = -1;     //!< Shell thickness property ID
    double Thick = 0.1;    //!< Shell thickness
    double Emod  = 2.1e11; //!< Young's modulus
    double Rny   = 0.3;    //!< Poisson's ratio
    double Rho   = 7850.0; //!< Mass density

    //! \brief Equality operator.
    bool operator==(const ShellProps& b)
    {
      return id1 == b.id1 && id2 == b.id2;
    }
  };

  //! \brief Data type for beam element properties.
  struct BeamProps
  {
    double Emod  = 2.1e11; //!< Young's modulus
    double Gmod  = 8.0e9;  //!< Shear modulus
    double Rho   = 7850.0; //!< Mass density

    Vec3 Zaxis; //!< Vector defining the local Z-axis of the element

    std::array<Vec3,2> eccN{}; //!< Nodal eccentricity vectors
    std::array<double,9> cs{}; //!< Cross section parameters
  };

private:
  std::vector<ShellProps>  myProps;  //!< Shell element property container
  std::map<int,BeamProps>  myBprops; //!< Beam element property container
  std::map<int,Matrix>     myStiff;  //!< Mass-less spring elements
  std::map<int,Matrix>     myMass;   //!< Concentrated mass elements
  std::map<int,Vec3Vec>    myLoads;  //!< Surface pressures

  double massMax; //!< The largets point mass in the model (for scaling)
  IntMat spiders; //!< Constraint element topologies

  std::vector<unsigned char> elmProp; //!< Element property index cache
  std::vector<unsigned char> elmPres; //!< Element pressure index cache

  using ElementBlockID = std::pair<int,ElementBlock*>; //!< Convenience type
  std::vector<ElementBlockID> myBlocks; //!< Geometries for masses and spiders
  std::array<size_t,2> nGnod; //!< Total number of additional geometry nodes

  ASMu1DLag* beamPatch; //!< Separate patch for beam elements
};


/*!
  \brief Driver for assembly of beam elements.
*/

class ASMuBeam : public ASMu1DLag
{
  using BeamProps = ASMu2DNastran::BeamProps; //!< Convenience alias

public:
  //! \brief The constructor initializes the mesh data.
  ASMuBeam(const Vec3Vec& coord, const IntMat& mmnpc,
           const IntVec& mlgn, const IntVec& mlge,
           const std::map<int,BeamProps>& props,
           unsigned char n, unsigned char n_f);
  //! \brief Disable default copy constructor.
  ASMuBeam(const ASMuBeam&) = delete;
  //! \brief Empty destructor.
  virtual ~ASMuBeam() {}

  //! \brief Sets \ref nGauss to 1 due to explicit matrices.
  virtual void setGauss(int) { nGauss = 1; }

  //! \brief Retrieves the properties for element with index \a id.
  bool getProps(int eId, double& E, double& G, double& rho,
                BeamProperty& bprop) const;

protected:
  //! \brief Initializes the local element axes for a patch of beam elements.
  virtual bool initLocalElementAxes(const Vec3& Zaxis);

private:
  //! Reference to the beam element property container
  const std::map<int,BeamProps>& myProps;
};

#endif
