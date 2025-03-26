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
  ASMu2DNastran(unsigned char n, unsigned char n_f, bool ns, char b)
    : ASMu2DLag(n,n_f,'N'), useBeams(b), useSets(!ns), beamPatch(nullptr) {}
  //! \brief Disable default copy constructor.
  ASMu2DNastran(const ASMu2DNastran&) = delete;
  //! \brief Empty destructor.
  virtual ~ASMu2DNastran() {}

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream& is);

  //! \brief Retrieves the shell thickness for element with external ID \a eId.
  bool getThickness(int eId, double& t) const;
  //! \brief Retrieves the properties for element with external ID \a eId.
  bool getProps(int eId, double& E, double& nu, double& rho, double& t) const;
  //! \brief Retrieves the stiffness matrix for element with external ID \a eId.
  bool getStiffnessMatrix(int eId, Matrix& K) const;
  //! \brief Retrieves the mass matrix for element with external ID \a eId.
  bool getMassMatrix(int eId, Matrix& M) const;
  //! \brief Retrieves the load vector for element with external ID \a eId.
  bool getLoadVector(int eId, const Vec3& g, Vector& S) const;

  //! \brief Evaluates the surface pressure at current integration point.
  bool addPressureAt(Vec3& p, int eId, const RealArray& N) const;

  //! \brief Evaluates the secondary solution field at all nodal points.
  //! \param[out] sField Solution field
  //! \param[in] integr Object with problem-specific data and methods
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integr,
                            const int*, char) const;
  //! \brief Evaluates the secondary solution field at all element centers.
  //! \param[out] sField Solution field
  //! \param[in] integr Object with problem-specific data and methods
  //! \param[in] atElmCenters If \e false, evaluate at nodal points instead
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integr,
                            const RealArray*, bool atElmCenters) const;

  //! \brief Checks if an external element ID is within a predefined set.
  //! \todo Maybe put this in ASMbase later, or extend isInElementSet()
  //! to handle both internal indices and external element numbers.
  bool isElementInSet(int elmId, int iset) const
  {
    return this->isInElementSet(iset,this->getElmIndex(elmId));
  }

  //! \brief Returns \e true if this patch has element-wise surface pressures.
  bool haveLoads() const { return !myLoads.empty(); }

  //! \brief Returns the sub-patch with beam elements, if any.
  ASMbase* haveBeams() const { return beamPatch; }

protected:
  //! \brief Adds a beam element to this patch.
  void addBeamElement(FFlElementBase* elm, int eId, const IntVec& mnpc,
                      IntMat& beamMNPC, IntVec& beamElms, IntVec& beamNodes,
                      int& nErr, bool useEcc = true);
  //! \brief Adds a shell element to this patch.
  void addShellElement(FFlElementBase* elm, int eId, const IntVec& mnpc);
  //! \brief Adds a mass element to this patch.
  void addMassElement(FFlElementBase* elm, int eId, int inod);
  //! \brief Adds a spring element to this patch.
  void addSpringElement(FFlElementBase* elm, int eId, const IntVec& mnpc,
                        int& nErr);
  //! \brief Adds MPCs representing a flexible coupling to this patch.
  void addFlexibleCouplings(FFlElementBase* elm, int eId, const IntVec& mnpc);
  //! \brief Adds MPCs representing a flexible coupling to this patch.
  void addFlexibleCoupling(int eId, int lDof, const int* indC,
                           const RealArray& weights,
                           const IntVec& mnpc, const Matrix& Xnod);

public:
  //! \brief Data type for shell element properties.
  struct ShellProps
  {
    double Thick = 0.1;    //!< Shell thickness
    double Emod  = 2.1e11; //!< Young's modulus
    double Rny   = 0.3;    //!< Poisson's ratio
    double Rho   = 7850.0; //!< Mass density
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
  std::map<int,ShellProps> myProps;  //!< Shell element property container
  std::map<int,BeamProps>  myBprops; //!< Beam element property container
  std::map<int,Matrix>     myStiff;  //!< Mass-less spring elements
  std::map<int,Matrix>     myMass;   //!< Concentrated mass elements
  std::map<int,Vec3Vec>    myLoads;  //!< Surface pressures

  char useBeams; //!< If nonzero include beam elements as a separate patch
  bool useSets; //!< If \e true, read Nastran SET definitions

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
