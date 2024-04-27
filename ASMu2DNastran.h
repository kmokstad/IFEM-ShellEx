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


/*!
  \brief Driver for assembly of unstructured 2D %Lagrange FE models.
  \details This class overrides the read() method of its parent class
  assuming the mesh data is stored in the Nastran Bulk Data File format.
*/

class ASMu2DNastran : public ASMu2DLag
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  ASMu2DNastran(unsigned char n, unsigned char n_f)
    : ASMu2DLag(n,n_f,'N') { beamPatch = nullptr; }
  //! \brief Disable default copy constructor.
  ASMu2DNastran(const ASMu2DNastran&) = delete;
  //! \brief Empty destructor.
  virtual ~ASMu2DNastran() {}

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream& is);

  //! \brief Retrieves the properties for element with index \a eId.
  bool getProps(int eId, double& E, double& nu, double& rho, double& t) const;
  //! \brief Retrieves the mass matrix for element with index \a id.
  bool getMassMatrix(int eId, Matrix& eM) const;
  //! \brief Retrieves the load vector for mass element with index \a id.
  bool getLoadVector(int eId, const Vec3& g, Vector& eS) const;

  //! \brief Evaluates the surface pressure at current integration point.
  Vec3 getPressureAt(int iel, const RealArray& N) const;

  //! \brief Returns the sub-patch with beam elements, if any.
  ASMbase* haveBeams() const { return beamPatch; }

protected:
  //! \brief Adds MPCs representing a flexible coupling to this patch.
  void addFlexibleCoupling(int iel, int lDof, const int* indC,
                           const std::vector<double>& weights,
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
  std::map<int,Matrix>     myMass;   //!< Concentrated mass elements
  std::map<int,Vec3Vec>    myLoads;  //!< Surface pressures

  ASMu1DLag* beamPatch; //!< Separate patch for beam elements
};


/*!
  \brief Driver for assembly of beam elements.
*/

class ASMuBeam : public ASMu1DLag
{
public:
  //! \brief The constructor initializes the mesh data.
  ASMuBeam(const Vec3Vec& coord, const IntMat& mmnpc,
           const IntVec& mlgn, const IntVec& mlge,
           const std::map<int,ASMu2DNastran::BeamProps>& props,
           unsigned char n, unsigned char n_f);
  //! \brief Disable default copy constructor.
  ASMuBeam(const ASMuBeam&) = delete;
  //! \brief Empty destructor.
  virtual ~ASMuBeam() {}

  //! \brief Retrieves the properties for element with index \a id.
  bool getProps(int eId, double& E, double& G, double& rho,
                BeamProperty& bprop) const;

protected:
  //! \brief Initializes the local element axes for a patch of beam elements.
  virtual bool initLocalElementAxes(const Vec3&);

private:
  //! Reference to the beam element property container
  const std::map<int,ASMu2DNastran::BeamProps>& myProps;
};

#endif
