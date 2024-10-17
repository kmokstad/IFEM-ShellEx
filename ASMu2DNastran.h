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


/*!
  \brief Driver for assembly of unstructured 2D %Lagrange FE models.
  \details This class overrides the read() method of its parent class
  assuming the mesh data is stored in the Nastran Bulk Data File format.
*/

class ASMu2DNastran : public ASMu2DLag
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  ASMu2DNastran(unsigned char n, unsigned char n_f, char fType = 'N')
    : ASMu2DLag(n,n_f,fType) {}
  //! \brief Disable default copy constructor.
  ASMu2DNastran(const ASMu2DNastran& pch) = delete;
  //! \brief Empty destructor.
  virtual ~ASMu2DNastran() {}

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream& is);

  //! \brief Retrieves the shell thickness for element with external ID \a eId.
  bool getThickness(int eId, double& t) const;
  //! \brief Retrieves the properties for element with external ID \a eId.
  bool getProps(int eId, double& E, double& nu, double& rho, double& t) const;
  //! \brief Retrieves the mass matrix for element with with external ID \a eId.
  bool getMassMatrix(int eId, Matrix& eM) const;
  //! \brief Retrieves the load vector for element with with external ID \a eId.
  bool getLoadVector(int eId, const Vec3& g, Vector& eS) const;

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
  //! \TODO Maybe put this in ASMbase later, or extend isInElementSet()
  //! to handle both internal indices and external element numbers.
  bool isElementInSet(int elmId, int idx) const
  {
    return this->isInElementSet(idx,this->getElmIndex(elmId));
  }

protected:
  //! \brief Adds MPCs representing a flexible coupling to this patch.
  void addFlexibleCoupling(int eId, int lDof, const int* indC,
                           const std::vector<double>& weights,
                           const IntVec& mnpc, const Matrix& Xnod);

private:
  //! \brief Data type for shell element properties.
  struct ShellProps
  {
    double Thick = 0.1;    //!< Shell thickness
    double Emod  = 2.1e11; //!< Young's modulus
    double Rny   = 0.3;    //!< Poisson's ratio
    double Rho   = 7850.0; //!< Mass density
  };

  std::map<int,ShellProps> myProps;  //!< Shell element property container
  std::map<int,Matrix>     myMass;   //!< Concentrated mass elements
  std::map<int,Vec3Vec>    myLoads;  //!< Surface pressures
};

#endif
