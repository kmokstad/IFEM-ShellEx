// $Id$
//==============================================================================
//!
//! \file ASMu2DMastran.h
//!
//! \date Feb 27 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of unstructured 2D Lagrange FE models from Nastran BDF-file.
//!
//==============================================================================

#ifndef _ASM_U2D_NASTRAN_H
#define _ASM_U2D_NASTRAN_H

#include "ASMu2DLag.h"


/*!
  \brief Driver for assembly of unstructured 2D Lagrange FE models.
  \details This class overrides the read() method of its parent class
  assuming the mesh data is stored in the Nastran Bulk Data Format.
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

  //! \brief Retrieves the properties for element with index \a id.
  bool getProps(int eId, double& E, double& nu, double& rho, double& t) const;
  //! \brief Retrieves the mass matrix for element with index \a id.
  bool getMassMatrix(int eId, Matrix& eM) const;
  //! \brief Retrieves the load vector for mass element with index \a id.
  bool getLoadVector(int eId, const Vec3& g, Vector& eS) const;

protected:
  void addFlexibleCoupling(int iel, int lDof, const int* indC,
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

  //! \brief Output stream operator.
  friend std::ostream& operator<<(std::ostream& os, const ShellProps& p)
  {
    return os <<" t="<< p.Thick <<" E="<< p.Emod <<" nu="<< p.Rny
              <<" rho="<< p.Rho;
  }

  std::map<int,ShellProps> myProps; //!< Shell element property container
  std::map<int,Matrix>     myMass;  //!< Concentrated mass elements
};

#endif
