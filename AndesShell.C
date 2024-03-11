// $Id$
//==============================================================================
//!
//! \file AndesShell.C
//!
//! \brief Class representing a linear elastic shell with rotational DOFs.
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \date Feb 15 2024
//!
//==============================================================================

#include "AndesShell.h"
#include "ASMu2DNastran.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include "tinyxml2.h"


#ifdef HAS_ANDES
extern "C" {
  //! \brief Interface to ANDES 3-noded shell element routine (FORTRAN-90 code).
  void ifem_andes3_(const int& iel, const double* X0, const double& Thick,
                    const double& Emod, const double& Rny, const double& Rho,
                    const double* Press, double* Ekt, double* Em, double* Es,
                    int& iERR);
  //! \brief Interface to ANDES 4-noded shell element routine (FORTRAN-90 code).
  void ifem_andes4_(const int& iel, const double* X0, const double& Thick,
                    const double& Emod, const double& Rny, const double& Rho,
                    double* Ekt, double* Em, int& iERR);
}
#endif


AndesShell::AndesShell (unsigned short int n)
{
  nsd = 3; // Number of spatial dimenstions
  npv = 6; // Number of primary unknowns per node
  nSV = n; // Number of solution vectors in core

  // Default material properties
  Emod  = 2.1e11;
  Rny   = 0.3;
  Thick = 0.1;

  currentPatch = nullptr;
}


void AndesShell::printLog () const
{
  IFEM::cout <<"Formulation: ANDES shell"<< std::endl;
}


void AndesShell::initForPatch (const ASMbase* pch)
{
  currentPatch = dynamic_cast<const ASMu2DNastran*>(pch);
}


LocalIntegral* AndesShell::getLocalIntegral (size_t nen, size_t, bool) const
{
  ElmMats* result = new ElmMats();

  switch (m_mode)
  {
    case SIM::STATIC:
    case SIM::MASS_ONLY:
      result->resize(1,1);
      break;

    case SIM::VIBRATION:
      result->resize(2,0);
      break;

    case SIM::STIFF_ONLY:
      result->resize(1,0);
      break;

    case SIM::RHS_ONLY:
      result->resize(0,1);
      result->rhsOnly = true;
      result->withLHS = false;
      break;

    default:
      ;
  }

  result->redim(6*nen);
  return result;
}


bool AndesShell::initElement (const std::vector<int>& MNPC,
                              const FiniteElement& fe, const Vec3&, size_t,
                              LocalIntegral& elmInt)
{
  if (!this->initElement(MNPC,elmInt))
    return false;
  else if (fe.Xn.cols() == 1)
    return true;
  else if (currentPatch)
    return currentPatch->getProps(fe.iel,Emod,Rny,Rho,Thick);
  else
    return false;
}


bool AndesShell::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                          const Vec3& X) const
{
  if (eS <= 0) return true; // No external load vector

  if (Rho <= 0.0 || gravity.isZero()) return true; // No gravity load

  // Equvivalent pressure load due to gravity
  Vec3 p = gravity*(Rho*Thick);

  // Integrate the external load vector
  Vector& Svec = static_cast<ElmMats&>(elmInt).b[eS-1];
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (unsigned short int i = 1; i <= 3; i++)
      Svec(npv*(a-1)+i) += fe.N(a)*p(i)*fe.detJxW;

  return true;
}


bool AndesShell::finalizeElement (LocalIntegral& elmInt,
                                  const FiniteElement& fe,
                                  const TimeDomain&, size_t)
{
  if (eKm <= 0) return false;

  int iERR = 0;
#ifdef HAS_ANDES
  // Invoke the Fortran wrapper
  Matrix& Kmat = static_cast<ElmMats&>(elmInt).A[eKm-1];
  Matrix& Mmat = static_cast<ElmMats&>(elmInt).A[eM > 0 ? eM-1 : eKm-1];
  size_t nenod = fe.Xn.cols();
  if (nenod == 1) // 1-noded concentrated mass element
  {
    if (currentPatch && eM > 0)
      currentPatch->getMassMatrix(fe.iel, Mmat);
    if (currentPatch && eS > 0)
      currentPatch->getLoadVector(fe.iel, gravity,
                                  static_cast<ElmMats&>(elmInt).b[eS-1]);
  }
  else if (nenod == 3) // 3-noded shell element
  {
    Matrix Press(3,nenod);
    if (eS > 0 && Rho > 0.0 && !gravity.isZero())
    {
      // Equvivalent pressure load due to gravity
      Vec3 p = gravity*(Rho*Thick);
      for (size_t i = 1; i <= nenod; i++)
        Press.fillColumn(i,p.ptr());
    }
    Vector vDummy;
    Vector& Svec = eS > 0 ? static_cast<ElmMats&>(elmInt).b[eS-1] : vDummy;
    ifem_andes3_(fe.iel, fe.Xn.ptr(), Thick, Emod, Rny, eM > 0 ? Rho : 0.0,
                 Press.ptr(), Kmat.ptr(), Mmat.ptr(), Svec.ptr(), iERR);
  }
  else if (nenod == 4) // 4-noded shell element
  {
    ifem_andes4_(fe.iel, fe.Xn.ptr(), Thick, Emod, Rny, eM > 0 ? Rho : 0.0,
                 Kmat.ptr(), Mmat.ptr(), iERR);
  }
  else
  {
    iERR = -98;
    std::cerr <<" *** AndesShell: Invalid element, nenod="<< nenod << std::endl;
  }
#else
  iERR = -99;
  std::cerr <<" *** AndesShell: Built without this element."<< std::endl;
#endif
  return iERR >= 0;
}
