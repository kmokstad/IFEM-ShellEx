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
                    double* Ekt, double* Em, int& iERR);
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
  IFEM::cout <<"Formulation: ANDES shell";
  IFEM::cout << std::endl;
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
      std::cerr <<" *** AndesShell::getLocalIntegral: Mode flag "<< m_mode
                <<" is not yet supported."<< std::endl;
      delete result;
      return nullptr;
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


bool AndesShell::finalizeElement (LocalIntegral& elmInt,
                                  const FiniteElement& fe,
                                  const TimeDomain&, size_t)
{
  int iERR = -99;
#ifdef HAS_ANDES
  Vector vDummy(1);
  Matrix mDummy(1,1);
  Vector& Svec = eS  > 0 ? static_cast<ElmMats&>(elmInt).b[eS-1]  : vDummy;
  Matrix& Kmat = eKm > 0 ? static_cast<ElmMats&>(elmInt).A[eKm-1] : mDummy;
  Matrix& Mmat = eM  > 0 ? static_cast<ElmMats&>(elmInt).A[eM-1 ] : mDummy;
  size_t nenod = fe.Xn.cols();
  if (currentPatch && nenod == 1) // 1-noded concentrated mass element
  {
    iERR = 0;
    if (eM > 0)
      currentPatch->getMassMatrix(fe.iel, Mmat);
    if (eS > 0)
      currentPatch->getLoadVector(fe.iel, gravity, Svec);
  }
  else if (nenod == 3)
    ifem_andes3_(fe.iel, fe.Xn.ptr(), Thick, Emod, Rny, eM > 0 ? Rho : 0.0,
                 Kmat.ptr(), Mmat.ptr(), iERR);
  else if (nenod == 4)
    ifem_andes4_(fe.iel, fe.Xn.ptr(), Thick, Emod, Rny, eM > 0 ? Rho : 0.0,
                 Kmat.ptr(), Mmat.ptr(), iERR);
  else
  {
    iERR = -98;
    std::cerr <<" *** AndesShell: Invalid element, nenod="<< nenod << std::endl;
  }
#endif
  if (iERR == -99)
    std::cerr <<" *** AndesShell: Built without this element."<< std::endl;
  return iERR >= 0;
}
