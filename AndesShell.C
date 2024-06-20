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
#include "NewmarkMats.h"
#include "TimeDomain.h"
#include "Function.h"
#include "Vec3Oper.h"
#include "Utilities.h"
#include "VTF.h"
#include "IFEM.h"
#include "tinyxml2.h"
#include <fstream>


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
  //! \brief Interface to 4-noded shell stress routine (FORTRAN-90 code).
  void ifem_strs24_(const int& iel, const double* X0, const double& Thick,
                    const double& Emod, const double& Rny, const double* Ev,
                    double* SR, double* Sigma, int& iERR);
}
#endif


AndesShell::AndesShell (unsigned short int n, bool modal)
{
  nsd = 3; // Number of spatial dimenstions
  npv = 6; // Number of primary unknowns per node
  nSV = n; // Number of solution vectors in core

  // Default material properties
  Emod  = 2.1e11;
  Rny   = 0.3;
  Thick = 0.1;
  Rho   = 7.85e3;
  ovrMat = false;

  trInside = trOutside = 0.0;

  isModal = modal;

  currentPatch = nullptr;
}


AndesShell::~AndesShell ()
{
  if (!degenerated.empty())
  {
    std::ofstream os("degenerated_elements.ftl");
    os <<"GROUP{1"; for (int e : degenerated) os <<" "<< e;
    os <<" {NAME \"degenerated elements\"}}\n";
  }

  if (!straightline.empty())
  {
    std::ofstream os("straightline_elements.ftl");
    os <<"GROUP{2"; for (int e : straightline) os <<" "<< e;
    os <<" {NAME \"straight line elements\"}}\n";
  }
}


Material* AndesShell::parseMatProp (const tinyxml2::XMLElement* elem, bool)
{
  if (utl::getAttribute(elem,"override",ovrMat) && ovrMat)
  {
    IFEM::cout <<"\tPatch-level material properties are overridden:";
    if (utl::getAttribute(elem,"E",Emod))
      IFEM::cout <<" E="<< Emod;
    if (utl::getAttribute(elem,"nu",Rny))
      IFEM::cout <<" nu="<< Rny;
    if (utl::getAttribute(elem,"rho",Rho))
      IFEM::cout <<" rho="<< Rho;
    IFEM::cout << std::endl;
  }

  const tinyxml2::XMLElement* child = elem->FirstChildElement("thickloss");
  if (child && child->FirstChild())
  {
    utl::getAttribute(child,"t1",trOutside);
    utl::getAttribute(child,"t2",trInside);
    std::istringstream(child->FirstChild()->Value()) >> Xlow >> Xupp;
    IFEM::cout <<"\tThickness loss: t1="<< trOutside <<" t2="<< trInside
               <<"  Xlower = "<< Xlow <<"  Xupper = "<< Xupp << std::endl;
  }

  return nullptr;
}


void AndesShell::printLog () const
{
  IFEM::cout <<"Formulation: ANDES shell";
  IFEM::cout << std::endl;
}


void AndesShell::setMode (SIM::SolutionMode mode)
{
  if (isModal && mode == SIM::DYNAMIC)
    mode = SIM::RHS_ONLY;

  this->ElasticBase::setMode(mode);

  if (mode == SIM::STATIC || (isModal && mode == SIM::RHS_ONLY)) iS  = 0;
}


void AndesShell::initLHSbuffers (size_t nEl)
{
  if (nEl > 1)
  {
    myKmats.resize(nEl);
    myMmats.resize(nEl);
  }
  else if (nEl == 0 && !myKmats.empty())
  {
    if (eKm > 0) eKm = -eKm;
    if (eM  > 0) eM  = -eM;
  }
}


void AndesShell::initForPatch (const ASMbase* pch)
{
  currentPatch = dynamic_cast<const ASMu2DNastran*>(pch);
}


bool AndesShell::setPressure (RealFunc* pf, int code,
                              const std::string& sName, const ASMbase* pch)
{
  size_t sIdx = 0;
  if (sName.empty())
    presFld[-code] = pf;
  else if ((sIdx = static_cast<const ASMu2DLag*>(pch)->getElementSetIdx(sName)))
    presFld[sIdx] = pf;
  else
    return false;

  return true;
}


void AndesShell::initIntegration (size_t nGp, size_t)
{
  presVal.clear();
  if (this->havePressure())
    presVal.resize(nGp,{Vec3(),Vec3()});
}


LocalIntegral* AndesShell::getLocalIntegral (size_t nen, size_t, bool) const
{
  ElmMats* result = nullptr;

  if (!isModal && m_mode == SIM::DYNAMIC)
    result = new NewmarkMats(intPrm[0],intPrm[1],intPrm[2],intPrm[3]);
  else
    result = new ElmMats();

  switch (m_mode)
  {
    case SIM::STATIC:
    case SIM::MASS_ONLY:
      result->resize(1,1);
      break;

    case SIM::DYNAMIC:
      result->resize(3,1);
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


int AndesShell::getIntegrandType () const
{
  return trInside > 0.0 && trOutside > 0.0 ? ELEMENT_CENTER : STANDARD;
}


/*!
  In case that element stiffness- and mass-matrices are cached, this method
  will copy the element matrices for current element from the cache
  and skip the property initialisation which then are not needed.
*/

bool AndesShell::initElement (const std::vector<int>& MNPC,
                              const FiniteElement& fe, const Vec3& Xc, size_t,
                              LocalIntegral& elmInt)
{
  if (fe.iel > 0)
  {
    size_t iel = fe.iel - 1;
    if (iel < myKmats.size() && eKm < 0)
      static_cast<ElmMats&>(elmInt).A[-eKm-1] = myKmats[iel];
    if (iel < myMmats.size() && eM  < 0)
      static_cast<ElmMats&>(elmInt).A[-eM-1]  = myMmats[iel];
  }

  if (!this->initElement(MNPC,elmInt))
    return false;
  else if (fe.Xn.cols() == 1 || (eKm+eM <= 0 && fe.Xn.cols() != 3))
    return true;
  else if (!currentPatch)
    return false;

  bool ok = true;
  if (ovrMat) // Override the patch-level material properties
  {
    double dum1, dum2, dum3;
    ok = currentPatch->getProps(fe.iel,dum1,dum2,dum3,Thick);
  }
  else // Use patch-level material properties
    ok = currentPatch->getProps(fe.iel,Emod,Rny,Rho,Thick);

  // Scale the thickness depending on location inside or outside given box
  if (trInside > 0.0 && trOutside > 0.0)
    Thick *= 1.0 - (Xlow < Xc && Xc < Xupp ? trInside : trOutside);

  return ok;
}


/*!
  This method only deals with the element load vector due to pressure and
  gravity loads on the 4-noded shell element. The element stiffness- and mass-
  matrices, as well as the load vector for the 3-noded shell, are all calculated
  in the finalizeElement() method by invoking the Fortran wrapper.
*/

bool AndesShell::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                          const Vec3& X) const
{
  if (eS <= 0) return true; // No external load vector

  Vec3 p, n;
  bool havePressure = Rho > 0.0 && !gravity.isZero();
  if (havePressure)
    p = gravity*(Rho*Thick); // Equivalent pressure load due to gravity

  if (fe.G.cols() >= 2 && this->havePressure(fe.iel))
  {
    // Shell normal vector
    n.cross(fe.G.getColumn(1),fe.G.getColumn(2));
    n.normalize();
    // Evaluate the pressure at this point
    this->addPressure(p,X,n,fe.iel);
    havePressure = true;
  }

  // Add surface loads from the FE model, if any
  if (currentPatch && currentPatch->addPressureAt(p,fe.iel,fe.N))
    havePressure = true;

  // Store pressure value for visualization
  if (havePressure && fe.iGP < presVal.size())
    presVal[fe.iGP] = { X, p };

  if (p.isZero()) return true; // No pressure load

  // Integrate the external load vector
  Vector& Svec = static_cast<ElmMats&>(elmInt).b[eS-1];
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (unsigned short int i = 1; i <= 3; i++)
      Svec(npv*(a-1)+i) += fe.N(a)*p(i)*fe.detJxW;

  return true;
}


/*!
  This method constructs the element stiffness and mass matrices for the 3- and
  4-noded ANDES shell elements, by invoking the appropriate Fortran wrapper.
  The internal forces are also calculated from the element stiffness matrix and
  the current element displacement vector stored in \a elmInt.

  The calculated element matrices are cached in the internal buffers
  \ref myKmats and \ref myMmats, if those buffers have been allocated
  by invoking initLHSbuffers() with \a nEl > 1 as argument.
*/

bool AndesShell::finalizeElement (LocalIntegral& elmInt,
                                  const FiniteElement& fe,
                                  const TimeDomain& time, size_t)
{
  int iERR = -99;
  Vector vDummy(1);
  Matrix mDummy(1,1);
  Vector& Svec = eS  > 0 ? static_cast<ElmMats&>(elmInt).b[eS-1]  : vDummy;
  Matrix& Kmat = eKm > 0 ? static_cast<ElmMats&>(elmInt).A[eKm-1] : mDummy;
  Matrix& Mmat = eM  > 0 ? static_cast<ElmMats&>(elmInt).A[eM-1 ] : mDummy;
  size_t nenod = fe.Xn.cols();
  if (currentPatch && nenod == 1) // 1-noded concentrated mass element
  {
    iERR = 0;
    if (eKm > 0)
      Kmat.clear(); // no stiffness
    if (eM > 0)
      currentPatch->getMassMatrix(fe.iel, Mmat);
    if (eS > 0)
      currentPatch->getLoadVector(fe.iel, gravity, Svec);
  }
  else if (nenod == 3) // 3-noded shell element
  {
    Matrix Press(3,nenod);
    if (eS > 0 && Rho > 0.0 && !gravity.isZero())
    {
      // Equivalent pressure load due to gravity
      Vec3 p = gravity*(Rho*Thick);
      for (size_t i = 1; i <= nenod; i++)
        Press.fillColumn(i,p.ptr());
    }
    if (eS > 0 && this->havePressure(fe.iel))
    {
      // Find the shell normal vector
      std::array<Vec3,3> X = { fe.Xn.ptr(0), fe.Xn.ptr(1), fe.Xn.ptr(2) };
      Vec3 n(X[1]-X[0],X[2]-X[0]);
      n.normalize();
      // Calculate nodal pressure intensities
      for (size_t i = 1; i <= 3; i++)
      {
        Vec3 p;
        this->addPressure(p,Vec4(X[i-1],time.t),n,fe.iel);
        for (size_t j = 1; j <= 3; j++)
          Press(j,i) += p(j);
      }
    }
    if (eS > 0 && currentPatch) // Add surface loads from the FE model, if any
      for (size_t i = 1; i <= 3; i++)
      {
        Vec3 p;
        if (currentPatch->addPressureAt(p,fe.iel,{-static_cast<double>(i)}))
          for (size_t j = 1; j <= 3; j++)
            Press(j,i) += p(j);
      }
#ifdef HAS_ANDES
    // Invoke Fortran wrapper for the 3-noded ANDES element
    ifem_andes3_(fe.iel, fe.Xn.ptr(),
                 eKm > 0 ? Thick : 0.0, Emod, Rny, eM > 0 ? Rho : 0.0,
                 Press.ptr(), Kmat.ptr(), Mmat.ptr(), Svec.ptr(), iERR);
#endif
  }
  else if (nenod == 4) // 4-noded shell element
  {
#ifdef HAS_ANDES
    // Invoke Fortran wrapper for the 4-noded ANDES element
    ifem_andes4_(fe.iel, fe.Xn.ptr(),
                 eKm > 0 ? Thick : 0.0, Emod, Rny, eM > 0 ? Rho : 0.0,
                 Kmat.ptr(), Mmat.ptr(), iERR);
#endif
  }
  else
  {
    iERR = -98;
    std::cerr <<" *** AndesShell: Invalid element, nenod="<< nenod << std::endl;
  }

  if (iERR >= 1 && iERR <= 3)
    degenerated.insert(fe.iel);
  else if (iERR == 4)
    straightline.insert(fe.iel);
  else if (iERR == -99)
    std::cerr <<" *** AndesShell: Built without this element."<< std::endl;

  if (fe.iel > 0 && iERR >= 0)
  {
    size_t iel = fe.iel - 1;
    if (iel < myKmats.size() && Kmat.size() > 1) myKmats[iel] = Kmat;
    if (iel < myMmats.size() && Mmat.size() > 1) myMmats[iel] = Mmat;
  }

  if (iS > 0 && !elmInt.vec.empty() && nenod > 1 && iERR >= 0)
  {
    size_t iel = fe.iel - 1;
    Matrix& Sm = fe.iel > 0 && iel < myKmats.size() ? myKmats[iel] : Kmat;
    Vector& Sv = static_cast<ElmMats&>(elmInt).b[iS-1];
    if (!Sm.multiply(elmInt.vec.front(),Sv,false,-1))
      iERR = -97;
  }

  return iERR >= 0 && this->finalizeElement(elmInt,time);
}


bool AndesShell::evalSol2 (Vector& s, const Vectors& eV,
                           const FiniteElement& fe, const Vec3&) const
{
  int iERR = 0;
  size_t nenod = fe.Xn.cols();
  if (nenod == 1) // 1-noded concentrated mass element (no stresses)
    s.clear();
  else if (nenod == 3) // 3-noded shell element
    s.clear();
  else if (nenod == 4) // 4-noded shell element
  {
    s.resize(18,0.0);
#ifdef HAS_ANDES
    // Invoke Fortran wrapper for the 4-noded ANDES element
    ifem_strs24_(fe.iel, fe.Xn.ptr(), Thick, Emod, Rny,
                 eV.front().data(), s.data(), s.data()+6, iERR);
#endif
  }
  else
  {
    iERR = -98;
    std::cerr <<" *** AndesShell: Invalid element, nenod="<< nenod << std::endl;
  }

  if (s.size() >= 18)
  {
    // Calculate principal and von Mises stresses at the top and bottom surfaces
    Vec3 sigma_p;
    size_t j = 12;
    for (size_t i = 9; i >= 6; i -= 3, j -= 6)
    {
      SymmTensor sigma({s[i],s[i+1],s[i+2]});
      sigma.principal(sigma_p);
      for (size_t k = 0; k < 3 && i > 6; k++)
        s[j+k] = s[i+k];
      s[j+3] = sigma_p.x;
      s[j+4] = sigma_p.y;
      s[j+5] = sigma.vonMises();
    }
  }

  return iERR == 0;
}


bool AndesShell::havePressure (int iel) const
{
  for (const std::pair<const int,RealFunc*>& press : presFld)
    if (press.first < 0 || iel < 0)
      return true;
    else if (!currentPatch)
      break;
    else
    {
      const IntVec& eSet = currentPatch->getElementSet(press.first);
      if (std::find(eSet.begin(),eSet.end(),iel) != eSet.end())
        return true;
    }

  return false;
}


void AndesShell::addPressure (Vec3& p, const Vec3& X,
                              const Vec3& n, int iel) const
{
  for (const std::pair<const int,RealFunc*>& press : presFld)
    if (press.first < 0)
      p += (*press.second)(X)*n;
    else
    {
      const IntVec& eSet = currentPatch->getElementSet(press.first);
      if (std::find(eSet.begin(),eSet.end(),iel) != eSet.end())
        p += (*press.second)(X)*n;
    }
}


bool AndesShell::hasTractionValues () const
{
  return !presVal.empty();
}


bool AndesShell::writeGlvT (VTF* vtf, int iStep,
                            int& geoBlk, int& nBlock) const
{
  if (presVal.empty())
    return true;
  else if (!vtf)
    return false;

  // Write surface pressures as discrete point vectors to the VTF-file
  return vtf->writeVectors(presVal,geoBlk,++nBlock,"Pressure",iStep);
}


std::string AndesShell::getField2Name (size_t i, const char* prefix) const
{
  static const char* s[12] = { "n_xx", "n_yy", "n_xy", "m_xx", "m_yy", "m_xy",
                               "sigma_x", "sigma_y", "tau_xy",
                               "sigma_1", "sigma_2", "sigma_m" };

  std::string name(s[i < 6 ? i : 6 + i%6]);
  if (i >= 12)
    name = "Top " + name;
  else if (i >= 6)
    name = "Bottom " + name;

  if (!prefix)
    return name;

  return prefix + std::string(" ") + name;
}
