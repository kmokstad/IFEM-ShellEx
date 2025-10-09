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
#include "ElasticBeam.h"
#include "ASMu2DNastran.h"
#include "GlobalIntegral.h"
#include "FiniteElement.h"
#include "NewmarkMats.h"
#include "TimeDomain.h"
#include "Functions.h"
#include "Vec3Oper.h"
#include "Utilities.h"
#include "VTF.h"
#include "IFEM.h"
#include "tinyxml2.h"
#include <fstream>
#include <functional>


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


AndesShell::AndesShell (unsigned short int ns, bool modal, bool withBeams)
{
  nsd =  3; // Number of spatial dimensions
  npv =  6; // Number of primary unknowns per node
  n2v = 18; // Number of secondary variables for output
  nCS = ns; // Number of consecutive solution states in core
  nSV = ns; // Total number of solution vectors in core

  // Default material properties
  Emod  = 2.1e11;
  GorNu = 0.3;
  Thick = Thck0 = 0.1;
  Rho   = 7.85e3;
  ovrMat = false;

  trInside = trOutside = 0.0;
  thickLoss = nullptr;

  isModal = modal;

  currentPatch = nullptr;
  beamPatch = nullptr;
  myReacI = nullptr;

  if (withBeams)
  {
    beamProblem = new ElasticBeam(nCS);
    beamProblem->setProperty(&myBeamProps);
  }
  else
    beamProblem = nullptr;
}


AndesShell::~AndesShell ()
{
  delete thickLoss;
  delete beamProblem;

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

  if (!failedElements.empty())
  {
    std::ofstream os("failed_elements.ftl");
    os <<"GROUP{3"; for (int e : failedElements) os <<" "<< e;
    os <<" {NAME \"failed elements\"}}\n";
  }
}


void AndesShell::printLog () const
{
  IFEM::cout <<"Formulation: ANDES shell";
  if (beamProblem)
  {
    IFEM::cout <<" with beam elements\n";
    double E, G, Rho;
    if (beamPatch->getProps(beamPatch->getElmID(1),E,G,Rho,
                            const_cast<AndesShell*>(this)->myBeamProps))
    {
      ElasticBeam* beam = const_cast<AndesShell*>(this)->beamProblem;
      beam->setStiffness(E,G);
      beam->setMass(Rho);
    }
    beamProblem->printLog();
  }
  else
    IFEM::cout << std::endl;
}


bool AndesShell::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Parent()->Value(),"postprocessing"))
  {
    if (!strcasecmp(elem->Value(),"vonMises_only"))
      n2v = 2; // only output von Mises as secondary variables
  }
  else if (!strcasecmp(elem->Value(),"lumpedBeamMass") && beamProblem)
  {
    char version = 1;
    utl::getAttribute(elem,"version",version);
    beamProblem->useLumpedMass(version);
    if (version)
      IFEM::cout <<"\tUsing lumped mass matrix for beams, version "
                 << static_cast<int>(version) << std::endl;
  }

  else if (strcasecmp(elem->Value(),"material"))
    return this->ElasticBase::parse(elem);

  IFEM::cout <<"  Parsing <material>"<< std::endl;
  if (utl::getAttribute(elem,"override",ovrMat) && ovrMat)
  {
    IFEM::cout <<"\tPatch-level material properties are overridden:";
    if (utl::getAttribute(elem,"E",Emod))
      IFEM::cout <<" E="<< Emod;
    if (utl::getAttribute(elem,"nu",GorNu))
      IFEM::cout <<" nu="<< GorNu;
    if (utl::getAttribute(elem,"rho",Rho))
      IFEM::cout <<" rho="<< Rho;
    IFEM::cout << std::endl;
  }

  const tinyxml2::XMLElement* child = elem->FirstChildElement("thickness");
  const char* value = child ? utl::getValue(child,"thickness") : nullptr;
  if (value)
  {
    Thck0 = Thick = atof(value);
    IFEM::cout <<"\tConstant thickness: "<< Thck0 << std::endl;
  }

  child = elem->FirstChildElement("thickloss");
  if (child && child->FirstChild())
  {
    std::string ctrInside;
    utl::getAttribute(child,"t1",trOutside);
    utl::getAttribute(child,"t2",ctrInside);
    std::istringstream(child->FirstChild()->Value()) >> Xlow >> Xupp;
    if (ctrInside.find("t") != std::string::npos)
    {
      IFEM::cout <<"\tThickness loss function: ";
      thickLoss = utl::parseTimeFunc(ctrInside.c_str(),"expression");
    }
    else if ((trInside = atoi(ctrInside.c_str())) < 0.0)
      trInside = 0.0;
    if (trOutside < 0.0)
      trOutside = 0.0;
    IFEM::cout <<"\tThickness loss: t1="<< trOutside <<" t2="<< ctrInside
               <<"  Xlower = "<< Xlow <<"  Xupper = "<< Xupp << std::endl;
  }

  return true;
}


void AndesShell::setIntegrationPrm (unsigned short int i, double prm)
{
  this->ElasticBase::setIntegrationPrm(i,prm);
  if (beamProblem)
    beamProblem->setIntegrationPrm(i,prm);
}


void AndesShell::setMode (SIM::SolutionMode mode)
{
  if (isModal && mode == SIM::DYNAMIC)
    mode = SIM::RHS_ONLY;

  this->ElasticBase::setMode(mode);

  if (mode == SIM::STATIC || (isModal && mode == SIM::RHS_ONLY))
    iS = 0;
  else if (mode == SIM::DYNAMIC && eS == 1 && intPrm[4] == 0.5)
  {
    eS = 3; // Store external forces separately, for visualization
    vecNames.push_back("(empty)");
    vecNames.push_back("external forces");
  }

  if (beamProblem)
  {
    beamProblem->setMode(mode);
    beamProblem->setGravity(gravity.x,gravity.y,gravity.z);
  }
}


void AndesShell::initLHSbuffers (size_t newLHS)
{
  if (!newLHS)
  {
    if (eKm > 0 && !myKmats.empty() && !myKmats.front().empty())
      eKm = -eKm;
    if (eM  > 0 && !myMmats.empty() && !myMmats.front().empty())
      eM  = -eM;
  }
}


void AndesShell::initMatrixBuffers (size_t nEl, size_t)
{
  myKmats.resize(nEl);
  myMmats.resize(nEl);
}


void AndesShell::initForPatch (const ASMbase* pch)
{
  currentPatch = dynamic_cast<const ASMu2DNastran*>(pch);
  if (beamProblem)
    beamPatch = dynamic_cast<const ASMuBeam*>(pch);
}


bool AndesShell::setPressure (RealFunc* pf, int code,
                              const std::string& sName, const ASMbase* pch)
{
  size_t sIdx = 0;
  if (sName.empty())
    presFld[-code] = pf;
  else if ((sIdx = pch->getElementSetIdx(sName)))
    presFld[sIdx] = pf;
  else
    return false;

  return true;
}


RealFunc* AndesShell::getPressure (const ASMbase* pch) const
{
  for (const std::pair<const int,RealFunc*>& pf : presFld)
    if (pf.first < 0 || !(pch && pch->getElementSet(pf.first).empty()))
      return pf.second;

  return nullptr;
}


void AndesShell::initIntegration (size_t nGp, size_t)
{
  presVal.clear();
  if (this->havePressure())
    presVal.resize(nGp,{Vec3(),Vec3()});
}


LocalIntegral* AndesShell::getLocalIntegral (size_t nen, size_t iEl, bool) const
{
  if (beamPatch)
    return beamProblem->getLocalIntegral(nen,iEl,false);

  ElmMats* result = nullptr;
  if (this->inActive(iEl))
    return result; // element is not in current material group

  if (isModal)
    result = new ElmMats();
  else if (this->getMode(true) == SIM::DYNAMIC)
    result = new NewmarkMats(intPrm[0],intPrm[1],intPrm[2],intPrm[3]);
  else
    result = new ElmMats();

  switch (m_mode)
  {
    case SIM::STATIC:
    case SIM::MASS_ONLY:
      result->resize(1,1,nsd);
      break;

    case SIM::DYNAMIC:
      result->resize(3,eS,nsd);
      break;

    case SIM::VIBRATION:
      result->resize(2,0);
      break;

    case SIM::STIFF_ONLY:
      result->resize(1,0);
      break;

    case SIM::RHS_ONLY:
      result->resize(nSV > nCS ? 3 : 0, nSV > nCS ? eS : 1, nsd);
      result->rhsOnly = true;
      result->withLHS = false;
      break;

    default:
      std::cerr <<" *** AndesShell::getLocalIntegral: Mode flag "<< m_mode
                <<" is not yet supported."<< std::endl;
      delete result;
      return nullptr;
  }

  result->redim(npv*nen);
  return result;
}


ElmMats* AndesShell::getDofMatrices () const
{
  ElmMats* result = nullptr;

  if (isModal)
    return result;
  else if (this->getMode(true) == SIM::DYNAMIC)
    result = new NewmarkMats(0.0,0.0,intPrm[2],intPrm[3]);
  else
    result = new ElmMats();

  switch (m_mode)
  {
    case SIM::STATIC:
      result->resize(1,1);
      break;

    case SIM::DYNAMIC:
      result->resize(3,1);
      break;

    case SIM::VIBRATION:
    case SIM::STIFF_ONLY:
      result->resize(1,0);
      break;

    case SIM::RHS_ONLY:
      result->resize(nSV > nCS ? 3 : 0, 1);
      result->rhsOnly = true;
      result->withLHS = false;
      break;

    default:
      delete result;
      return nullptr;
  }

  result->redim(1);
  if (result->A.size() > 1)
    result->A[1].clear(); // No mass matrix
  return result;
}


void AndesShell::setSecondaryInt (GlobalIntegral* gq)
{
  delete myReacI;
  myReacI = gq;
}


GlobalIntegral& AndesShell::getGlobalInt (GlobalIntegral* gq) const
{
  if (m_mode == SIM::RHS_ONLY && myReacI)
    return *myReacI;

  return this->ElasticBase::getGlobalInt(gq);
}


int AndesShell::getIntegrandType () const
{
  if (beamPatch)
    return beamProblem->getIntegrandType();

  // Basis function derivatives (Jacobian, etc.) are needed only if surface
  // pressures (and/or gravity forces) exist, since the element matrices
  // are evaluated externally (from the finalizeElement() method).
  int itgType = eS > 0 && this->havePressure() ? STANDARD : NO_DERIVATIVES;
  if (this->getMode(true) == SIM::DYNAMIC) itgType |= POINT_DEFORMATION;
  if (thickLoss || trInside+trOutside > 0.0) itgType |= ELEMENT_CENTER;
  return itgType;
}


/*!
  In the case that element stiffness- and mass-matrices are cached,
  this method will copy the element matrices for current element from
  the cache and skip the property initialisation, which then is not needed
  (unless the element has gravity loads which require the shell thickness).
*/

bool AndesShell::initElement (const std::vector<int>& MNPC,
                              const FiniteElement& fe, const Vec3& Xc, size_t,
                              LocalIntegral& elmInt)
{
  if (beamPatch)
  {
    if (!beamPatch->getProps(fe.iel,Emod,GorNu,Rho,myBeamProps))
      return false;

    beamProblem->setStiffness(Emod,GorNu);
    beamProblem->setMass(Rho);
    return beamProblem->initElement(MNPC,fe,Vec3(),0,elmInt);
  }

  if (fe.idx < myKmats.size() && eKm < 0)
    static_cast<ElmMats&>(elmInt).A[-eKm-1] = myKmats[fe.idx];
  if (fe.idx < myMmats.size() && eM  < 0)
    static_cast<ElmMats&>(elmInt).A[-eM-1]  = myMmats[fe.idx];

  if (!this->initElement(MNPC,elmInt))
    return false;
  else if (fe.Xn.cols() < 3)
    return true; // No properties for 1-noded mass elements
  else if (eKm+eM <= 0 && !this->havePressure())
    return true; // Nothing to integrate for this shell element

  bool ok = true;
  if (!currentPatch)
    Thick = Thck0;
  else if (ovrMat) // Override the patch-level material properties
  {
    double dum1, dum2, dum3;
    ok = currentPatch->getProps(fe.iel,dum1,dum2,dum3,Thick);
  }
  else // Use the patch-level material properties
    ok = currentPatch->getProps(fe.iel,Emod,GorNu,Rho,Thick);

  // Scale the thickness depending on location inside or outside given box
  const Vec4* Xt = dynamic_cast<const Vec4*>(&Xc);
  if (Xt && thickLoss)
    trInside = std::min(std::max((*thickLoss)(Xt->t),0.0),1.0);
  if (trInside+trOutside > 0.0)
    Thick *= 1.0 - (Xc.inside(Xlow,Xupp) ? trInside : trOutside);

  return ok;
}


bool AndesShell::isInLossArea (const Vec3& Xc) const
{
  return thickLoss || trInside+trOutside > 0.0 ? Xc.inside(Xlow,Xupp) : false;
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
  if (beamPatch) // 2-noded beam element
    return beamProblem->evalInt(elmInt,fe,X);

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
  p *= fe.detJxW;
  Vector& Svec = static_cast<ElmMats&>(elmInt).b[eS-1];
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (unsigned short int i = 1; i <= 3; i++)
      Svec(npv*(a-1)+i) += fe.N(a)*p[i-1];

  // Integrate the total external load
  RealArray& sumLoad = static_cast<ElmMats&>(elmInt).c;
  for (size_t i = 0; i < sumLoad.size() && i < 3; i++)
    sumLoad[i] += p[i];

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
  if (beamPatch) // 2-noded beam element
    return beamProblem->finalizeElement(elmInt,fe,time);

  int iERR = -99;
  Vector vDummy(1);
  Matrix mDummy(1,1);
  Vector& Svec = eS  > 0 ? static_cast<ElmMats&>(elmInt).b[eS-1]  : vDummy;
  Matrix& Kmat = eKm > 0 ? static_cast<ElmMats&>(elmInt).A[eKm-1] : mDummy;
  Matrix& Mmat = eM  > 0 ? static_cast<ElmMats&>(elmInt).A[eM-1 ] : mDummy;
  RealArray& F = static_cast<ElmMats&>(elmInt).c;
  size_t nenod = fe.Xn.cols();
  if (currentPatch && nenod == 1) // 1-noded concentrated mass element
  {
    iERR = 0;
    if (eKm > 0)
      Kmat.clear(); // no stiffness
    if (eM > 0)
      currentPatch->getMassMatrix(fe.iel, Mmat);
    if (eS > 0)
    {
      currentPatch->getLoadVector(fe.iel, gravity, Svec);
      for (size_t i = 0; i < F.size() && i < Svec.size(); i++)
        F[i] += Svec[i];
    }
  }
  else if (currentPatch && nenod == 2) // 2-noded mass-less spring element
  {
    iERR = 0;
    if (eKm > 0)
      currentPatch->getStiffnessMatrix(fe.iel, Kmat);
    if (eM > 0)
      Mmat.clear(); // no mass
    if (eS > 0)
      Svec.clear(); // no external load
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
                 eKm > 0 ? Thick : 0.0, Emod, GorNu, eM > 0 ? Rho : 0.0,
                 Press.ptr(), Kmat.ptr(), Mmat.ptr(), Svec.ptr(), iERR);
#endif
    if (eS > 0)
      for (int a = 0; a < 3; a++)
        for (size_t i = 0; i < F.size() && npv*a+i < Svec.size(); i++)
          F[i] += Svec[npv*a+i];
  }
  else if (nenod == 4) // 4-noded shell element
  {
#ifdef HAS_ANDES
    // Invoke Fortran wrapper for the 4-noded ANDES element
    ifem_andes4_(fe.iel, fe.Xn.ptr(),
                 eKm > 0 ? Thick : 0.0, Emod, GorNu, eM > 0 ? Rho : 0.0,
                 Kmat.ptr(), Mmat.ptr(), iERR);
#endif
  }
  else
  {
    iERR = -98;
    std::cerr <<" *** AndesShell: Invalid element, nenod="<< nenod << std::endl;
  }

  const char* groupName = nullptr;
  if (iERR >= 1 && iERR <= 3)
  {
    groupName = "Degenerated elements";
    degenerated.insert(fe.iel);
  }
  else if (iERR == 4)
  {
    groupName = "Straight lines";
    straightline.insert(fe.iel);
  }
  else if (iERR == 9)
  {
    groupName = "Failed elements";
    failedElements.insert(fe.iel);
  }
  else if (iERR == -99)
    std::cerr <<" *** AndesShell: Built without this element."<< std::endl;

  if (currentPatch && groupName)
    const_cast<ASMu2DNastran*>(currentPatch)->addToElemSet(groupName,fe.iel);

  if (iERR >= 0)
  {
    if (fe.idx < myKmats.size() && Kmat.size() > 1) myKmats[fe.idx] = Kmat;
    if (fe.idx < myMmats.size() && Mmat.size() > 1) myMmats[fe.idx] = Mmat;
  }

  if (iS > 0 && !elmInt.vec.empty() && nenod > 1 && iERR >= 0)
  {
    Matrix& Sm = fe.idx < myKmats.size() ? myKmats[fe.idx] : Kmat;
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
  if (nenod <= 2) // 1-noded concentrated mass or 2-noded beam element
    s.clear();
  else if (nenod == 3) // 3-noded shell element
    s.clear();
  else if (nenod == 4) // 4-noded shell element
  {
    s.resize(n2v > 12 ? n2v : 12, 0.0);
#ifdef HAS_ANDES
    // Invoke Fortran wrapper for the 4-noded ANDES element
    ifem_strs24_(fe.iel, fe.Xn.ptr(), Thick, Emod, GorNu,
                 eV.front().ptr(), s.ptr(), s.ptr()+6, iERR);
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
  else if (s.size() >= 12)
  {
    // Calculate von Mises stresses only
    Vector vms;
    vms.reserve(n2v);
    for (size_t i = 6; i < 12; i += 3)
      vms.push_back(SymmTensor({s[i],s[i+1],s[i+2]}).vonMises());
    s.swap(vms);
  }

  return iERR == 0;
}


/*!
  If \a iEl is negative, this method returns \e true if gravity loads exists,
  or any element have surface pressure loads. Otherwise, it returns \e true
  only if the element \a iEl has surface pressure loads.
*/

bool AndesShell::havePressure (int iEl) const
{
  if (iEl < 0 && Rho > 0.0 && !gravity.isZero())
    return true;

  for (const std::pair<const int,RealFunc*>& press : presFld)
    if (press.first < 0 || iEl < 0)
      return true;
    else if (currentPatch && currentPatch->isInElementSet(press.first,-iEl))
      return true;

  return currentPatch ? currentPatch->haveLoads() : false;
}


void AndesShell::addPressure (Vec3& p, const Vec3& X,
                              const Vec3& n, int iEl) const
{
  for (const std::pair<const int,RealFunc*>& press : presFld)
    if (press.first < 0 || currentPatch->isInElementSet(press.first,-iEl))
      p += (*press.second)(X)*n;
}


void AndesShell::setParam (const std::string& name, const Vec3& value)
{
  for (const std::pair<const int,RealFunc*>& press : presFld)
    press.second->setParam(name,value);
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


void AndesShell::primaryScalarFields (Matrix& field)
{
  if (field.rows() != npv) return;

  // Lambda function returning the magnitude of a displacement vector.
  std::function<double(const double*)> absDis = [](const double* u) -> double
  {
    return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
  };

  // Insert the absolute value as the 7th solution component
  field.expandRows(1);
  for (size_t c = 1; c <= field.cols(); c++)
    field(npv+1,c) = absDis(field.ptr(c-1));
}


std::string AndesShell::getField2Name (size_t i, const char* prefix) const
{
  static const char* s[12] = { "n_xx", "n_yy", "n_xy", "m_xx", "m_yy", "m_xy",
                               "sigma_x", "sigma_y", "tau_xy",
                               "sigma_1", "sigma_2", "sigma_m" };

  if (n2v == 2) // output von Mises stresses only
    i = 6*i + 11;

  std::string name(s[i < 6 ? i : 6 + i%6]);
  if (i >= 12)
    name = "Top " + name;
  else if (i >= 6)
    name = "Bottom " + name;

  if (!prefix)
    return name;

  return prefix + std::string(" ") + name;
}
