// $Id$
//==============================================================================
//!
//! \file ASMu2DNastran.C
//!
//! \date Feb 27 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of 2D %Lagrange FE models from Nastran Bulk Data File (BDF).
//!
//==============================================================================

#include "ASMu2DNastran.h"
#include "FiniteElement.h"
#include "IntegrandBase.h"
#include "MPC.h"
#include "Vec3Oper.h"
#ifdef HAS_FFLLIB
#include "IFEM.h"
#include "FFlLinkHandler.H"
#include "FFlNastranReader.H"
#include "FFlElementBase.H"
#include "FFlLoadBase.H"
#include "FFlFEParts/FFlNode.H"
#include "FFlFEParts/FFlPMAT.H"
#include "FFlFEParts/FFlPMASS.H"
#include "FFlFEParts/FFlPTHICK.H"
#include "FFlFEParts/FFlPWAVGM.H"


/*!
  \brief Class for reading Nastran bulk data into a FFlLinkHandler object.
*/

class MyNastranReader : public FFlNastranReader
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  MyNastranReader(FFlLinkHandler& fePart, int lCount)
    : FFlNastranReader(&fePart,lCount) {}
  //! \brief Empty destructor.
  virtual ~MyNastranReader() {}

  //! \brief Reads the FE model and resolves all topological references.
  bool readAndResolve(std::istream& is)
  {
    if (!this->resolve(this->read(is)))
      myLink->deleteGeometry(); // Parsing failure, delete all FE data
    else if (nWarnings+nNotes > 0)
      IFEM::cout <<"\n  ** Parsing FE data succeeded."
                 <<"\n     However, "<< nWarnings
                 <<" warning(s) and "<< nNotes <<" note(s) were reported.\n"
                 <<"     Review the messages and check the FE data file.\n"
                 << std::endl;
    return myLink->hasGeometry() && myLink->resolve();
  }
};
#endif


bool ASMu2DNastran::read (std::istream& is)
{
  myMLGE.clear();
  myMLGN.clear();
  myCoord.clear();
  myLoads.clear();

  // Fast-forward until "BEGIN BULK"
  int lCount = 0;
  char cline[256];
  while (is.getline(cline,255))
    if (strncmp(cline,"BEGIN BULK",10))
      ++lCount;
    else
      break;

  if (!is) return false; // No bulk data file

#ifdef HAS_FFLLIB
#define GET_ATTRIBUTE(el,att) dynamic_cast<FFl##att*>((*el)->getAttribute(#att))

  FFlLinkHandler  fem;
  MyNastranReader reader(fem,lCount);
  if (!reader.readAndResolve(is))
  {
    std::cerr <<"\n *** Parsing/resolving FE data failed.\n"
              <<"     The FE model is probably not consistent and has not been"
              <<" resolved completely.\n";
    return false;
  }

  nnod = fem.getNodeCount(FFlLinkHandler::FFL_FEM);
  nel  = fem.getElementCount(FFlTypeInfoSpec::SHELL_ELM);
#ifdef INT_DEBUG
  fem.dump();
#else
  IFEM::cout <<"\nTotal number of nodes:          "<< nnod
             <<"\nNumber of shell elements:       "<< nel
             <<"\nNumber of constraint elements:  "
             << fem.getElementCount(FFlTypeInfoSpec::CONSTRAINT_ELM)
             <<"\nNumber of other elements:       "
             << fem.getElementCount(FFlTypeInfoSpec::OTHER_ELM)
             << std::endl;
#endif
  size_t allNodes = fem.getNodeCount(FFlLinkHandler::FFL_ALL);
  if (allNodes > nnod)
    IFEM::cout <<"\n  ** Warning: This model contains "<< allNodes-nnod
               <<" node(s) without any element connections (ignored)."
               <<"\n     Please check the FE data file.\n"<< std::endl;

  myMLGN.reserve(nnod);
  myMLGE.reserve(nel);
  myCoord.reserve(nnod);

  // Extract the nodal points
  for (size_t inod = 1; inod <= nnod; inod++)
  {
    FFlNode* node = fem.getFENode(inod);
    int nid = node->getID();
    const FaVec3& X = node->getPos();
#if INT_DEBUG > 10
    std::cout << myMLGN.size() <<" "<< nid <<": "<< X << std::endl;
#endif
    myMLGN.push_back(nid);
    myCoord.push_back(Vec3(X.x(),X.y(),X.z()));

    if (node->isExternal()) // Create a node set for the supernodes
      this->addToNodeSet("ASET",myMLGN.size());
    else if (node->isFixed())
    {
      // Create a node set for prescribed nodes,
      // one for each DOF constellation
      std::string cstat = std::to_string(-node->getStatus(-64));
      this->addToNodeSet("SPC"+cstat,myMLGN.size());
    }
  }

  // Extract the element data
  for (ElementsCIter e = fem.elementsBegin(); e != fem.elementsEnd(); ++e)
  {
    int eid = (*e)->getID();

    // Nodal connectivities
    IntVec mnpc((*e)->getNodeCount());
    for (size_t i = 0; i < mnpc.size(); i++)
      mnpc[i] = this->getNodeIndex((*e)->getNodeID(1+i)) - 1;

    if ((*e)->getCathegory() == FFlTypeInfoSpec::SHELL_ELM)
    {
      if (mnpc.size() == 4)
      {
        // Need to swap nodes 3 and 4 into IFEM order
        std::swap(mnpc[2],mnpc[3]);
        swapNode34 = true;
      }

#if INT_DEBUG > 10
      std::cout <<"Shell element "<< myMLGE.size() <<" "<< eid <<":";
      for (int node : mnpc) std::cout <<" "<< MLGN[node];
#endif
      myMLGE.push_back(eid);
      myMNPC.push_back(mnpc);

      // Extract material properties
      ShellProps& sprop = myProps[eid];
      FFlPMAT* mat = GET_ATTRIBUTE(e,PMAT);
      if (mat)
      {
        sprop.Emod = mat->youngsModule.getValue();
        sprop.Rny  = mat->poissonsRatio.getValue();
        sprop.Rho  = mat->materialDensity.getValue();
      }
      else
        IFEM::cout <<"  ** No material attached to element "<< eid
                   <<", using default properties"<< std::endl;

      // Extract shell thickness
      FFlPTHICK* thk = GET_ATTRIBUTE(e,PTHICK);
      if (thk)
        sprop.Thick = thk->thickness.getValue();
      else
        IFEM::cout <<"  ** No shell thickness attached to element "<< eid
                   <<", using default value "<< sprop.Thick << std::endl;

#if INT_DEBUG > 10
      std::cout <<" t="<< sprop.Thick <<" E="<< sprop.Emod <<" nu="<< sprop.Rny
                <<" rho="<< sprop.Rho << std::endl;
#endif
    }
    else if ((*e)->getTypeName() == "RGD" && mnpc.size() > 1)
    {
#if INT_DEBUG > 5
      std::cout <<"Rigid element "<< eid <<": master = "<< MLGN[mnpc.front()]
                <<" slaves =";
      for (size_t i = 1; i < mnpc.size(); i++) std::cout <<" "<< MLGN[mnpc[i]];
      std::cout << std::endl;
#endif
      for (int& n : mnpc) ++n; // Need 1-based node indices
      this->addRigidCouplings((*e)->getNodeID(1),myCoord[mnpc.front()-1],
                              IntVec(mnpc.begin()+1,mnpc.end()));
    }
    else if ((*e)->getTypeName() == "WAVGM" && mnpc.size() > 1)
    {
#if INT_DEBUG > 5
      std::cout <<"Weighted average motion element "<< eid
                <<": reference node = "<< MLGN[mnpc.front()]
                <<"\n\tmasters =";
      for (size_t i = 1; i < mnpc.size(); i++) std::cout <<" "<< MLGN[mnpc[i]];
      std::cout << std::endl;
#endif
      int indC[6] = { -1, -1, -1, 0, 0, 0 };
      std::set<int> refC;
      std::vector<double> weights;
      FFlPWAVGM* wavgm = GET_ATTRIBUTE(e,PWAVGM);
      if (wavgm)
      {
        int div = 100000;
        int dofIds = wavgm->refC.getValue();
        if (dofIds > 0)
          for (int i = 0; i < 6; i++)
          {
            int dof = dofIds/div;
            if (dof > 0 && dof < 7) refC.insert(dof);
            dofIds -= dof*div;
            div    /= 10;
          }
        for (int i = 0; i < 6; i++)
          indC[i] = wavgm->indC[i].getValue();
        weights = wavgm->weightMatrix.getValue();
      }
      else
        for (int lDof = 1; lDof <= 6; lDof++)
          refC.insert(lDof);

      size_t icol = 0;
      Matrix Xnod(nsd,mnpc.size());
      for (int inod : mnpc)
        Xnod.fillColumn(++icol,this->getCoord(1+inod).ptr());
      for (int lDof : refC)
        this->addFlexibleCoupling(eid,lDof,indC,weights,mnpc,Xnod);
    }
    else if ((*e)->getTypeName() == "CMASS" && !mnpc.empty())
    {
      FFlPMASS* mass = GET_ATTRIBUTE(e,PMASS);
      if (mass)
      {
        const std::vector<double>& Mvec = mass->M.getValue();
        std::vector<double>::const_iterator m = Mvec.begin();

#if INT_DEBUG > 5
        std::cout <<"Mass element "<< myMLGE.size() <<" "<< eid
                  <<": node = "<< MLGN[mnpc.front()] << std::endl;
#endif
        myMLGE.push_back(eid);
        myMNPC.push_back({mnpc.front()});
        Matrix& M = myMass[eid];
        M.resize(6,6);
        for (int i = 1; i <= 6; i++)
          for (int j = 1; j <= i && m != Mvec.end(); j++)
          {
            M(i,j) = *(m++);
            if (j < i) M(j,i) = M(i,j);
          }
      }
      else
        IFEM::cout <<"  ** No mass property attached to mass element "<< eid
                   <<" (ignored)"<< std::endl;
    }
    else
      IFEM::cout <<"  ** Ignored element "<< (*e)->getTypeName()
                 <<" "<< eid << std::endl;
  }

  // Extract the pressure loads, if any
  std::set<int> loadCases;
  fem.getLoadCases(loadCases);
  for (int lc : loadCases)
  {
    LoadsVec loads;
    fem.getLoads(lc,loads);
    for (FFlLoadBase* load : loads)
    {
      std::vector<FaVec3> Pe;
      int iel, iface = 0;
      while ((iel = load->getLoad(Pe,iface)) > 0)
        if (iface >= 0)
        {
          Vec3Vec& elLoad = myLoads[iel];
          elLoad.clear();
          elLoad.reserve(Pe.size());
          for (const FaVec3& P : Pe)
            elLoad.push_back(P.getPt());
          if (elLoad.size() == 4)
            std::swap(elLoad[2],elLoad[3]);
#if INT_DEBUG > 5
          std::cout <<"Surface pressure on element "<< iel <<":";
          for (const Vec3& p : elLoad) std::cout <<"  "<< p;
          std::cout << std::endl;
#endif
        }
    }
  }
#endif

  return true;
}


bool ASMu2DNastran::getProps (int eId, double& E, double& nu,
                              double& rho, double& t) const
{
  std::map<int,ShellProps>::const_iterator it = myProps.find(eId);
  if (it == myProps.end())
  {
    std::cerr <<" *** No properties for shell element "<< eId << std::endl;
    return false;
  }

  E   = it->second.Emod;
  nu  = it->second.Rny;
  rho = it->second.Rho;
  t   = it->second.Thick;

  return true;
}


bool ASMu2DNastran::getThickness (int eId, double& t) const
{
  if (myMass.find(eId) != myMass.end())
    t = 0.0; // Silently ignore mass elements
  else
  {
    std::map<int,ShellProps>::const_iterator it = myProps.find(eId);
    if (it != myProps.end())
      t = it->second.Thick;
    else
    {
      std::cerr <<" *** No properties for shell element "<< eId << std::endl;
      return false;
    }
  }

  return true;
}


bool ASMu2DNastran::getMassMatrix (int eId, Matrix& M) const
{
  std::map<int,Matrix>::const_iterator it = myMass.find(eId);
  if (it == myMass.end())
  {
    std::cerr <<" *** No mass matrix for element "<< eId << std::endl;
    return false;
  }

  M = it->second;

  return true;
}


bool ASMu2DNastran::getLoadVector (int eId, const Vec3& g, Vector& S) const
{
  std::map<int,Matrix>::const_iterator it = myMass.find(eId);
  if (it == myMass.end())
  {
    std::cerr <<" *** No mass matrix for element "<< eId << std::endl;
    return false;
  }

  for (size_t i = 1; i <= 3 && i < it->second.rows() && i <= S.size(); i++)
    S(i) = it->second(i,1)*g.x + it->second(i,2)*g.y + it->second(i,3)*g.z;

  return true;
}


bool ASMu2DNastran::addPressureAt (Vec3& p, int eId, const RealArray& N) const
{
  std::map<int,Vec3Vec>::const_iterator it = myLoads.find(eId);
  if (it == myLoads.end() || it->second.empty()) return false;

  size_t n = it->second.size();
  if (n == 1 || N.size() > n)
    p += it->second.front();
  else if (N.size() == 1 && N.front() < 0.0)
  {
    // N[0] contains the nodal index
    size_t idx = -N.front();
    if (idx > 0 && idx <= n)
      p += it->second[idx-1];
  }

  for (size_t i = 0; i < N.size(); i++)
    p += N[i]*it->second[i];

  return true;
}


#ifdef HAS_ANDES
extern "C" void wavgmconstreqn_(const int& iel, const int& lDof,
                                const int& nM, const int& nW, const int* indC,
                                const double* tenc, const double* weight,
                                const double& epsX, double* dX,
                                double* work, double* omega,
                                const int& ipsw, const int& lpu);
#endif


void ASMu2DNastran::addFlexibleCoupling (int eId, int lDof, const int* indC,
                                         const std::vector<double>& weights,
                                         const IntVec& mnpc, const Matrix& Xnod)
{
#ifdef HAS_ANDES
  const int nM = Xnod.cols() - 1;
  const int nW = weights.size();
#if INT_DEBUG > 2
  const int ips = INT_DEBUG;
#else
  const int ips = 0;
#endif
  const int lpu = 6;
  const double epsX = 1.0e-4;
  const double Zero = 1.0e-8;
  double* rwork = new double[10*nM+3];
  double* omega = rwork+4*nM+3;

  //TODO: Convert this Fortran subroutine to C++
  wavgmconstreqn_(eId,lDof,nM,nW,indC,Xnod.ptr(),weights.data(),epsX,
                  rwork+nM,rwork,omega,ips,lpu);

  MPC* cons = new MPC(MLGN[mnpc.front()],lDof);
  if (this->addMPC(cons) && cons)
    for (int iM = 1; iM <= nM; iM++)
      for (int mDof = 1; mDof <= 6; mDof++, omega++)
        if (*omega < -Zero || *omega > Zero)
          cons->addMaster(MLGN[mnpc[iM]],mDof,*omega);
#ifdef INT_DEBUG
        else
          std::cout <<"  ** Ignoring small coupling coefficient "<< *omega
                    <<" to local dof "<< mDof
                    <<" of master node "<< MLGN[mnpc[iM]]
                    <<" in RBE3 element "<< eId << std::endl;
#endif
  delete[] rwork;
#endif
}


/*!
  This method overrides the parent class method to always evaluate the secondary
  solution at the nodal points of the patch.
*/

bool ASMu2DNastran::evalSolution (Matrix& sField, const IntegrandBase& integr,
                                  const int*, char) const
{
  return this->evalSolution(sField,integr,nullptr,false);
}


/*!
  This method overrides the parent class method to always evaluate the secondary
  solution at the element nodes and then perform nodal averaging to obtain the
  unique nodal values, or perform direct evaluation at the element centers.
  It is assumed that all calculations are performed by the
  IntegrandBase::evalSol() call, therefore no basis function evaluations here.
*/

bool ASMu2DNastran::evalSolution (Matrix& sField, const IntegrandBase& integr,
                                  const RealArray*, bool atElmCenters) const
{
  sField.clear();

  FiniteElement fe;
  Vector        solPt;
  Vectors       globSolPt(atElmCenters ? 0 : nnod);
  IntVec        checkPt(atElmCenters ? 0 : nnod,0);

  // Evaluate the secondary solution field at each element node or center
  for (size_t iel = 1; iel <= nel; iel++)
    if ((fe.iel = MLGE[iel-1]) > 0) // ignore the zero-area elements
    {
      if (!this->getElementCoordinates(fe.Xn,iel))
        return false;

      const IntVec& mnpc = MNPC[iel-1];
      const size_t nenod = atElmCenters ? 1 : mnpc.size();
      for (size_t loc = 0; loc < nenod; loc++)
      {
        if (nenod == 3)
          switch (1+loc) {
          case 1: fe.xi = 1.0; fe.eta = 0.0; break;
          case 2: fe.xi = 0.0; fe.eta = 1.0; break;
          case 3: fe.xi = 0.0; fe.eta = 0.0; break;
          }
        else if (nenod == 4)
        {
          fe.xi  = -1.0 + 2.0*(loc%2);
          fe.eta = -1.0 + 2.0*(loc/2);
        }
        else if (mnpc.size() == 3)
          fe.xi = fe.eta = 1.0/3.0;
        else
          fe.xi = fe.eta = 0.0;

        if (!integr.evalSol(solPt,fe,Vec3(),mnpc))
          return false;
        else if (solPt.empty())
          break; // a valid element with no secondary solution

        if (sField.empty())
          sField.resize(solPt.size(), atElmCenters ? nel : nnod, true);

        if (atElmCenters)
          sField.fillColumn(iel,solPt);
        else if (++checkPt[mnpc[loc]] == 1)
          globSolPt[mnpc[loc]] = solPt;
        else
          globSolPt[mnpc[loc]] += solPt;
      }
    }

  // Nodal averaging
  for (size_t i = 0; i < checkPt.size(); i++)
    if (checkPt[i])
      sField.fillColumn(1+i, globSolPt[i] /= static_cast<double>(checkPt[i]));

  return true;
}
