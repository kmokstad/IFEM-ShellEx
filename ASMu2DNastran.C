// $Id$
//==============================================================================
//!
//! \file ASMu2DMastran.C
//!
//! \date Feb 27 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of unstructured 2D Lagrange FE models from Nastran BDF-file.
//!
//==============================================================================

#include "ASMu2DNastran.h"
#include "MPC.h"
#ifdef HAS_FFLLIB
#include "FFlLinkHandler.H"
#include "FFlNastranReader.H"
#include "FFlElementBase.H"
#include "FFlFEParts/FFlNode.H"
#include "FFlFEParts/FFlPMAT.H"
#include "FFlFEParts/FFlPMASS.H"
#include "FFlFEParts/FFlPTHICK.H"
#include "FFlFEParts/FFlPWAVGM.H"


class MyNastranReader : public FFlNastranReader
{
public:
  MyNastranReader(FFlLinkHandler& fePart, int lCount)
    : FFlNastranReader(&fePart,lCount) {}
  virtual ~MyNastranReader() {};

  bool readAndResolve(std::istream& is)
  {
    if (!this->resolve(this->read(is)))
      myLink->deleteGeometry(); // Parsing failure, delete all FE data
    else if (nWarnings+nNotes > 0)
      std::cout <<"\n  ** Parsing FE data succeeded."
                <<"\n     However, "<< nWarnings
                <<" warning(s) and "<< nNotes <<" note(s) were reported.\n"
                <<"     Review the messages and check the FE data file.\n";
    return myLink->hasGeometry();
  }
};
#endif


bool ASMu2DNastran::read (std::istream& is)
{
  myMLGE.clear();
  myMLGN.clear();
  myCoord.clear();

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
  FFlLinkHandler  fem;
  MyNastranReader reader(fem,lCount);
  if (!reader.readAndResolve(is) || !fem.resolve())
  {
    std::cout <<"\n *** Parsing/resolving FE data failed.\n"
              <<"     The FE model is probably not consistent and has not been"
              <<" resolved completely.\n";
    return false;
  }

  nnod = fem.getNodeCount();
  nel  = fem.getElementCount(FFlTypeInfoSpec::SHELL_ELM);
#ifdef INT_DEBUG
  fem.dump();
#else
  std::cout <<"\nTotal number of nodes:          "<< nnod
            <<"\nNumber of shell elements:       "<< nel
            <<"\nNumber of constraint elements:  "
            << fem.getElementCount(FFlTypeInfoSpec::CONSTRAINT_ELM)
            <<"\nNumber of other elements:       "
            << fem.getElementCount(FFlTypeInfoSpec::OTHER_ELM)
            << std::endl;
#endif

  myMLGN.reserve(nnod);
  myMLGE.reserve(nel);
  myCoord.reserve(nnod);

  for (NodesCIter n = fem.nodesBegin(); n != fem.nodesEnd(); ++n)
  {
    int nid = (*n)->getID();
    const FaVec3& X = (*n)->getPos();
#if INT_DEBUG > 1
    std::cout << myMLGN.size() <<" "<< nid <<": "<< X << std::endl;
#endif
    myMLGN.push_back(nid);
    myCoord.push_back(Vec3(X.x(),X.y(),X.z()));
    if ((*n)->isExternal())
      this->getNodeSet("ASET",lCount).push_back(myMLGN.size());
    else if ((*n)->isFixed())
      this->getNodeSet("SPC",lCount).push_back(myMLGN.size());
  }

  for (ElementsCIter e = fem.elementsBegin(); e != fem.elementsEnd(); ++e)
  {
    int eid = (*e)->getID();
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

#if INT_DEBUG > 1
      std::cout <<"Shell element "<< myMLGE.size() <<" "<< eid <<":";
      for (int node : mnpc) std::cout <<" "<< MLGN[node];
#endif
      myMLGE.push_back(eid);
      myMNPC.push_back(mnpc);

      ShellProps& sprop = myProps[eid];
      FFlPMAT* mat = dynamic_cast<FFlPMAT*>((*e)->getAttribute("PMAT"));
      if (mat)
      {
        sprop.Emod = mat->youngsModule.getValue();
        sprop.Rny  = mat->poissonsRatio.getValue();
        sprop.Rho  = mat->materialDensity.getValue();
      }
      else
        std::cout <<"  ** No material attached to element "<< eid
                  <<", using default properties"<< std::endl;

      FFlPTHICK* thk = dynamic_cast<FFlPTHICK*>((*e)->getAttribute("PTHICK"));
      if (thk)
        sprop.Thick = thk->thickness.getValue();
      else
        std::cout <<"  ** No shell thickness attached to element "<< eid
                  <<", using default value "<< sprop.Thick << std::endl;

#if INT_DEBUG > 1
      std::cout << sprop << std::endl;
#endif
    }
    else if ((*e)->getTypeName() == "RGD" && mnpc.size() > 1)
    {
#if INT_DEBUG > 1
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
#if INT_DEBUG > 1
      std::cout <<"Weighted average motion element "<< eid
                <<": reference node = "<< MLGN[mnpc.front()]
                <<"\n\tmasters =";
      for (size_t i = 1; i < mnpc.size(); i++) std::cout <<" "<< MLGN[mnpc[i]];
      std::cout << std::endl;
#endif
      int indC[6] = { -1, -1, -1, 0, 0, 0 };
      std::set<int> refC;
      std::vector<double> weights;
      FFlPWAVGM* wavgm = dynamic_cast<FFlPWAVGM*>((*e)->getAttribute("PWAVGM"));
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
      FFlPMASS* mass = dynamic_cast<FFlPMASS*>((*e)->getAttribute("PMASS"));
      if (mass)
      {
        const std::vector<double>& Mvec = mass->M.getValue();
        std::vector<double>::const_iterator m = Mvec.begin();

#if INT_DEBUG > 1
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
        std::cout <<"  ** No mass property attached to mass element "<< eid
                  <<" (ignored)"<< std::endl;
    }
    else
      std::cout <<"  ** Ignored element "<< (*e)->getTypeName()
                <<" "<< eid << std::endl;
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


#ifdef HAS_ANDES
extern "C" void wavgmconstreqn_(const int& iel, const int& lDof,
                                const int& nM, const int& nW, const int* indC,
                                const double* tenc, const double* weight,
                                const double& epsX, double* dX,
                                double* work, double* omega,
                                const int& ipsw, const int& lpu);
#endif


void ASMu2DNastran::addFlexibleCoupling (int iel, int lDof, const int* indC,
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
  const double Zero = 1.0e-6;
  double* rwork = new double[10*nM+3];
  double* omega = rwork+4*nM+3;

  //TODO: Convert this Fortran subroutine to C++
  wavgmconstreqn_(iel,lDof,nM,nW,indC,Xnod.ptr(),weights.data(),epsX,
                  rwork+nM,rwork,omega,ips,lpu);

  MPC* cons = new MPC(MLGN[mnpc.front()],lDof);
  if (this->addMPC(cons) && cons)
    for (int iM = 1; iM <= nM; iM++)
      for (int mDof = 1; mDof <= 6; mDof++, omega++)
        if (*omega < -Zero || *omega > Zero)
          cons->addMaster(MLGN[mnpc[iM]],mDof,*omega);

  delete[] rwork;
#endif
}
