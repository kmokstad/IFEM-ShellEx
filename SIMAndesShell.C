// $Id$
//==============================================================================
//!
//! \file SIMAndesShell.C
//!
//! \date Feb 20 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for FE analysis using the ANDES shell elements.
//!
//==============================================================================

#include "SIMAndesShell.h"
#include "ASMu2DNastran.h"
#include "AndesShell.h"
#include "ElasticBeam.h"
#include "AlgEqSystem.h"
#include "ElmMats.h"
#include "SystemMatrix.h"
#include "SAM.h"
#include "DataExporter.h"
#include "HDF5Writer.h"
#include "Functions.h"
#include "WaveFunc.h"
#include "Utilities.h"
#include "ElementBlock.h"
#include "Vec3Oper.h"
#include "VTF.h"
#include "IFEM.h"
#include "tinyxml2.h"
#include <fstream>

namespace ASM { extern char useBeam; extern IntVec fixRBE3; }
namespace Elastic { extern double time; }

using Parent = SIMElasticity<SIM2D>; //!< Convenience type alias


SIMAndesShell::SIMAndesShell (unsigned short int n, bool m) : nss(n), modal(m)
{
  nsd = 3;
  nf.front() = 6;
  seaLx = seaLy = seaGridSize = 0.0;
  seaBlock = 0;
  seasurf = nullptr;
  shellp = nullptr;
}


SIMAndesShell::~SIMAndesShell ()
{
  for (PointLoad& load : myLoads)
    delete load.p;
}


bool SIMAndesShell::printProblem () const
{
  for (const ASMbase* pch : myModel)
    if (dynamic_cast<const ASMu1DLag*>(pch))
    {
      // To enable correct print of beam element properties (for first element)
      if (ElasticBase* p = const_cast<SIMAndesShell*>(this)->getIntegrand(); p)
        p->initForPatch(pch);
      break;
    }

  return this->Parent::printProblem();
}


ElasticBase* SIMAndesShell::getIntegrand ()
{
  if (!myProblem)
    myProblem = shellp = new AndesShell(nss,modal,ASM::useBeam);

  return shellp;
}


bool SIMAndesShell::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),"geometry"))
  {
    const tinyxml2::XMLElement* child = elem->FirstChildElement("fixRBE3");
    if (child && child->FirstChild())
      utl::parseIntegers(ASM::fixRBE3,child->FirstChild()->Value());

    std::string fName;
    child = elem->FirstChildElement("patchfile");
    if (child && child->FirstChild())
      if (utl::getAttribute(child,"type",fName,true) && fName == "nastran")
      {
        // Extract the path of specified nastran file
        fName = child->FirstChild()->Value();
        size_t ipos = fName.find_last_of("/\\");
        if (ipos < std::string::npos)
          myPath = fName.substr(0,ipos);
        // Check if the 1-noded mass elements should be placed
        // in a separate group for multi-threaded assembly
        utl::getAttribute(child,"separatePointMasses",lagMTOK);
      }
  }

  if (!this->Parent::parse(elem))
    return false;

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(elem->Value(),"postprocessing"))
    {
      if (!strcasecmp(child->Value(),"seasurface"))
      {
        utl::getAttribute(child,"X0",seaX0);
        utl::getAttribute(child,"Lx",seaLx);
        utl::getAttribute(child,"Ly",seaLy);
        utl::getAttribute(child,"gridsize",seaGridSize);
        IFEM::cout <<"\tSea surface ["<< seaX0
                   <<"] to ["<< seaX0 + Vec3(seaLx,seaLy,0.0)
                   <<"], gridSize="<< seaGridSize;
        const tinyxml2::XMLElement* zeta = child->FirstChildElement("zeta");
        if (zeta && zeta->FirstChild() && !seasurf)
        {
          std::string type("expression");
          utl::getAttribute(zeta,"type",type);
          IFEM::cout <<" ("<< type <<")";
          seasurf = utl::parseRealFunc(zeta->FirstChild()->Value(),type);
        }
        IFEM::cout << std::endl;
      }
      else if (!strcasecmp(child->Value(),"reactions"))
      {
        if (!utl::getAttribute(child,"set",myRFset))
          myRFset = "(all)";
      }
      else if (ElasticBase* problem = this->getIntegrand(); problem)
        problem->parse(child);
    }
    else if (strcasecmp(elem->Value(),"elasticity"))
      continue; // The remaining should be within the elasticity context

    else if (!strcasecmp(child->Value(),"seasurface") && !seasurf)
    {
      IFEM::cout <<"  Parsing <seasurface>"<< std::endl;

      std::string type("expression");
      utl::getAttribute(child,"type",type,true);
      IFEM::cout <<"\tSea surface function ("<< type <<")";
      if (type == "spectrum")
      {
        double g = this->getIntegrand()->getGravity().length();
        seasurf = WaveSpectrum::parse(child,g);
      }
      else
        seasurf = utl::parseRealFunc(child->FirstChild()->Value(),type);
      IFEM::cout << std::endl;
    }

    else if (!strcasecmp(child->Value(),"pressure"))
    {
      IFEM::cout <<"  Parsing <pressure>"<< std::endl;

      std::string set, type;
      utl::getAttribute(child,"set",set);
      int code = this->getUniquePropertyCode(set,1);
      if (code == 0) utl::getAttribute(child,"code",code);
      if (code > 0)
      {
        utl::getAttribute(child,"type",type,true);
        IFEM::cout <<"\tPressure code "<< code;
        if (!type.empty()) IFEM::cout <<" ("<< type <<")";
        if (type == "hydrostatic" && seasurf)
        {
          double g = this->getIntegrand()->getGravity().length();
          double z0 = 0.0, rhow = 1000.0;
          utl::getAttribute(child,"waterline",z0);
          utl::getAttribute(child,"waterdensity",rhow);
          IFEM::cout <<": z0 = "<< z0 <<" rhow = "<< rhow;
          myScalars[code] = new HydroStaticPressure(*seasurf,z0,g,rhow);
        }
        else if (child->FirstChild())
          myScalars[code] = utl::parseRealFunc(child->FirstChild()->Value(),
                                               type);

        this->setPropertyType(code,Property::BODYLOAD);
        IFEM::cout << std::endl;
        AndesShell* shellInt = dynamic_cast<AndesShell*>(this->getIntegrand());
        if (!shellInt) continue; // Empty integrand - skip

        if (set.empty()) // Applies to all elements in the model
          shellInt->setPressure(myScalars[code],code);
        else for (const ASMbase* pch : myModel)
          if (shellInt->setPressure(myScalars[code],code,set,pch))
            break; // Note: This assumes a set has elements from one patch only
      }
    }
    else if (!strcasecmp(child->Value(),"spring") && child->FirstChild())
    {
      IFEM::cout <<"  Parsing <spring>"<< std::endl;

      int inod = 0, ldof = 0;
      utl::getAttribute(child,"node",inod);
      utl::getAttribute(child,"dof",ldof);
      double coeff = atof(child->FirstChild()->Value());

      std::string set;
      utl::getAttribute(child,"set",set);
      IntVec nodes = this->getNodeSet(set);
      if (ldof < 1 || ldof > 6)
        nodes.clear();
      else if (nodes.empty() && inod > 0)
        nodes = { inod };

      for (int n : nodes)
      {
        IFEM::cout <<"\tNode "<< n <<" dof "<< ldof
                   <<" Spring stiffness: "<< coeff << std::endl;
        mySprings.emplace_back(n,ldof,coeff);
      }
    }
    else if (!strcasecmp(child->Value(),"nodeload") && child->FirstChild())
    {
      IFEM::cout <<"  Parsing <nodeload>"<< std::endl;

      int inod = 0, ldof = 0;
      ScalarFunc* f = nullptr;
      utl::getAttribute(child,"node",inod);
      utl::getAttribute(child,"dof",ldof);

      if (inod > 0 && ldof > 0 && ldof <= 6)
      {
        std::string type("constant");
        utl::getAttribute(child,"type",type);

        IFEM::cout <<"\tNode "<< inod <<" dof "<< ldof <<" Load: ";
        if (type == "constant")
        {
          f = new ConstantFunc(atof(child->FirstChild()->Value()));
          IFEM::cout << (*f)(0.0) << std::endl;
        }
        else
          f = utl::parseTimeFunc(child->FirstChild()->Value(),type);

        myLoads.emplace_back(inod,ldof,f);
      }
    }

  return true;
}


ASMbase* SIMAndesShell::readPatch (std::istream& isp, int, const CharVec&,
                                   const char*) const
{
  ASMbase* pch = NULL;
  ASMu2DNastran* shell = NULL;
  if (nf.size() == 2 && nf[1] == 'n') // Nastran bulk data file
    pch = shell = new ASMu2DNastran(nsd,nf.front(),myPath);
  else if (!(pch = ASM2D::create(opt.discretization,nsd,nf)))
    return pch;

  if (!pch->read(isp))
  {
    delete pch;
    return nullptr;
  }

  if (shell)
  {
    // Check if we also have beam elements in the model.
    // They will be kept in a separate patch of 1D elements.
    ASMbase* bpch = shell->haveBeams();
    if (bpch)
    {
      bpch->idx = myModel.size();
      const_cast<SIMAndesShell*>(this)->myModel.push_back(bpch);
    }
    else
      ASM::useBeam = 0;
  }

  pch->idx = myModel.size();
  if (shell && ASM::useBeam && pch->getNoElms() > 0)
    IFEM::cout <<"\tCreated shell patch "<< pch->idx+1
               <<" with "<< pch->getNoElms() <<" elements"<< std::endl;
  else
    IFEM::cout <<"\tReading patch "<< pch->idx+1 << std::endl;

  if (shell)
  {
    // Add topology items for the predefined node- and element sets
    std::string name;
    SIMinput* sim = const_cast<SIMAndesShell*>(this);
    for (int iset = 1; shell->getNodeSet(iset,name); iset++)
      sim->topology(name).emplace(pch->idx+1,iset,4);
    for (int iset = 1; shell->getElementSet(iset,name); iset++)
      sim->topology(name).emplace(pch->idx+1,iset,5);
  }

  return pch;
}


bool SIMAndesShell::preprocessB ()
{
  for (ASMbase* pch : myModel)
    if (ASMu2DNastran* shl = dynamic_cast<ASMu2DNastran*>(pch); shl)
      shl->initPressureCache();

  return true;
}


bool SIMAndesShell::dumpShellMesh (const char* fname) const
{
  for (const ASMbase* pch : myModel)
    if (const ASMu2DLag* shell = dynamic_cast<const ASMu2DLag*>(pch); shell)
      return shell->writeXML(fname);

  return false; // No shell element patch
}


size_t SIMAndesShell::getNoShellElms () const
{
  size_t nShells = 0;
  for (const ASMbase* pch : myModel)
    if (dynamic_cast<const ASMu2DLag*>(pch))
      for (size_t iel = 1; iel <= pch->getNoElms(true); iel++)
        if (pch->getElementNodes(iel).size() > 2) nShells++;

  return nShells;
}


void SIMAndesShell::getShellThicknesses (RealArray& elmThick) const
{
  // Include also the collapsed and non-shell elements,
  // because that is what the VTF-writer expects
  elmThick.clear();
  elmThick.reserve(this->getNoElms(false,true));

  int missing = 0;
  for (const ASMbase* pch : myModel)
    if (const ASMu2DNastran* shl = dynamic_cast<const ASMu2DNastran*>(pch); shl)
    {
      for (size_t iel = 1; iel <= pch->getNoElms(true); iel++)
        if (pch->getElementNodes(iel).size() > 1) // skip 1-noded mass elements
        {
          elmThick.push_back(0.0);
          if (!shl->getThickness(iel-1,elmThick.back()))
            ++missing;
        }
    }
    else
      elmThick.insert(elmThick.end(),pch->getNoElms(true),0);

  if (missing > 1)
    IFEM::cout <<" *** A total of "<< missing <<" elements lack thickness.\n"
               <<"     Please check your model for inconsistency."<< std::endl;
}


bool SIMAndesShell::getElementGroup (int iset, std::string& name,
                                     RealArray& elGroup) const
{
  elGroup.clear();
  elGroup.reserve(this->getNoElms(false,true));

  for (const ASMbase* pch : myModel)
    if (const ASMu2DLag* shell = dynamic_cast<const ASMu2DLag*>(pch);
        shell && shell->getElementSet(iset,name))
    {
      for (size_t iel = 1; iel <= pch->getNoElms(true); iel++)
        if (pch->getElementNodes(iel).size() > 1) // skip 1-noded mass elements
          elGroup.push_back(pch->getElmID(iel) > 0 &&
                            shell->isInElementSet(iset,iel));
    }
    else
      elGroup.insert(elGroup.end(),pch->getNoElms(true),0);

  return !elGroup.empty();
}


bool SIMAndesShell::renumberNodes (const std::map<int,int>& nodeMap)
{
  bool ok = this->Parent::renumberNodes(nodeMap);

  for (DOFspring& spr : mySprings)
    if (spr.inod > 0)
      ok &= utl::renumber(spr.inod,nodeMap,true);

  for (PointLoad& load : myLoads)
    if (load.inod > 0)
      ok &= utl::renumber(load.inod,nodeMap,true);

  return ok;
}


/*!
  This method is overridden in case a beam patch exists.
*/

bool SIMAndesShell::extractPatchSolution (IntegrandBase* itg,
                                          const Vectors& sol,
                                          size_t pindx) const
{
  if (itg == myProblem && dynamic_cast<ASMuBeam*>(this->getPatch(pindx+1)))
    if (shellp) itg = shellp->hasBeamProblem();

  return this->Parent::extractPatchSolution(itg,sol,pindx);
}


/*!
  This method is invoked after the assembly loop over the elements is finished.
  It will therefore check if there were any failed elements during the assembly,
  and return \e false in that case such that the simulation will be aborted.
*/

bool SIMAndesShell::assembleDiscreteItems (const IntegrandBase* itg,
                                           const TimeDomain& time,
                                           const Vectors& sol)
{
  bool ok = !shellp || shellp->allElementsOK();
  if (itg != myProblem || !myEqSys)
    return ok;

  if (!mySprings.empty() && shellp)
    if (ElmMats* sprMat = shellp->getDofMatrices(); sprMat)
    {
      // Assemble discrete DOF spring stiffnesses and associated internal forces
      SystemMatrix* K = (sprMat->A.empty() ? nullptr : myEqSys->getMatrix());
      SystemVector* S = (sprMat->b.empty() || sol.front().zero(1e-15) ?
                         nullptr : myEqSys->getVector());
      sprMat->setStepSize(time.dt,time.it);
      for (const DOFspring& spr : mySprings)
      {
        if (K)
          sprMat->A.back().fill(spr.coeff);
        if (S)
        {
          double v = mySam->getDofVal(sol.front(),spr.inod,spr.ldof);
          sprMat->b.back().fill(-spr.coeff*v);
        }
#if INT_DEBUG > 1
        std::cout <<"\nSpring at DOF "<< spr.inod <<","<< spr.ldof
                  <<": "<< spr.coeff;
        if (S)
          std::cout <<" force = "<< -sprMat->b.back()(1);
#endif
        if (K)
          ok &= K->add(sprMat->getNewtonMatrix()(1,1),
                       mySam->getEquation(spr.inod,spr.ldof));
        if (S)
          ok &= mySam->assembleSystem(*S,sprMat->getRHSVector()(1),
                                      {spr.inod,spr.ldof});
      }
    }

  if (!myLoads.empty())
  {
    // Assemble external nodal point loads
    const size_t nrhs = myEqSys->getNoRHS();
    SystemVector* R = nrhs > 0 ? myEqSys->getVector(nrhs-1) : nullptr;
    for (const PointLoad& load : myLoads)
    {
      double P = (*load.p)(time.t);
      int ldof = load.ldof;
      if (R)
      {
        myEqSys->addScalar(P,ldof-1);
        ok &= mySam->assembleSystem(*R,P,{load.inod,ldof});
      }
#if INT_DEBUG > 1
      std::cout <<"\nExternal load at DOF "<< load.inod <<","<< load.ldof
                <<": "<< P;
#endif
    }
  }

#if INT_DEBUG > 1
  if (!mySprings.empty() || !myLoads.empty())
    std::cout << std::endl;
#endif
  return ok;
}


/*!
  This method is overridden to optionally also compute the reaction forces.
  This requires an additional assembly loop calculating the internal forces,
  since we only are doing a linear solve here.
*/

bool SIMAndesShell::solveSystem (Vector& solution, int printSol, double* rCond,
                                 const char*, size_t)
{
  if (!this->solveEqSystem(solution,0,rCond,printSol))
    return false;
  else if (myRFset.empty())
    return true;

  // Assemble the reaction forces. Strictly, we only need to assemble those
  // elements that have nodes on the Dirichlet boundaries, but...
  int oldlevel = msgLevel;
  msgLevel = 1;
  bool oki = this->assembleForces({solution},Elastic::time,&myReact);
  msgLevel = oldlevel;
  if (!oki) return false;

  IFEM::cout <<"\n >>> Nodal reaction forces <<<\n"
             <<"\n   Internal External               Fx             Fy"
             <<"             Fz             Mx             My             Mz";
  this->printNRforces(this->getNodeSet(myRFset));
  return true;
}


const RealArray* SIMAndesShell::getReactionForces() const
{
  if (!myRFset.empty())
    return myReact.empty() ? nullptr : &myReact;

  return this->Parent::getReactionForces();
}


DataExporter* SIMAndesShell::getHDF5writer (const Vector& psol,
                                            double dumpNodeMap) const
{
  IFEM::cout <<"\nWriting HDF5 file "<< opt.hdf5 <<".hdf5"<< std::endl;

  DataExporter* writer = new DataExporter(true,opt.saveInc);
  writer->registerWriter(new HDF5Writer(opt.hdf5,adm));

  int result = DataExporter::PRIMARY | DataExporter::DISPLACEMENT;
  if (!opt.pSolOnly) result |= DataExporter::SECONDARY;
  if (dumpNodeMap)   result |= DataExporter::L2G_NODE;
  writer->registerField("u","solution",DataExporter::SIM,result);
  writer->setFieldValue("u",this,&psol);

  return writer;
}


bool SIMAndesShell::writeGlvLoc (std::vector<std::string>& locfiles,
                                 bool nodalR, int& nBlock) const
{
  for (const std::string& fName : locfiles)
  {
    ElementBlock* sensorBlock = new ElementBlock(8);
    std::ifstream locs(fName);
    while (locs.good())
    {
      Vec3 XYZloc;
      int idx = 0;
      locs >> idx;
      for (const ASMbase* pch : myModel)
        if (const ASMu2DNastran* p = dynamic_cast<const ASMu2DNastran*>(pch); p)
          if (ElementBlock* sensor = p->sensorGeometry(idx,nodalR); sensor)
          {
            sensorBlock->merge(*sensor,false);
            break;
          }
    }

    if (sensorBlock->getNoElms() < 1)
      delete sensorBlock;
    else if (!this->getVTF()->writeGrid(sensorBlock,fName.c_str(),++nBlock))
      return false;
  }

  return true;
}


bool SIMAndesShell::writeGlvNormal (int& geoBlk, int& nBlock) const
{
  VTF* vtf = this->getVTF();
  if (!vtf || myProblem->hasTractionValues()) return true;

  for (const ASMbase* pch : myModel)
    if (const ASMu2DNastran* shl = dynamic_cast<const ASMu2DNastran*>(pch); shl)
      if (std::vector<Vec3Pair> normals; shl->getShellNormals(normals))
      {
        if (msgLevel > 1)
          IFEM::cout <<"Writing shell normal vectors"<< std::endl;

        // Write shell normals as discrete point vectors to the VTF-file
        return vtf->writeVectors(normals,geoBlk,++nBlock,"Normal vectors",1);
      }

  return true; // No shell elements
}


/*!
  This method is overridden to also write out the sea surface, if any.
*/

bool SIMAndesShell::writeGlvG (int& nBlock, const char* inpFile, bool doClear)
{
  if (!this->Parent::writeGlvG(nBlock,inpFile,doClear))
    return false;

  if (seaGridSize < 1.0e-8 || seaLx < seaGridSize || seaLy < seaGridSize)
    return true; // No sea surface visuzlization

  const size_t nx = seaLx / seaGridSize;
  const size_t ny = seaLy / seaGridSize;
  const double dx = seaLx / nx;
  const double dy = seaLy / ny;

  ElementBlock* seaSurf = new ElementBlock(4);
  seaSurf->resize(nx+1,ny+1);

  size_t i, j, n = 0;
  for (j = 0; j <= ny; j++)
    for (i = 0; i <= nx; i++)
      seaSurf->setCoor(n++, seaX0.x+i*dx, seaX0.y+j*dy, seaX0.z);

  for (j = n = 0; j < ny; j++)
    for (i = 0; i < nx; i++)
    {
      seaSurf->setNode(n++, j*(nx+1)+i);
      seaSurf->setNode(n++, j*(nx+1)+i+1);
      seaSurf->setNode(n++, (j+1)*(nx+1)+i+1);
      seaSurf->setNode(n++, (j+1)*(nx+1)+i);
    }

  seaBlock = ++nBlock;
  return this->getVTF()->writeGrid(seaSurf,"Sea surface",seaBlock);
}


bool SIMAndesShell::writeGlvA (int& nBlock, int iStep, double time, int) const
{
  if (seaBlock < 1 || !seasurf) return true; // no sea surface visualization

  VTF* vtf = this->getVTF();
  const ElementBlock* sea = vtf ? vtf->getBlock(seaBlock) : nullptr;
  if (!sea) return false;

  // Evaluate the sea surface elevation at all grid points
  size_t npt = sea->getNoNodes();
  Matrix zeta(3,npt);
  for (size_t i = 0; i < npt; i++)
    zeta(3,i+1) = (*seasurf)(Vec4(sea->getCoord(i),time));

  // Output as vector field
  if (!vtf->writeVres(zeta,++nBlock,seaBlock))
    return false;

  const_cast<SIMAndesShell*>(this)->addDisBlk[seaBlock] = nBlock;
  return true;
}


/*!
  This method is overridden to also write out the von Mises stresses as a
  piecewise constant element field.
*/

int SIMAndesShell::writeGlvS2 (const Vector& psol, int iStep, int& nBlock,
                               double time, int idBlock, int psolComps)
{
  int idB = this->Parent::writeGlvS2(psol,iStep,nBlock,time,idBlock,psolComps);
  if (idB < idBlock) return idB;

  size_t nf = myProblem->getNoFields(2);
  if (nf != 2 && nf != 18) return idB; // no von Mises stress output

  std::array<Vector,2> elmRes;
  for (Vector& v : elmRes)
    v.reserve(this->getNoElms(false,true));

  for (const ASMbase* pch : myModel)
    if (dynamic_cast<const ASMu2DNastran*>(pch))
    {
      // Direct evaluation of secondary solution variables at element centers
      Matrix field;
      pch->extractNodeVec(psol,myProblem->getSolution());
      if (!pch->evalSolution(field,*myProblem,nullptr,true))
        return -3;

      for (size_t iel = 1; iel <= pch->getNoElms(true); iel++)
        if (pch->getElementNodes(iel).size() > 1) // skip 1-noded mass elements
          for (size_t i = 0; i < 2; i++) // extract bottom/top von Mises stress
            elmRes[i].push_back(field(nf-1+i,iel));
    }
    else for (Vector& v : elmRes)
      v.resize(v.size()+pch->getNoElms(true),utl::RETAIN);

  const char* vms[2] = { "Bottom von Mises stress", "Top von Mises stress" };
  const char** q = vms;
  bool ok = true;
  for (int i = 0; i < 2 && ok; i++, q++)
    ok = this->writeGlvE(elmRes[i],iStep,nBlock,*q,idB++,true);

  return ok ? idB : -5;
}


bool SIMAndesShell::writeAddFuncs (int& nBlock, int& idBlock,
                                   const Vector& psol,
                                   int iStep, double time)
{
  // Find the shell surface pressure function, if any
  for (const ASMbase* pch : myModel)
    if (dynamic_cast<const ASMu2DLag*>(pch) && shellp)
      if (const RealFunc* pressure = shellp->getPressure(pch); pressure)
        // Write the pressure function as a scalar field
        if (!this->writeGlvF(*pressure, "Pressure", iStep, nBlock,
                             psol.empty() ? nullptr : &psol,
                             idBlock++, time, pch))
          return false;

  return this->Parent::writeAddFuncs(nBlock,idBlock,psol,iStep,time);
}


bool SIMAndesShell::writeField (const Vector& field, const std::string& fldName,
                                int iStep, int& nBlock, int idBlock,
                                bool isNodal, bool isLocal) const
{
  if (isLocal && (isNodal || field.size() != this->getNoElms(true)))
  {
    // This field is already assumed to be associated with the
    // nodal points (or elements) of the shell element patch
    int geomID = this->getStartGeo();
    for (const ASMbase* pch : myModel)
      if (dynamic_cast<const ASMu2DLag*>(pch))
      {
        if (msgLevel > 1)
          IFEM::cout <<"Writing "<< (isNodal ? "nodal" : "element")
                     <<" scalar field \""<< fldName
                     <<"\" for patch "<< pch->idx+1 << std::endl;

        VTF* vtf = this->getVTF();
        if (isNodal)
        {
          Matrix sField;
          if (!pch->evalSolution(sField,field,nullptr,false))
            return false;

          if (!vtf->writeNres(sField,++nBlock,++geomID))
            return false;
        }
        else
          if (!vtf->writeEres(field,++nBlock,++geomID))
            return false;

        return vtf->writeSblk(nBlock,fldName.c_str(),idBlock,iStep,!isNodal);
      }
      else
        ++geomID;

    return true; // no shell elements (silently ignore)
  }
  else if (isNodal) // global nodal field
    return this->writeGlvS(field,fldName.c_str(),iStep,nBlock,idBlock);
  else // global element field
    return this->writeGlvE(field,iStep,nBlock,fldName.c_str(),idBlock,true);
}
