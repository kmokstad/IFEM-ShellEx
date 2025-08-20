// $Id$
//==============================================================================
//!
//! \file main.C
//!
//! \date Feb 20 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Main program for the linear elastic shell solver.
//!
//==============================================================================

#include "IFEM.h"
#include "SIMenums.h"
#include "SIMShellModal.h"
#include "SIMAndesSplit.h"
#include "ElasticityUtils.h"
#include "ElasticityArgs.h"
#include "MultiLoadCaseSim.h"
#include "DynamicSim.h"
#include "DataExporter.h"
#include "Profiler.h"
#ifdef INT_DEBUG
#include "SAM.h"
#endif
#ifdef HAS_FFLLIB
#include "FFlLib/FFlFEParts/FFlAllFEParts.H"
#endif
#include <array>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

namespace ASM { extern double Ktra, Krot; extern bool skipVTFmass; }


/*!
  \brief Main program for the linear elastic shell solver.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -eig \a iop : Eigenproblem solver to use (1...6)
  \arg -nev \a nev : Number of eigenvalues to compute
  \arg -ncv \a ncv : Number of Arnoldi vectors to use in the eigenvalue analysis
  \arg -shift \a shf : Shift value to use in the eigenproblem solver
  \arg -free : Ignore all boundary conditions (use in free vibration analysis)
  \arg -time : Time for evaluation of possible time-dependent functions
  \arg -mlc : Solve the linear static problem as a multi-load-case problem
  \arg -qstatic : Solve the linear dynamics problem as quasi-static
  \arg -dynamic : Solve the linear dynamics problem using modal transformation
  \arg -modes : Create time history as a linear combination of mode shapes
  \arg -check : Data check only, read model and output to VTF (no solution)
  \arg -ignoreSol : Assembly only and output to VTF (skip solution)
  \arg -fixDup \a tol : Resolve co-located nodes by merging them into one
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -vtfres \a file1 \a file2 ... : Extra files for direct VTF output
  \arg -vtfgrp \a file1 \a file2 ... : Extra files for element set visualisation
  \arg -vtfloc \a file1 \a file2 ... : Extra files for sensor visualisation
  \arg -refsol \a file1 \a file2 ... : Files with reference solution
  \arg -vizRHS : Save the right-hand-side load vector on the VTF-file
  \arg -hdf5 : Write primary and secondary solution to HDF5 file
  \arg -noSets : Ignore Nastran SET definitions
  \arg -noBeams : Ignore beam elements
  \arg -noEccs : Ignore beam end offsets
  \arg -noRBE3 : Replace RBE3 elements by equivalent RBE2 elements
  \arg -dumpNodeMap : Dump Local-to-global node number mapping to HDF5
  \arg -split : Split the model into two material regions
  \arg -keep-previous-state : Use previous state when evaluating the
       state-dependent property functions
  \arg -no-vtfmass : Skip the sphere geometries for 1-noded mass elements
  \arg -Kbush : Spring stiffness(es) connecting RBE3 reference nodes to ground
*/

int main (int argc, char** argv)
{
  Profiler* prof = new Profiler(argv[0]);
  utl::profiler->start("Initialization");

  int iop = 0;
  bool vizRHS = false;
  bool fixDup = false;
  bool mlcase = false;
  char nodalR = false;
  bool splitM = false;
  char dynSol = false;
  bool eigSim = false;
  unsigned short int nstates = 0;
  bool dumpNodeMap = false;
  char* infile = nullptr;
  char* bdfile = nullptr;
  ElasticityArgs args;
  std::vector<std::string> resfiles, grpfiles, locfiles, disfiles;

  IFEM::Init(argc,argv,"Linear Elastic Shell solver");
  ASM::cachePolicy = ASM::NO_CACHE;

  // Lambda function checking for valid Nastran bulkk data file extension.
  auto&& isBDF = [](const char* fname)
  {
    static const char* bdfext[] = { ".dat", ".nas", ".bdf", NULL };
    for (const char** ext = bdfext; *ext; ++ext)
      if (strcasestr(fname,*ext)) return true;
    return false;
  };

  for (int i = 1; i < argc; i++)
    if (argv[i] == infile || args.parseArg(argv[i]))
      ; // ignore the input file on the second pass
    else if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-time") && i < argc-1)
      Elastic::time = atof(argv[++i]);
    else if (!strcmp(argv[i],"-free"))
      SIMbase::ignoreDirichlet = true;
    else if (!strcmp(argv[i],"-check"))
      iop = 100;
    else if (!strcmp(argv[i],"-ignoreSol"))
      iop = 200;
    else if (!strcmp(argv[i],"-vizRHS"))
      vizRHS = true;
    else if (!strcmp(argv[i],"-noSets"))
      SIMAndesShell::readSets = false;
    else if (!strcmp(argv[i],"-noBeams"))
      SIMAndesShell::useBeams = 0;
    else if (!strcmp(argv[i],"-noEccs"))
      SIMAndesShell::useBeams = 2;
    else if (!strcmp(argv[i],"-noRBE3"))
      SIMAndesShell::replRBE3 = true;
    else if (!strcmp(argv[i],"-Kbush"))
    {
      if (i+1 < argc && argv[i+1][0] != '-')
        ASM::Ktra = ASM::Krot = atof(argv[++i]);
      if (i+1 < argc && argv[i+1][0] != '-')
        ASM::Krot = atof(argv[++i]);
    }
    else if (!strcmp(argv[i],"-fixDup"))
    {
      fixDup = true;
      if (i+1 < argc && argv[i+1][0] != '-')
        Vec3::comparisonTolerance = atof(argv[++i]);
    }
    else if (!strncmp(argv[i],"-mlc",4))
      mlcase = true;
    else if (!strncmp(argv[i],"-modes",6))
      eigSim = true;
    else if (!strncmp(argv[i],"-qstat",6))
    {
      dynSol = 's';
      if (nstates == 0) nstates = 1;
    }
    else if (!strncmp(argv[i],"-dyn",4))
    {
      dynSol = 'd';
      if (nstates == 0) nstates = 1;
    }
    else if (!strcmp(argv[i],"-keep-previous-state"))
      nstates = 2; // both current and previous state will reside in core
    else if (!strcmp(argv[i],"-no-vtfmass"))
      ASM::skipVTFmass = true;
    else if (!strncmp(argv[i],"-vtfres",6))
    {
      while (i+1 < argc && argv[i+1][0] != '-')
        if (!strcmp(argv[++i],"nodal"))
          nodalR = 'd';
        else if (!strcmp(argv[i],"modes"))
          nodalR = 'm';
        else
          resfiles.push_back(argv[i]);
    }
    else if (!strncmp(argv[i],"-vtfgrp",6))
      while (i+1 < argc && argv[i+1][0] != '-')
        grpfiles.push_back(argv[++i]);
    else if (!strncmp(argv[i],"-vtfloc",6))
      while (i+1 < argc && argv[i+1][0] != '-')
        locfiles.push_back(argv[++i]);
    else if (!strncmp(argv[i],"-refsol",6))
      while (i+1 < argc && argv[i+1][0] != '-')
        disfiles.push_back(argv[++i]);
    else if (!strncmp(argv[i],"-dumpNod",8))
      dumpNodeMap = true;
    else if (!strncmp(argv[i],"-split",6))
      splitM = true;
    else if (!infile && strcasestr(argv[i],".xinp") && !bdfile)
    {
      infile = argv[i];
      if (!args.readXML(infile,false))
        return 1;
      i = 0; // start over and let command-line options override input file
    }
    else if (!bdfile && isBDF(argv[i]) && !infile)
      bdfile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (bdfile)
  {
    // Create a temporary input file referring to the Nastran bulk data file
    static const char* tmpname = "/tmp/tmp.xinp";
    infile = const_cast<char*>(tmpname);
    std::ofstream xinp(infile);
    xinp <<"<simulation><geometry dim=\"3\">\n"
         <<"  <patchfile type=\"Nastran\">"<< bdfile <<"</patchfile>\n"
         <<"</geometry></simulation>\n";

    if (IFEM::getOptions().format >= 0)
    {
      // Set default vtf file name
      std::string& vtf = IFEM::getOptions().vtf;
      vtf = bdfile;
      vtf.replace(vtf.find_last_of('.'),std::string::npos,".vtf");
    }
  }

  if (!infile)
  {
    // Lambda function for nicely print of usage.
    auto&& showUsage = [argv](const std::vector<const char*>& args)
    {
      const size_t width = 80;
      size_t col = 7 + strlen(argv[0]);
      std::cout <<"usage: "<< argv[0];
      for (const char* arg : args)
      {
        size_t w = strlen(arg);
        if (col+w <= width)
        {
          std::cout <<" "<< arg;
          col += 1+w;
        }
        else
        {
          std::cout <<"\n       "<< arg;
          col = 7+w;
        }
      }
      std::cout << std::endl;
    };

    showUsage({"<inputfile>","[-dense|-spr|-superlu[<nt>]|-samg|-petsc]",
               "[-eig <iop> [-nev <nev>] [-ncv <ncv] [-shift <shf>]]",
               "[-free]","[-time <t>]","[-check]","[-ignoreSol]",
               "[-mlc|-qstatic|-dynamic|-modes]","[-keep-previous-state]",
               "[-vtf <format> [-vtfres <files>] [-vtfgrp <files>] [-vizRHS] "
               "[-no-vtfmass]]",
               "[-hdf5 [<filename>] [-dumpNodeMap]]","[-fixDup [<tol>]]",
               "[-refsol <files>]","[-noBeams]","[-noEccs]","[-noSets]",
               "[-noRBE3]","[-split]"});
    delete prof;
    return 0;
  }

  if (IFEM::getOptions().eig < 0)
    args.eig = 0;
  else if (IFEM::getOptions().eig > 0)
    args.eig = IFEM::getOptions().eig;

  // Boundary conditions can be ignored only in generalized eigenvalue analysis
  if (args.eig < 3) SIMbase::ignoreDirichlet = false;

  bool modal = dynSol && args.eig >= 3 && args.eig != 5; // Modal solution
  if (dynSol || args.eig < 3) eigSim = false;

  IFEM::cout <<"\nInput file: "<< infile;
  if (!mlcase && !dynSol)
    IFEM::cout <<"\nEvaluation time for property functions: "<< Elastic::time;
  else if (Elastic::time > 1.0)
    IFEM::cout <<"\nSimulation stop time: "<< Elastic::time;
  if (SIMbase::ignoreDirichlet)
    IFEM::cout <<"\nSpecified boundary conditions are ignored";
  if (fixDup)
    IFEM::cout <<"\nCo-located nodes will be merged,"
               <<" using comparison tolerance "<< Vec3::comparisonTolerance
               << std::endl;

#ifdef HAS_FFLLIB
  FFl::initAllElements();
#endif
  utl::profiler->stop("Initialization");
  utl::profiler->start("Model input");

  // Create the simulation model
  std::vector<Mode> modes;
  DataExporter* writer = nullptr;
  SIMAndesShell* model = nullptr;
  if (modal)
    model = new SIMShellModal(modes);
  else if (splitM)
    model = new SIMAndesSplit(nstates);
  else
    model = new SIMAndesShell(nstates);

  // Lambda function for cleaning the heap-allocated objects on termination.
  // To ensure that their destructors are invoked also on simulation failure.
  auto&& terminate = [model,writer,prof,dynSol](int status, bool relFFl = false)
  {
    if (status > 10 && !dynSol)
      utl::profiler->stop("Postprocessing");
#ifdef HAS_FFLLIB
    if (relFFl)
      FFl::releaseAllElements();
#endif
    delete model;
    delete writer;
    delete prof;
    return status;
  };

  if (mlcase) // Solve the multi-load-case linear static problem
    return terminate(mlcSim(infile,model,fixDup,dumpNodeMap),true);

  if (dynSol && !modal) // Invoke the linear Newmark time integration simulator
    return terminate(dynamicSim(infile,model,fixDup),true);

  if (eigSim) // Create a time history from the mode shapes
    return terminate(modeHistSim(infile,model));

  // Read in model definitions
  if (!model->read(infile))
    return terminate(1,true);

  utl::profiler->stop("Model input");
#ifdef HAS_FFLLIB
  FFl::releaseAllElements();
#endif

  model->opt.print(IFEM::cout,true) << std::endl;

  // Establish the FE data structures
  if (!model->preprocess({},fixDup))
    return terminate(2);

  std::array<Vector,2> displ;
  if (!disfiles.empty())
  {
    // Read reference solution from external file(s) into displ.back()
    size_t incd = disfiles.size() == 1 ? 1 : 6;
    size_t ndof = model->getNoDOFs();
    displ.back().resize(ndof);
    for (size_t ldof = 0; ldof < disfiles.size() && ldof < 6; ldof++)
    {
      std::ifstream ifs(disfiles[ldof]);
      for (size_t idof = ldof; ifs.good() && idof < ndof; idof += incd)
        ifs >> displ.back()[idof];
    }
  }

  // Open HDF5 result database, if requested
  if (model->opt.dumpHDF5(infile))
    writer = model->getHDF5writer(displ.front(),dumpNodeMap);

  Vector load;
  RealArray Fex;
  switch (iop+model->opt.eig) {
  case 0:
  case 200:
    // Static solution: Assemble [Km] and {R}
    model->initForSingleStep();
    model->setMode(SIM::STATIC);
    model->setQuadratureRule(2,true);
    model->initSystem(model->opt.solver,1,1,3);
    if (!model->assembleSystem(Elastic::time))
    {
      if (model->opt.format < 0)
        return terminate(4);
      else
        break;
    }
    if (vizRHS)
    {
      model->extractLoadVec(load,0,"external load");
#ifdef INT_DEBUG
      model->getSAM()->printVector(std::cout,load,"\nLoad vector");
#endif
    }
    if (model->extractScalars(Fex))
    {
      IFEM::cout <<"\nTotal external load: Sum(Fex) =";
      for (double f : Fex) IFEM::cout <<" "<< utl::trunc(f);
      IFEM::cout << std::endl;
    }

    // Solve the linear system of equations
    if (iop == 200)
      model->dumpEqSys(); // No solution, just dump the system matrices to file
    else if (!model->solveSystem(displ.front(),1))
      return terminate(5);
    break;

  case 100:
#if INT_DEBUG > 2
    model->getSAM()->print(std::cout);
#endif
    break; // Model check

  case 1:
  case 2:
    // Assemble and solve the regular eigenvalue problem
    model->setMode(SIM::STIFF_ONLY);
    model->setQuadratureRule(2);
    model->initSystem(model->opt.solver,1,0);
    if (!model->assembleSystem())
      return terminate(8);

    if (!model->systemModes(modes))
      return terminate(9);
    break;

  default:
    // Free vibration: Assemble [Km] and [M]
    model->setMode(SIM::VIBRATION);
    model->setQuadratureRule(2);
    model->initSystem(model->opt.solver,2,0);
    if (!model->assembleSystem())
      return terminate(8);

    // Solve the generalized eigenvalue problem
    if (iop == 200)
      model->dumpEqSys(); // No solution, just dump the system matrices to file
    else if (!model->systemModes(modes))
      return terminate(9);
  }

  if (modal) // Solve the dynamics problem using modal transformation
    return terminate(modalSim(infile,modes.size(),false,dynSol=='s',model));

  utl::profiler->start("Postprocessing");

  if (!displ.front().empty() && !displ.back().empty())
  {
    // Calculate the L2-norm of the two solutions and the difference
    size_t ndof = model->getNoDOFs();
    std::array<Vector,3> tra, rot;
    for (int i = 0; i < 2; i++)
    {
      tra[i].resize(ndof/2);
      rot[i].resize(ndof/2);
      size_t jdof = 1;
      for (size_t idof = 1; idof <= ndof; idof += 6)
        for (int j = 0; j < 3; j++, jdof++)
        {
          tra[i](jdof) = displ[i](idof+j);
          rot[i](jdof) = displ[i](idof+3+j);
        }
    }
    tra[2] = tra[0] - tra[1];
    rot[2] = rot[0] - rot[1];
    IFEM::cout <<"\nL2-norm of this solution: "
               << tra[0].norm2() <<" "<< rot[0].norm2()
               <<"\nL2-norm of ref. solution: "
               << tra[1].norm2() <<" "<< rot[1].norm2()
               <<"\nL2-norm of difference:  "
               << tra[2].norm2() <<" "<< rot[2].norm2() <<" ("
               << 100.0*tra[2].norm2()/tra[1].norm2() <<"% "
               << 100.0*rot[2].norm2()/rot[1].norm2() <<"%)"<< std::endl;
  }

  if (model->opt.format >= 0)
  {
    int geoBlk = 0, nBlock = 0;
    size_t iStep = 1, nStep = 0;
    double time = 0.0;
    Vector data;

    // Write VTF-file with model geometry
    if (!model->writeGlvG(geoBlk,infile))
      return terminate(12);

    // Write sensor locations, if any
    if (nodalR != 'm' && !model->writeGlvLoc(locfiles,nodalR,geoBlk))
      return terminate(12);

    // Write surface pressures, if any
    if (!model->writeGlvT(iStep,geoBlk,nBlock))
      return terminate(13);

    // Write Dirichlet boundary conditions
    if (!model->writeGlvBC(nBlock))
      return terminate(13);

    // Write global node and element numbers as scalar fields
    int idBlock = 7;
    if (!model->writeGlvNo(nBlock,idBlock,18))
      return terminate(14);

    // Write shell thickness as scalar field
    model->getShellThicknesses(data);
    if (!model->writeGlvE(data,iStep,nBlock,"Shell thickness",idBlock++,true))
      return terminate(15);

    // Write load vector to VTF-file
    if (vizRHS && !model->writeGlvV(load,"Load vector",iStep,nBlock,1))
      return terminate(20);

    // Write solution fields to VTF-file
    model->setMode(SIM::RECOVERY);
    if (!model->writeGlvS(displ.front(),iStep,nBlock,time,nullptr,idBlock))
      return terminate(16);

    // Write reference solution, if any
    if (model->writeGlvS1(displ.back(),iStep,nBlock,time,
                          "Reference",idBlock,-1) < 0)
      return terminate(17);

    if (modes.empty() && resfiles.size() == 1 && nodalR == 'm')
    {
      // Read mode shapes from external file
      std::string firstLine;
      std::ifstream ifs(resfiles.front());
      if (std::getline(ifs,firstLine))
      {
        // Use the first line to count the number of modes
        RealArray values;
        char* endPt = const_cast<char*>(firstLine.c_str());
        while (endPt && strlen(endPt) > 0)
          values.push_back(strtod(endPt,&endPt));
        modes.resize(values.size());
        size_t iDof = 0;
        for (Mode& mode : modes)
        {
          mode.eigVec.resize(model->getNoDOFs());
          mode.eigVec[0] = values[iDof++];
          mode.eigVal = mode.eigNo = iDof;
        }
        // Now read the rest of the file - one line for each DOF
        for (iDof = 1; ifs.good() && iDof < modes.front().eigVec.size(); iDof++)
          for (Mode& mode : modes)
            ifs >> mode.eigVec[iDof];
      }
      resfiles.clear();
      nodalR = false;
    }

    // Write eigenmodes
    for (const Mode& mode : modes)
      if (!model->writeGlvM(mode,true,nBlock))
        return terminate(18);

    if (!grpfiles.empty() || !resfiles.empty())
      data.resize(nodalR ? model->getNoNodes() : model->getNoElms(true));

    bool ok = true;
    idBlock = 100;
    for (const std::string& fName : grpfiles)
    {
      // Write node/element set definition
      data.fill(0.0);
      std::ifstream ifs(fName);
      while (ifs.good())
      {
        int iel = 0;
        ifs >> iel;
        if (ifs.good() && iel > 1 && iel <= static_cast<int>(data.size()))
          data[iel-1] = 1.0;
      }

      if (nodalR)
        ok = model->writeGlvS(data,fName.c_str(),iStep,nBlock,++idBlock);
      else
        ok = model->writeGlvE(data,iStep,nBlock,fName.c_str(),++idBlock,true);
      if (!ok) return terminate(19);
    }

    model->writeGlvStep(iStep, time, resfiles.empty() ? -1 : 0);

    if (!resfiles.empty())
    {
      // Write external results (from OSP calculations)
      std::vector<Vectors> extResults(resfiles.size());
      std::vector<double> times;
      for (size_t i = 0; i < resfiles.size(); i++)
      {
        std::ifstream ifs(resfiles[i]);
        ifs >> time;
        while (ifs.good())
        {
          if (i == 0) times.push_back(time);
          for (double& val : data) ifs >> val;
          extResults[i].push_back(data);
          ifs >> time;
        }
        if (!nStep || nStep > extResults[i].size())
          nStep = extResults[i].size();
      }

      IFEM::cout <<"Writing "<< nStep <<" steps"<< std::endl;
      for (iStep = 1; iStep <= nStep; iStep++)
      {
        int jdBlock = idBlock;
        for (size_t j = 0; j < extResults.size() && ok; j++)
          if (nodalR)
            ok = model->writeGlvS(extResults[j][iStep-1],
                                  resfiles[j].c_str(),iStep+1,nBlock,++jdBlock);
          else
            ok = model->writeGlvE(extResults[j][iStep-1],iStep+1,nBlock,
                                  resfiles[j].c_str(),++jdBlock,true);
        if (!ok) return terminate(19);

        model->writeGlvStep(iStep+1,times[iStep-1]);
      }
    }
    model->closeGlv();
  }

  if (writer)
    writer->dumpTimeLevel();

  utl::profiler->stop("Postprocessing");
  return terminate(0);
}
