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
#include "NonlinearDriver.h"
#include "ElasticityUtils.h"
#include "ElasticityArgs.h"
#include "DynamicSim.h"
#include "DataExporter.h"
#include "HDF5Writer.h"
#include "Profiler.h"
#if INT_DEBUG > 2
#include "SAM.h"
#endif
#ifdef HAS_FFLLIB
#include "FFlLib/FFlFEParts/FFlAllFEParts.H"
#endif
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


/*!
  \brief Reads the input file and invokes the multi-load-case simulation driver.
*/

int mlcSim (char* infile, SIMbase* model, bool fixDup)
{
  IFEM::cout <<"\nUsing the multi-load-case simulation driver."<< std::endl;
  NonlinearDriver simulator(*model,true);

  // Read in solver and model definitions
  if (!simulator.read(infile))
    return 1;

  // Let the stop time specified on command-line override input file setting
  if (Elastic::time > 1.0)
    simulator.setStopTime(Elastic::time);

  model->opt.print(IFEM::cout,true) << std::endl;

  utl::profiler->stop("Model input");

  // Preprocess the model and establish data structures for the algebraic system
  if (!model->preprocess({},fixDup))
    return 2;

  // Save FE model to VTF file for visualization
  if (model->opt.format >= 0 && !simulator.saveModel(infile))
    return 3;

  // Initialize the solution vectors
  simulator.initPrm();
  simulator.initSol();

  // Initialize the linear equation solver
  if (!simulator.initEqSystem())
    return 3;

  // Open HDF5 result database, if requested
  DataExporter* writer = nullptr;
  if (model->opt.dumpHDF5(infile))
  {
    const std::string& fileName = model->opt.hdf5;
    IFEM::cout <<"\nWriting HDF5 file "<< fileName <<".hdf5"<< std::endl;

    writer = new DataExporter(true,model->opt.saveInc);
    writer->registerWriter(new HDF5Writer(fileName,model->getProcessAdm()));

    int results = DataExporter::PRIMARY | DataExporter::DISPLACEMENT;
    writer->registerField("u","solution",DataExporter::SIM,results);
    writer->setFieldValue("u",model,&simulator.getSolution());
  }

  // Now invoke the main solution driver
  int status = simulator.solveProblem(writer);
  delete writer;
  return status;
}


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
  \arg -check : Data check only, read model and output to VTF (no solution)
  \arg -fixDup \a tol : Resolve co-located nodes by merging them into one
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -vtfres \a file1 \a file2 ... : Extra files for direct VTF output
  \arg -vtfgrp \a file1 \a file2 ... : Extra files for element set visualisation
  \arg -refsol \a file1 \a file2 ... : Files with reference solution
*/

int main (int argc, char** argv)
{
  Profiler* prof = new Profiler(argv[0]);
  utl::profiler->start("Initialization");

  int iop = 0;
  bool fixDup = false;
  bool mlcase = false;
  bool nodalR = false;
  char dynSol = false;
  char* infile = nullptr;
  ElasticityArgs args;
  std::vector<std::string> resfiles, grpfiles, disfiles;

  IFEM::Init(argc,argv,"Linear Elastic Shell solver");

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
    else if (!strcmp(argv[i],"-fixDup"))
    {
      fixDup = true;
      if (i+1 < argc && argv[i+1][0] != '-')
        Vec3::comparisonTolerance = atof(argv[++i]);
    }
    else if (!strncmp(argv[i],"-mlc",4))
      mlcase = true;
    else if (!strncmp(argv[i],"-qstat",6))
      dynSol = 's';
    else if (!strncmp(argv[i],"-dyn",4))
      dynSol = 'd';
    else if (!strncmp(argv[i],"-vtfres",6))
    {
      while (i+1 < argc && argv[i+1][0] != '-')
        if (!strcmp(argv[++i],"nodal"))
          nodalR = true;
        else
          resfiles.push_back(argv[i]);
    }
    else if (!strncmp(argv[i],"-vtfgrp",6))
      while (i+1 < argc && argv[i+1][0] != '-')
        grpfiles.push_back(argv[++i]);
    else if (!strncmp(argv[i],"-refsol",6))
      while (i+1 < argc && argv[i+1][0] != '-')
        disfiles.push_back(argv[++i]);
    else if (!infile)
    {
      infile = argv[i];
      if (strcasestr(infile,".xinp"))
      {
        if (!args.readXML(infile,false))
          return 1;
        i = 0; // start over and let command-line options override input file
      }
    }
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-eig <iop> [-nev <nev>] [-ncv <ncv] [-shift <shf>]]\n"
              <<"       [-free] [-time <t>] [-mlc|-qstatic|-dynamic] [-check]\n"
              <<"       [-vtf <format> [-vtfres <files>] [-vtfgrp <files>]]\n"
              <<"       [-fixDup [<tol>]] [-refsol <files>]\n";
    delete prof;
    return 0;
  }

  if (IFEM::getOptions().eig < 0)
    args.eig = 0;
  else if (IFEM::getOptions().eig > 0)
    args.eig = IFEM::getOptions().eig;

  // Boundary conditions can be ignored only in generalized eigenvalue analysis
  if (args.eig < 3) SIMbase::ignoreDirichlet = false;

  bool modalS = dynSol && args.eig >= 3 && args.eig != 5; // Modal solution

  IFEM::cout <<"\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
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
  SIMoutput* model = modalS ? new SIMShellModal(modes) : new SIMAndesShell();

  // Lambda function for cleaning the heap-allocated objects on termination.
  // To ensure that their destructors are invoked also on simulation failure.
  auto&& terminate = [model,prof,dynSol](int status)
  {
    if (status > 10 && !dynSol)
      utl::profiler->stop("Postprocessing");
    delete model;
    delete prof;
    exit(status);
  };

  if (mlcase) // Solve the multi-load-case linear static problem
    terminate(mlcSim(infile,model,fixDup));

  if (dynSol && !modalS) // Invoke the linear Newmark time integration simulator
    terminate(dynamicSim(infile,model,fixDup));

  // Read in model definitions
  if (!model->read(infile))
    terminate(1);

  model->opt.print(IFEM::cout,true) << std::endl;

  utl::profiler->stop("Model input");

  // Establish the FE data structures
  if (!model->preprocess({},fixDup))
    terminate(2);

  Vectors displ(2);
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

  switch (iop+model->opt.eig) {
  case 0:
  case 200:
    // Static solution: Assemble [Km] and {R}
    model->setMode(SIM::STATIC);
    model->setQuadratureRule(2);
    model->initSystem(model->opt.solver);
    if (!model->assembleSystem(Elastic::time))
      terminate(4);

    // Solve the linear system of equations
    if (iop == 200)
      model->dumpEqSys(); // No solution, just dump the system matrices to file
    else if (!model->solveSystem(displ,1))
      terminate(5);
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
      terminate(8);

    if (!model->systemModes(modes))
      terminate(9);
    break;

  default:
    // Free vibration: Assemble [Km] and [M]
    model->setMode(SIM::VIBRATION);
    model->setQuadratureRule(2);
    model->initSystem(model->opt.solver,2,0);
    if (!model->assembleSystem())
      terminate(8);

    // Solve the generalized eigenvalue problem
    if (iop == 200)
      model->dumpEqSys(); // No solution, just dump the system matrices to file
    else if (!model->systemModes(modes))
      terminate(9);
  }

  if (modalS) // Solve the dynamics problem using modal transformation
    terminate(modalSim(infile,modes.size(),false,dynSol=='s',model));

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

    // Write VTF-file with model geometry
    if (!model->writeGlvG(geoBlk,infile))
      terminate(12);

    // Write Dirichlet boundary conditions
    if (!model->writeGlvBC(nBlock))
      terminate(13);

    // Write global node numbers as scalar fields
    if (!model->writeGlvNo(nBlock))
      terminate(14);

    Vector data(model->getNoElms());
    static_cast<SIMAndesShell*>(model)->getShellThicknesses(data);
    if (!model->writeGlvE(data,iStep,nBlock,"Shell thickness",11,true))
      terminate(15);

    // Write solution fields to VTF-file
    model->setMode(SIM::RECOVERY);
    if (!model->writeGlvS(displ.front(),iStep,nBlock,time,nullptr,12))
      terminate(16);

    // Write reference solution, if any
    if (model->writeGlvS1(displ.back(),iStep,nBlock,time,"Reference",20,-1) < 0)
      terminate(17);

    // Write eigenmodes
    for (const Mode& mode : modes)
      if (!model->writeGlvM(mode,true,nBlock))
        terminate(18);

    if (nodalR && !(grpfiles.empty() && resfiles.empty()))
      data.resize(model->getNoNodes());

    bool ok = true;
    int idBlock = 100;
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
      if (!ok) terminate(19);
    }

    model->writeGlvStep(iStep, time, resfiles.empty() ? -1 : 0);

    if (!resfiles.empty())
    {
      // Write external  results (from OSP calculations)
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
        if (!ok) terminate(19);

        model->writeGlvStep(iStep+1,times[iStep-1]);
      }
    }
  }

  utl::profiler->stop("Postprocessing");
  terminate(0);
}
