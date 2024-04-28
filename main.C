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
#include "ModalSim.h"
#include "Profiler.h"
#ifdef HAS_FFLLIB
#include "FFlLib/FFlFEParts/FFlAllFEParts.H"
#endif
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


/*!
  \brief Reads the input file and invokes the multi-load-case simulation driver.
*/

int mlcSim (SIMbase* model, char* infile)
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
  if (!model->preprocess())
    return 2;

  if (model->opt.format >= 0)
  {
    // Save FE model to VTF file for visualization
    model->opt.nViz[2] = 1;
    if (!simulator.saveModel(infile))
      return 3;
  }

  // Initialize the solution vectors
  simulator.initSol();

  // Initialize the linear equation solver
  if (!simulator.initEqSystem())
    return 3;

  // Now invoke the main solution driver
  return simulator.solveProblem(nullptr,nullptr,nullptr,false,0.0,1.0e-6,0);
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
  \arg -fixDup <tol> : Resolve co-located nodes by merging them into one
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
*/

int main (int argc, char** argv)
{
  Profiler* prof = new Profiler(argv[0]);
  utl::profiler->start("Initialization");

  int iop = 0;
  bool fixDup = false;
  bool mlcase = false;
  char dynSol = false;
  char* infile = nullptr;

  IFEM::Init(argc,argv,"Linear Elastic Shell solver");

  for (int i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-time") && i < argc-1)
      Elastic::time = atof(argv[++i]);
    else if (!strcmp(argv[i],"-free"))
      SIMbase::ignoreDirichlet = true;
    else if (!strcmp(argv[i],"-check"))
      iop = 100;
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
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-eig <iop> [-nev <nev>] [-ncv <ncv] [-shift <shf>]]\n"
              <<"       [-free] [-time <t>] [-mlc|-qstatic|-dynamic] [-check]\n"
              <<"       [-fixDup [<tol>]] [-vtf <format>]\n";
    delete prof;
    return 0;
  }

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
  SIMoutput* model = dynSol ? new SIMShellModal(modes) : new SIMAndesShell();

  // Lambda function for cleaning the heap-allocated objects on termination.
  // To ensure that their destructors are invoked also on simulation failure.
  auto&& terminate = [model,prof](int status)
  {
    delete model;
    delete prof;
    exit(status);
  };

  if (mlcase) // Solve the multi-load-case linear static problem
    terminate(mlcSim(model,infile));

  // Read in model definitions
  if (!model->read(infile))
    terminate(1);

  if (model->opt.eig != 3 && model->opt.eig != 4 && model->opt.eig != 6)
    dynSol = false; // Dynamics solution requires eigenmode calculation

  model->opt.print(IFEM::cout,true) << std::endl;

  utl::profiler->stop("Model input");

  // Establish the FE data structures
  if (!model->preprocess({},fixDup))
    terminate(2);

  Vectors displ(1);

  switch (iop+model->opt.eig) {
  case 0:
    // Static solution: Assemble [Km] and {R}
    model->setMode(SIM::STATIC);
    model->setQuadratureRule(2);
    model->initSystem(model->opt.solver);
    if (!model->assembleSystem(Elastic::time))
      terminate(4);

    // Solve the linear system of equations
    if (!model->solveSystem(displ,1))
      terminate(5);
    break;

  case 3:
  case 4:
  case 6:
    // Free vibration: Assemble [Km] and [M]
    model->setMode(SIM::VIBRATION);
    model->setQuadratureRule(2);
    model->initSystem(model->opt.solver,2,0);
    if (!model->assembleSystem())
      terminate(8);

    // Solve the generalized eigenvalue problem
    if (!model->systemModes(modes))
      terminate(9);
  }

  if (dynSol) // Solve the dynamics problem using modal transformation
    terminate(modalSim(infile,modes.size(),false,dynSol=='s',model));

  utl::profiler->start("Postprocessing");

  if (model->opt.format >= 0)
  {
    int geoBlk = 0, nBlock = 0;

    // Write VTF-file with model geometry
    if (!model->writeGlvG(geoBlk,infile))
      terminate(12);

    // Write Dirichlet boundary conditions
    if (!model->writeGlvBC(nBlock))
      terminate(13);

    // Write global node numbers as scalar fields
    if (!model->writeGlvNo(nBlock))
      terminate(14);

    // Write solution fields to VTF-file
    model->setMode(SIM::RECOVERY);
    if (!model->writeGlvS(displ[0],1,nBlock))
      terminate(16);

    // Write eigenmodes
    for (const Mode& mode : modes)
      if (!model->writeGlvM(mode,true,nBlock))
        terminate(18);

    model->writeGlvStep(1,0.0,-1);
  }
  model->closeGlv();

  utl::profiler->stop("Postprocessing");
  terminate(0);
}
