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
#include "SIMAndesShell.h"
#include "Profiler.h"
#ifdef HAS_FFLLIB
#include "FFlLib/FFlFEParts/FFlAllFEParts.H"
#endif
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


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
  char* infile = nullptr;

  IFEM::Init(argc,argv,"Linear Elastic Shell solver");

  for (int i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-check"))
      iop = 100;
    else if (!strcmp(argv[i],"-fixDup"))
    {
      fixDup = true;
      if (i+1 < argc && argv[i+1][0] != '-')
        Vec3::comparisonTolerance = atof(argv[++i]);
    }
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-eig <iop> [-nev <nev>] [-ncv <ncv] [-shift <shf>]]\n"
              <<"       [-fixDup [<tol>]] [-vtf <format>] [-check]\n";
    delete prof;
    return 0;
  }

  IFEM::cout <<"\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
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
  SIMoutput* model = new SIMAndesShell();
  std::vector<Mode> modes;

  // Lambda function for cleaning the heap-allocated objects on termination.
  // To ensure that their destructors are invoked also on simulation failure.
  auto&& terminate = [model,prof](int status)
  {
    delete model;
    delete prof;
    exit(status);
  };

  // Read in model definitions
  if (!model->read(infile))
    terminate(1);

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
    if (!model->assembleSystem())
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
