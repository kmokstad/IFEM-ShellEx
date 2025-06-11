// $Id$
//==============================================================================
//!
//! \file HasPointLoads.C
//!
//! \date Jun 12 2025
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class representing a nodal point load container.
//!
//==============================================================================

#include "HasPointLoads.h"
#include "AlgEqSystem.h"
#include "SystemMatrix.h"
#include "SAM.h"
#include "Functions.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml2.h"


HasPointLoads::~HasPointLoads ()
{
  for (PointLoad& load : myLoads)
    delete load.func;
}


bool HasPointLoads::parseLoad (const tinyxml2::XMLElement* elem)
{
  int inod = 0, ldof = 0;
  utl::getAttribute(elem,"node",inod);
  utl::getAttribute(elem,"dof",ldof);
  if (inod <= 0 || ldof <= 0 || !elem->FirstChild())
    return false;

  IFEM::cout <<"\tNode "<< inod <<" dof "<< ldof <<" Load: ";
  std::string type("constant");
  utl::getAttribute(elem,"type",type);
  if (type == "file")
  {
    IntSet dofs = utl::getDigits(ldof);
    char* funcv = strdup(elem->FirstChild()->Value());
    char* fname = strtok(funcv," ");
    char* colmn = nullptr;
    IFEM::cout <<"\""<< fname <<"\"";
    for (int ldof : dofs)
      if ((colmn = strtok(nullptr," ")))
      {
	IFEM::cout <<" ("<< ldof <<","<< colmn <<")";
	myLoads.emplace_back(inod, ldof, new LinearFunc(fname,atoi(colmn)));
      }
    free(funcv);
    IFEM::cout << std::endl;
  }
  else if (ldof <= 6)
  {
    ScalarFunc* f = nullptr;
    if (type == "constant")
    {
      f = new ConstantFunc(atof(elem->FirstChild()->Value()));
      IFEM::cout << (*f)(0.0) << std::endl;
    }
    else
      f = utl::parseTimeFunc(elem->FirstChild()->Value(),type);

    myLoads.emplace_back(inod,ldof,f);
  }
  else
    IFEM::cout <<" (invalid, ignored)."<< std::endl;

  return true;
}


bool HasPointLoads::renumberLoadedNodes (const std::map<int,int>& nodeMap)
{
  bool ok = true;
  for (PointLoad& load : myLoads)
    if (load.inod > 0)
      ok &= utl::renumber(load.inod,nodeMap,true);

  return ok;
}


bool HasPointLoads::assembleLoads (AlgEqSystem& eqSys, double t) const
{
  bool ok = true;
  const size_t nrhs = eqSys.getNoRHS();
  SystemVector* R = nrhs > 0 ? eqSys.getVector(nrhs-1) : nullptr;
  if (!R) return ok; // No external force vector

  for (const PointLoad& load : myLoads)
  {
    double P = (*load.func)(t);
    int ldof = load.ldof;
    eqSys.addScalar(P,ldof-1);
    ok &= eqSys.getSAM().assembleSystem(*R,P,{load.inod,ldof});
  }

  return ok;
}
