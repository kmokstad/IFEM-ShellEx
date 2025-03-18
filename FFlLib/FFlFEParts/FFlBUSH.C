// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlBUSH.H"
#include "FFlLib/FFlFEParts/FFlPBUSHECCENT.H"
#include "FFlLib/FFlFEParts/FFlPCOORDSYS.H"
#include "FFlLib/FFlFEParts/FFlPORIENT.H"
#include "FFaLib/FFaAlgebra/FFaMat34.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"

//! \brief Convenience macro for dynamic casting of an element attribute pointer
#define GET_ATTRIBUTE(att) dynamic_cast<FFl##att*>(this->getAttribute(#att))

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


void FFlBUSH::init()
{
  using TypeInfoSpec   = FFaSingelton<FFlTypeInfoSpec,FFlBUSH>;
  using AttributeSpec  = FFaSingelton<FFlFEAttributeSpec,FFlBUSH>;
  using ElementTopSpec = FFaSingelton<FFlFEElementTopSpec,FFlBUSH>;

  TypeInfoSpec::instance()->setTypeName("BUSH");
  TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::OTHER_ELM);

  ElementFactory::instance()->registerCreator(TypeInfoSpec::instance()->getTypeName(),
                                              FFaDynCB2S(FFlBUSH::create,int,FFlElementBase*&));

  AttributeSpec::instance()->addLegalAttribute("PBUSHCOEFF",false);
  AttributeSpec::instance()->addLegalAttribute("PBUSHECCENT",false);
  AttributeSpec::instance()->addLegalAttribute("PBUSHORIENT",false);
  AttributeSpec::instance()->addLegalAttribute("PORIENT",false);
  AttributeSpec::instance()->addLegalAttribute("PCOORDSYS",false);

  ElementTopSpec::instance()->setNodeCount(2);
  ElementTopSpec::instance()->setNodeDOFs(6);
}


int FFlBUSH::getNodalCoor(double* X, double* Y, double* Z) const
{
  int ierr = this->FFlElementBase::getNodalCoor(X,Y,Z);
  if (ierr < 0) return ierr;

  X[2] = X[1];
  Y[2] = Y[1];
  Z[2] = Z[1];
  X[1] = X[0];
  Y[1] = Y[0];
  Z[1] = Z[0];

  // Get bushing eccentricity, if any
  FFlPBUSHECCENT* curEcc = GET_ATTRIBUTE(PBUSHECCENT);
  if (!curEcc) return ierr;

  FaVec3 e1 = curEcc->offset.getValue();
  X[0] += e1[0];
  Y[0] += e1[1];
  Z[0] += e1[2];

  return ierr;
}


bool FFlBUSH::getLocalSystem(double* Tlg) const
{
  FaMat34 Tmat;
  FFlPCOORDSYS* cs = GET_ATTRIBUTE(PCOORDSYS);
  if (cs)
    Tmat.makeCS_Z_XZ(cs->Origo.getValue(),
                     cs->Zaxis.getValue(),
                     cs->XZpnt.getValue());
  else
  {
    FFlPORIENT* po = GET_ATTRIBUTE(PORIENT);
    if (!po) // If no PORIENT, try the equivalent old name also
      po = dynamic_cast<FFlPORIENT*>(this->getAttribute("PBUSHORIENT"));

    if (po && !po->directionVector.getValue().isZero())
    {
      int nnod = 0;
      std::array<FaVec3,2> X;
      for (const NodeRef& node : myNodes)
        if (node.isResolved() && nnod < 2)
          X[nnod++] = node->getPos();
        else
          ListUI <<" *** Element node "<< nnod+1 <<" ("<< node->getID()
                 <<") is not resolved for BUSH element "<< this->getID()
                 <<".\n";
      if (nnod < 2)
        return false;

      Tmat.makeCS_X_XZ(X[0],X[1],X[0]+po->directionVector.getValue());
    }
  }

  for (int i = 0; i < 3; i++)
  {
    Tlg[i  ] = Tmat[0][i];
    Tlg[i+3] = Tmat[1][i];
    Tlg[i+6] = Tmat[2][i];
  }

  return true;
}

#ifdef FF_NAMESPACE
} // namespace
#endif
