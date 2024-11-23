// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlBEAM2.H"
#include "FFlLib/FFlFEParts/FFlPMAT.H"
#include "FFlLib/FFlFEParts/FFlPNSM.H"
#include "FFlLib/FFlFEParts/FFlPBEAMSECTION.H"
#include "FFlLib/FFlFEParts/FFlPBEAMECCENT.H"
#include "FFlLib/FFlFEParts/FFlPORIENT.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFaLib/FFaAlgebra/FFaTensor3.H"
#include "FFaLib/FFaAlgebra/FFaMat34.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"


void FFlBEAM2::init()
{
  FFlBEAM2TypeInfoSpec::instance()->setTypeName("BEAM2");
  FFlBEAM2TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::BEAM_ELM);

  ElementFactory::instance()->registerCreator(FFlBEAM2TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlBEAM2::create,int,FFlElementBase*&));

  FFlBEAM2AttributeSpec::instance()->addLegalAttribute("PMAT");
  FFlBEAM2AttributeSpec::instance()->addLegalAttribute("PORIENT", false);
  FFlBEAM2AttributeSpec::instance()->addLegalAttribute("PBEAMORIENT", false);
  FFlBEAM2AttributeSpec::instance()->addLegalAttribute("PBEAMECCENT", false);
  FFlBEAM2AttributeSpec::instance()->addLegalAttribute("PBEAMSECTION");
  FFlBEAM2AttributeSpec::instance()->addLegalAttribute("PBEAMPIN", false);
  FFlBEAM2AttributeSpec::instance()->addLegalAttribute("PNSM", false);
  FFlBEAM2AttributeSpec::instance()->addLegalAttribute("PEFFLENGTH", false);

  FFlBEAM2ElementTopSpec::instance()->setNodeCount(2);
  FFlBEAM2ElementTopSpec::instance()->setNodeDOFs(6);

  FFlBEAM2ElementTopSpec::instance()->myExplEdgePattern = 0xffe4;
}


double FFlBEAM2::getMassDensity() const
{
  FFlPMAT* pmat = dynamic_cast<FFlPMAT*>(this->getAttribute("PMAT"));
  double rho = pmat ? pmat->materialDensity.getValue() : 0.0;

  // Check for non-structural mass
  FFlPNSM* pnsm = dynamic_cast<FFlPNSM*>(this->getAttribute("PNSM"));
  if (!pnsm) return rho;

  // Beam cross section
  FFlPBEAMSECTION* psec = dynamic_cast<FFlPBEAMSECTION*>
			  (this->getAttribute("PBEAMSECTION"));
  if (!psec) return rho;

  double area = psec->crossSectionArea.getValue();
  if (area < 1.0e-16) return rho; // avoid division by 0

  // Modify the mass density to account for the non-structural mass
  return rho + pnsm->NSM.getValue() / area;
}


bool FFlBEAM2::getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const
{
  FFlPBEAMSECTION* psec = dynamic_cast<FFlPBEAMSECTION*>
			  (this->getAttribute("PBEAMSECTION"));
  if (!psec) return false; // Should not happen

  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(2)->getPos());

  // Account for eccentricities
  FFlPBEAMECCENT* pecc = dynamic_cast<FFlPBEAMECCENT*>
			 (this->getAttribute("PBEAMECCENT"));
  if (pecc)
  {
    v1 += pecc->node1Offset.getValue();
    v2 += pecc->node2Offset.getValue();
  }

  double A = psec->crossSectionArea.getValue();
  double L = (v2-v1).length();
  cog      = (v1+v2)*0.5;
  volume   = A*L;

  // Compute local coordinate system
  FaMat33 Telm;
  FaVec3 Zaxis = this->getLocalZdirection();
  if (Zaxis.isZero())
    Telm.makeGlobalizedCS(v2-v1);
  else
  {
    Telm[0] = v2-v1;
    Telm[1] = Telm[0] ^ Zaxis;
    Telm[0].normalize();
    Telm[1].normalize();
    Telm[2] = Telm[0] ^ Telm[1];
  }

  // Set up the inertia tensor in local axes
  double Ixx = psec->Iy.getValue() + psec->Iz.getValue();
  inertia[0] = L*(Ixx > 0.0 ? Ixx : psec->It.getValue());
  inertia[1] = L*(psec->Iy.getValue() + A*L*L/12.0);
  inertia[2] = L*(psec->Iz.getValue() + A*L*L/12.0);
  inertia[3] = inertia[4] = inertia[5] = 0.0;

  // Transform to global axes
  inertia.rotate(Telm.transpose());
  return true;
}


int FFlBEAM2::getNodalCoor(double* X, double* Y, double* Z) const
{
  int ierr = this->FFlElementBase::getNodalCoor(X,Y,Z);
  if (ierr < 0) return ierr;

  // Get beam orientation
  FaVec3 Zaxis = this->getLocalZdirection();
  if (Zaxis.isZero())
  {
    // No beam orientation given, compute a "globalized" Z-axis
    FaMat34 Telm;
    Telm.makeGlobalizedCS(FaVec3(X[0],Y[0],Z[0]),FaVec3(X[1],Y[1],Z[1]));
    Zaxis = Telm[VZ];
    ierr = 1;
    ListUI <<"   * Note: Computing globalized Z-axis for beam element "
           << this->getID() <<" : "<< Zaxis <<"\n";
  }
  X[2] = X[0] + Zaxis[0];
  Y[2] = Y[0] + Zaxis[1];
  Z[2] = Z[0] + Zaxis[2];
  X[3] = X[0];
  Y[3] = Y[0];
  Z[3] = Z[0];
  X[4] = X[1];
  Y[4] = Y[1];
  Z[4] = Z[1];

  // Get beam eccentricities, if any
  FFlPBEAMECCENT* curEcc = dynamic_cast<FFlPBEAMECCENT*>(this->getAttribute("PBEAMECCENT"));
  if (!curEcc) return ierr;

  FaVec3 e1 = curEcc->node1Offset.getValue();
  FaVec3 e2 = curEcc->node2Offset.getValue();
  X[0] += e1[0];
  Y[0] += e1[1];
  Z[0] += e1[2];
  X[1] += e2[0];
  Y[1] += e2[1];
  Z[1] += e2[2];
  X[2] += e1[0];
  Y[2] += e1[1];
  Z[2] += e1[2];

  return ierr;
}


FaVec3 FFlBEAM2::getLocalZdirection() const
{
  FFlPORIENT* pori = dynamic_cast<FFlPORIENT*>(this->getAttribute("PORIENT"));
  if (!pori) // If no PORIENT, try the equivalent old name also
    pori = dynamic_cast<FFlPORIENT*>(this->getAttribute("PBEAMORIENT"));
  return pori ? pori->directionVector.getValue() : FaVec3();
}
