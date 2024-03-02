// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlFEParts/FFlBEAM3.H"
#include "FFlLib/FFlFEParts/FFlPMAT.H"
#include "FFlLib/FFlFEParts/FFlPNSM.H"
#include "FFlLib/FFlFEParts/FFlPBEAMSECTION.H"
#include "FFlLib/FFlFEParts/FFlPBEAMECCENT.H"
#include "FFlLib/FFlFEParts/FFlPBEAMPIN.H"
#include "FFlLib/FFlFEParts/FFlPORIENT.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlLinkHandler.H"
#include "FFaLib/FFaAlgebra/FFaTensor3.H"
#include "FFaLib/FFaAlgebra/FFaMat33.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"


void FFlBEAM3::init()
{
  FFlBEAM3TypeInfoSpec::instance()->setTypeName("BEAM3");
  FFlBEAM3TypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::BEAM_ELM);

  ElementFactory::instance()->registerCreator(FFlBEAM3TypeInfoSpec::instance()->getTypeName(),
					      FFaDynCB2S(FFlBEAM3::create,int,FFlElementBase*&));

  FFlBEAM3AttributeSpec::instance()->addLegalAttribute("PMAT");
  FFlBEAM3AttributeSpec::instance()->addLegalAttribute("PORIENT", false);
  FFlBEAM3AttributeSpec::instance()->addLegalAttribute("PORIENT3", false);
  FFlBEAM3AttributeSpec::instance()->addLegalAttribute("PBEAMECCENT", false);
  FFlBEAM3AttributeSpec::instance()->addLegalAttribute("PBEAMSECTION");
  FFlBEAM3AttributeSpec::instance()->addLegalAttribute("PBEAMPIN", false);
  FFlBEAM3AttributeSpec::instance()->addLegalAttribute("PNSM", false);

  FFlBEAM3ElementTopSpec::instance()->setNodeCount(3);
  FFlBEAM3ElementTopSpec::instance()->setNodeDOFs(6);

  FFlBEAM3ElementTopSpec::instance()->myExplEdgePattern = 0xffe4;
}


bool FFlBEAM3::split(Elements& newElem, FFlLinkHandler* owner, int)
{
  newElem.clear();
  newElem.reserve(2);

  int elm_id = owner->getNewElmID();
  int or3_id = this->getAttributeID("PORIENT3");
  int ec3_id = this->getAttributeID("PBEAMECCENT");
  int pin_id = this->getAttributeID("PBEAMPIN");

  FFlPORIENT*  orient = NULL;
  FFlPORIENT3* orien3 = NULL;
  if (or3_id)
    orien3 = dynamic_cast<FFlPORIENT3*>(owner->getAttribute("PORIENT3",or3_id));

  FFlPBEAMECCENT* ecc = NULL;
  FFlPBEAMECCENT* ec3 = NULL;
  if (ec3_id)
    ec3 = dynamic_cast<FFlPBEAMECCENT*>(owner->getAttribute("PBEAMECCENT",ec3_id));

  FFlPBEAMPIN* pin = NULL;
  if (pin_id)
    pin = dynamic_cast<FFlPBEAMPIN*>(owner->getAttribute("PBEAMPIN",pin_id));

  for (int i = 0; i < 2; i++)
  {
    newElem.push_back(ElementFactory::instance()->create("BEAM2",elm_id++));
    if (orien3)
    {
      FaVec3 univec = orien3->directionVector[i].getValue() + orien3->directionVector[i+1].getValue();
      orient = dynamic_cast<FFlPORIENT*>(AttributeFactory::instance()->create("PORIENT",elm_id-1));
      orient->directionVector.setValue(univec.normalize());
      if (newElem[i]->setAttribute(orient))
        owner->addAttribute(orient);
    }
    if (ec3)
    {
      ecc = dynamic_cast<FFlPBEAMECCENT*>(AttributeFactory::instance()->create("PBEAMECCENT",elm_id-1));
      if (i == 0)
      {
        ecc->node1Offset.setValue(ec3->node1Offset.getValue());
        ecc->node2Offset.setValue(ec3->node2Offset.getValue());
        this->clearAttribute("PBEAMECCENT");
      }
      else
      {
        ecc->node1Offset.setValue(ec3->node2Offset.getValue());
        ecc->node2Offset.setValue(ec3->node3Offset.getValue());
      }
      if (newElem[i]->setAttribute(ecc))
        owner->addAttribute(ecc);
    }
  }

  if (pin)
  {
    if (pin->PA.getValue() == 0)
      newElem[1]->setAttribute(pin);
    else if (pin->PB.getValue() == 0)
      newElem[0]->setAttribute(pin);
    else
      ListUI <<"\n  ** Warning: parabolic beam element "<< this->getID()
             <<" has pin flags at both ends, unsupported (ignored).";
    this->clearAttribute("PBEAMPIN");
  }

#ifdef FFL_DEBUG
  bool silence = false;
#else
  bool silence = true;
#endif

  if (or3_id == this->getID())
    owner->removeAttribute("PORIENT3",or3_id,silence);

  if (ec3_id == this->getID())
    owner->removeAttribute("PBEAMECCENT",ec3_id,silence);

  newElem[0]->setNode(1,this->getNodeID(1));
  newElem[0]->setNode(2,this->getNodeID(2));
  newElem[1]->setNode(1,this->getNodeID(2));
  newElem[1]->setNode(2,this->getNodeID(3));

  return true;
}


double FFlBEAM3::getMassDensity() const
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


bool FFlBEAM3::getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const
{
  FFlPBEAMSECTION* psec = dynamic_cast<FFlPBEAMSECTION*>
			  (this->getAttribute("PBEAMSECTION"));
  if (!psec) return false; // Should not happen

  FaVec3 v1(this->getNode(1)->getPos());
  FaVec3 v2(this->getNode(2)->getPos());
  FaVec3 v3(this->getNode(3)->getPos());

  // Account for eccentricities
  FFlPBEAMECCENT* pecc = dynamic_cast<FFlPBEAMECCENT*>
			 (this->getAttribute("PBEAMECCENT"));
  if (pecc)
  {
    v1 += pecc->node1Offset.getValue();
    v2 += pecc->node2Offset.getValue();
    v3 += pecc->node3Offset.getValue();
  }

  double A = psec->crossSectionArea.getValue();
  double L = (v2-v1).length() + (v3-v2).length();
  cog      = (v1+v3)*0.25 + v2*0.5;
  volume   = A*L;

  // Compute local coordinate system
  FFlPORIENT* pori = dynamic_cast<FFlPORIENT*>(this->getAttribute("PORIENT"));
  FFlPORIENT3* po3 = dynamic_cast<FFlPORIENT3*>(this->getAttribute("PORIENT3"));

  FaMat33 Telm;
  if (pori || po3)
  {
    const FaVec3& n = pori ? pori->directionVector.getValue() : po3->directionVector[1].getValue();
    Telm[0] = v3-v1;
    Telm[1] = Telm[0] ^ n;
    Telm[0].normalize();
    Telm[1].normalize();
    Telm[2] = Telm[0] ^ Telm[1];
  }
  else
    Telm.makeGlobalizedCS(v3-v1);

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
