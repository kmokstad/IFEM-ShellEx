// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlIOAdaptors/FFlNastranReader.H"
#include "FFlLib/FFlLinkHandler.H"
#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlFEParts/FFlPORIENT.H"
#include "FFlLib/FFlFEParts/FFlPBEAMSECTION.H"
#include "FFlLib/FFlFEParts/FFlPBEAMECCENT.H"
#include "FFlLib/FFlFEParts/FFlPBUSHECCENT.H"
#include "FFlLib/FFlFEParts/FFlPCOORDSYS.H"
#include "FFlLib/FFlFEParts/FFlPSPRING.H"
#include "FFlLib/FFlFEParts/FFlPMASS.H"
#include "FFlLib/FFlFEParts/FFlPMAT.H"
#include "FFlLib/FFlFEParts/FFlPMATSHELL.H"
#include "FFlLib/FFlFEParts/FFlPEFFLENGTH.H"
#include "FFlLib/FFlFEParts/FFlCLOAD.H"

#include "FFaLib/FFaAlgebra/FFaAlgebra.H"
#include "FFaLib/FFaAlgebra/FFaMat33.H"
#include "FFaLib/FFaAlgebra/FFaMath.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"

#include <algorithm>
#if FFL_DEBUG > 1
#include <iomanip>

static void printMatrix6 (const double A[6][6])
{
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
      std::cout << std::setw(13) << A[i][j];
    std::cout << std::endl;
  }
}
#endif

#define CREATE_ATTRIBUTE(type,name,id) \
  dynamic_cast<type*>(AttributeFactory::instance()->create(name,id))

#define GET_ATTRIBUTE(type,name,id) \
  dynamic_cast<type*>(myLink->getAttribute(name,id))

////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::resolveCoordinates ()
{
  // First, assign the default coordinate system to all grid points
  // that have not been assigned any coordinate system through the GRID entry

  const int defaultCP = gridDefault ? (gridDefault->CP ? gridDefault->CP:0):0;
  const int defaultCD = gridDefault ? (gridDefault->CD ? gridDefault->CD:0):0;
  for (NodesCIter n = myLink->nodesBegin(); n != myLink->nodesEnd(); ++n)
  {
    int ID = (*n)->getID();
    IntMap::iterator cpit = nodeCPID.find(ID);
    if (cpit == nodeCPID.end())
    {
      if (defaultCP > 0)
      {
#if FFL_DEBUG > 1
	std::cout <<"Node "<< ID <<" assigned default coordinate system "
		  << defaultCP << std::endl;
#endif
        nodeCPID[ID] = defaultCP;
      }
    }
    else if (cpit->second < 1)
      nodeCPID.erase(cpit);

    // Displacement coordinate systems (and beam orientation and eccentricities)
    IntMap::iterator cdit = nodeCDID.find(ID);
    if (cdit == nodeCDID.end())
    {
      if (defaultCD > 0)
      {
#if FFL_DEBUG > 1
	std::cout <<"Node "<< ID
		  <<" assigned default solution coordinate system "
		  << defaultCD << std::endl;
#endif
        nodeCDID[ID] = defaultCD;
      }
    }
    else if (cdit->second < 1)
      nodeCDID.erase(cdit);

    // Check external status
    IntMap::iterator xit = nodeStat.find(ID);
    if (xit != nodeStat.end())
    {
      if (xit->second > 0)
      {
#if FFL_DEBUG > 1
	std::cout <<"External Node "<< ID << std::endl;
#endif
        (*n)->setExternal();
      }
      else if (xit->second < 0)
      {
#if FFL_DEBUG > 1
	std::cout <<"Constrained Node "<< ID <<": "<< xit->second << std::endl;
#endif
        (*n)->setStatus(xit->second);
      }
      else
      {
        nWarnings++;
        ListUI <<"\n  ** Warning: Node "<< ID <<" has been specified on an ASET"
               <<" entry with a zero component number code (ignored).\n";
      }
    }
  }

  // Then, transform all nodes to the common "global" coordinate system

  bool ok = true;
  while (!nodeCPID.empty())
  {
    IntMap::iterator nit = nodeCPID.begin();
    ok &= this->transformNode(nit);
#if FFL_DEBUG > 1
    std::cout << std::endl;
#endif
  }

  // Now that all nodes have their global coordinates, round to 10 significant
  // digits to avoid checksum issues when saving and reopening the model

  for (NodesCIter n = myLink->nodesBegin(); n != myLink->nodesEnd(); ++n)
    (*n)->position().round(10);

  // Add all local coordinate systems which are not referred by any elements
  // to the link object. They may be used to aid the mechanism modeling later.

  for (const std::pair<const int,CORD*>& ics : cordSys)
    if (!myLink->getAttribute("PCOORDSYS",ics.first))
    {
      if (this->computeTmatrix(ics.first,*ics.second))
      {
        FFlPCOORDSYS* cs = CREATE_ATTRIBUTE(FFlPCOORDSYS,"PCOORDSYS",ics.first);
        cs->Origo = ics.second->Origo.round(10);
        cs->Zaxis = ics.second->Zaxis.round(10);
        cs->XZpnt = ics.second->XZpnt.round(10);
        myLink->addAttribute(cs);
      }
      else
        ok = false;
    }

  return ok;
}

////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::transformNode (IntMap::iterator& nit)
{
  int CP = nit->second;

#if FFL_DEBUG > 1
  std::cout <<" Transforming node "<< nit->first <<", CP = "<< CP;
#endif

  bool ok = this->transformPoint(myLink->getNode(nit->first)->position(),CP);
  nodeCPID.erase(nit);
  return ok;
}

////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::transformPoint (FaVec3& X, const int CID,
                                       bool orientationOnly)
{
  std::map<int,CORD*>::iterator cit = cordSys.find(CID);
  if (cit == cordSys.end())
  {
    ListUI <<"\n *** Error: Coordinate system "<< CID <<" does not exist.\n";
    return false;
  }

  // Compute the 3x4 transformation matrix for CID if not already computed
  if (!this->computeTmatrix(CID,*(cit->second)))
    return false;

#if FFL_DEBUG > 1
  std::cout <<" "<< X;
#endif

  if (cit->second->type == cylindrical)
  {
    // Transform from cylindrical to local cartesian coordinates
    const double PIo180 = M_PI/180.0;
    const double r      = X[0];
    const double theta  = X[1];
    X[0] = r*cos(PIo180*theta);
    X[1] = r*sin(PIo180*theta);
  }
  else if (cit->second->type == spherical)
  {
    // Transform from spherical to local cartesian coordinates
    const double PIo180 = M_PI/180.0;
    const double r      = X[0];
    const double theta  = X[1];
    const double phi    = X[2];
    X[0] = r*cos(PIo180*theta)*sin(PIo180*phi);
    X[1] = r*sin(PIo180*theta)*sin(PIo180*phi);
    X[2] = r*cos(PIo180*phi);
  }

  // Transform from local cartesian to global cartesian coordinates
  if (orientationOnly)
    X = cit->second->Tmat.direction() * X;
  else
    X = cit->second->Tmat * X;

#if FFL_DEBUG > 1
  std::cout <<" --> "<< X;
#endif

  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::computeTmatrix (const int CID, CORD& cs)
{
  if (cs.isComputed) return true;

  bool ok = true;
  if (cs.G[0] > 0)

    // This coordinate system is defined through three nodes
    for (int i = 0; i < 3; i++)
    {
      int GID = cs.G[i];
      IntMap::iterator nit = nodeCPID.find(GID);
      if (nit != nodeCPID.end())
      {
        // This node is again defined in (hopefully!) another coordinate system
        if (nit->second == CID)
        {
          ok = false;
          nodeCPID.erase(nit);
          ListUI <<"\n *** Error: Coordinate system "<< CID
                 <<" is defined through nodes whose coordinates"
                 <<" are given in the same coordinate system.\n";
        }
        else
          ok &= this->transformNode(nit);
      }

      // Get the (now global) coordinates of the nodal point
      switch (i)
      {
        case 0 : cs.Origo = myLink->getNode(GID)->getPos(); break;
        case 1 : cs.Zaxis = myLink->getNode(GID)->getPos(); break;
        case 2 : cs.XZpnt = myLink->getNode(GID)->getPos();
      }
    }

  else if (cs.RID == CID)
  {
    ok = false;
    ListUI <<"\n *** Error: Coordinate system "<< CID <<" is defined by points"
           <<" that are given in the same coordinate system.\n";
  }
  else if (cs.RID > 0)
    // This coordinate system is defined through three spatial points
    // given in the coordinate system RID (which must be different than CID)
    ok &= (this->transformPoint(cs.Origo,cs.RID) &&
	   this->transformPoint(cs.Zaxis,cs.RID) &&
	   this->transformPoint(cs.XZpnt,cs.RID));

  // Compute the 3x4 transformation matrix
  if (ok)
  {
    cs.Tmat.makeCS_Z_XZ(cs.Origo,cs.Zaxis,cs.XZpnt);
#ifdef FFL_DEBUG
    std::cout <<"\nCoordinate system "<< CID <<", type = "<< cs.type
	      <<", Tmat ="<< cs.Tmat << std::endl;
#endif
  }
  cs.isComputed = true;
  return ok;
}

////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::getTmatrixAtPt (const int CID, CORD& cs,
                                       const FaVec3& X, FaMat33& T)
{
  // Compute the 3x4 transformation matrix for CID if not already computed
  if (!this->computeTmatrix(CID,cs))
    return false;

  const double Zero = 1.0e-16;
  if (cs.type == cylindrical)
  {
    // Rotate the local axes such that global point X lies in the local XZ-plane
    FaVec3 Xloc = cs.Tmat.inverse() * X;
    double lXY  = hypot(Xloc[0],Xloc[1]);
    if (lXY > Zero)
      T = cs.Tmat.direction() * FaMat33::makeZrotation(acos(Xloc[0]/lXY));
    else
      T = cs.Tmat.direction();
  }
  else if (cs.type == spherical)
  {
    // Rotate the local axes such that global point X lies on the local X-axis
    FaVec3 Xloc = cs.Tmat.inverse() * X;
    double lXY  = hypot(Xloc[0],Xloc[1]);
    if (lXY > Zero)
      T = cs.Tmat.direction() * FaMat33::makeZrotation(acos(Xloc[0]/lXY));
    else
      T = cs.Tmat.direction();

    double lXYZ = hypot(lXY,Xloc[2]);
    if (lXYZ > Zero) T = T * FaMat33::makeYrotation(acos(Xloc[2]/lXYZ));
    T.shift(1);
  }
  else
    T = cs.Tmat.direction();

#if FFL_DEBUG > 1
  std::cout <<"\nTransformation matrix at global point "<< X <<" :\n"<< T;
#endif
  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::resolveAttributes ()
{
  FFlElementBase* curElm;
  std::string     curTyp;
  int             EID, PID, noProp = 0, noMat = 0;
  std::map<int,bool> invalidMat;

  bool ok = true;
  nWarnings += myLink->sortElementsAndNodes(true);
  int newEID = myLink->getNewElmID(); // Used for the auto-created elements

  // Note: We cannot use an iterator in this loop,
  // since elements may be added within the loop (in resolveWeldElement)
  for (int iel = 0; iel < myLink->getElementCount(); iel++)
  {
    // Get the property identifier for this element, they are temporarily stored
    // in different data structures depending on the element type

    curElm = myLink->getElement(iel,true);
    curTyp = curElm->getTypeName();
#if FFL_DEBUG > 2
    std::cout <<"\t"<< curTyp <<"\t"<< curElm->getID() << std::endl;
#endif

    if (curElm->getCathegory() == FFlTypeInfoSpec::SHELL_ELM)
    {
      // Shell elements
      IntMap::iterator pit = shellPID.find(curElm->getID());
      PID = pit == shellPID.end() ? 0 : pit->second;

      // Check for non-structural mass, add attribute if present
      if (shellPIDnsm.find(PID) != shellPIDnsm.end())
        ok &= curElm->setAttribute("PNSM",PID);

      // Check for shell thickness or composite properties
      if (PTHICKs.find(PID) != PTHICKs.end())
        ok &= curElm->setAttribute("PTHICK",PID);
      else if (PCOMPs.find(PID) != PCOMPs.end())
      {
        ok &= curElm->setAttribute("PCOMP",PID);
        PID = 0; // PCOMP contains the material information
      }
    }
    else if (curElm->getCathegory() == FFlTypeInfoSpec::BEAM_ELM)
    {
      PID = curElm->getAttributeID("PBEAMSECTION");
      if (std::find(myWelds.begin(),myWelds.end(),curElm) != myWelds.end())
        // This is a weld connection element, must resolve the orientation
        // vector and the associated constraint elements, if any
        ok &= this->resolveWeldElement(curElm,newEID,PID);

      else
      {
        // For beam elements, transform the eccentricity (if any)
        // and orientation vectors to the global coordinate system
        // For rod elements (PID > 0), nothing needs to be done here
        if (PID == 0) PID = this->resolveBeamAttributes(curElm,ok);

        // Check for non-structural mass, add attribute if present
        if (beamPIDnsm.find(PID) != beamPIDnsm.end())
          ok &= curElm->setAttribute("PNSM",PID);
      }
    }
    else if (curTyp == "SPRING" || curTyp == "RSPRING")
    {
      EID = curElm->getID();
      PID = 0; // No geometric properties for spring elements
      // Compute the stiffness matrix in the global coordinate system
      double K = 0.0;
      std::map<int,double>::iterator Kit = sprK.find(EID);
      if (Kit != sprK.end())
        K = Kit->second;
      else
      {
        IntMap::iterator pit = sprPID.find(EID);
        Kit = propK.find(pit == sprPID.end() ? EID : pit->second);
        if (Kit != propK.end()) K = Kit->second;
      }
      IntMap::iterator cit = sprComp.find(EID);
      int Comp = cit == sprComp.end() ? 0 : cit->second;
      this->resolveSpringAttributes(curElm,K,Comp,ok);
    }
    else if (curTyp == "BUSH")
    {
      EID = curElm->getID();
      PID = 0; // No geometric properties for bushing elements
      // Transform eccentricity and orientation vectors,
      // if any, to the global coordinate system
      IntMap::iterator cit = sprComp.find(EID);
      std::map<int,double>::iterator Sit = sprK.find(EID);
      int CID  = cit == sprComp.end() ? 0 : cit->second;
      double S = Sit == sprK.end() ? 0.5 : Sit->second;
      this->resolveBushAttributes(curElm,S,CID,ok);
    }
    else if (curTyp == "CMASS")
    {
      PID = 0; // No geometric properties for mass elements
      // Transform the specified mass matrix to the global coordinate system
      IntMap::iterator cit = massCID.find(curElm->getID());
      std::map<int,FaVec3*>::iterator xit = massX.find(curElm->getID());
      int CID = cit != massCID.end() ? cit->second : 0;
      FaVec3* X = xit != massX.end() ? xit->second : NULL;
      if (CID || X)
        ok &= this->transformMassMatrix(curElm,CID,X);
    }
    else if (curTyp == "RGD" || curTyp == "RBAR" || curTyp == "WAVGM")

      PID = 0; // No geometric properties for rigid and constraint elements

    else
    {
      // Solid elements
      IntMap::iterator pit = solidPID.find(curElm->getID());
      PID = pit == solidPID.end() ? 0 : pit->second;
    }

    if (PID > 0)
    {
      // Find the material identifier for this element
      IntMap::iterator mit = propMID[curElm->getCathegory()].find(PID);
      if (mit != propMID[curElm->getCathegory()].end())
      {
        FFlPMAT* theMat;
        const int MID = mit->second;
        if (myLink->getAttribute("PMATSHELL",MID))
        {
          if (!curElm->setAttribute("PMATSHELL",MID))
          {
            ok = false;
            if (++noMat < 11)
              ListUI <<"\n *** Error: Material PMATSHELL "<< MID
                     <<" referenced by \""<< curTyp
                     <<"\" element "<< curElm->getID()
                     <<" through property "<< PID
                     <<" is an illegal attribute.\n";
          }
        }
        else if ((theMat = GET_ATTRIBUTE(FFlPMAT,"PMAT",MID)))
        {
          if (!curElm->setAttribute("PMAT",MID))
            ok = false;
          else if (curTyp == "BEAM2")
          {
            if (theMat->shearModule.getValue() <= 0.0)
            {
              // Reset the zero shear-modulus for beam elements
              double E = theMat->youngsModule.getValue();
              double v = theMat->poissonsRatio.getValue();
              double G = v > -1.0 ? E/(2.0+v+v) : -1.0;
              if (G > 0.0)
              {
                nNotes++;
                ListUI <<"\n   * Note: Material "<< theMat->getID()
                       <<" has a zero shear modulus (G) but is used by beams.\n"
                       <<"           Resetting to "<< G <<" = E/(2+2*nu).\n";
                theMat->shearModule = round(G,10); // use 10 significant digits
              }
              else if (invalidMat.find(theMat->getID()) == invalidMat.end())
              {
                nWarnings++;
                ListUI <<"\n  ** Warning: Material "<< theMat->getID()
                       <<" has invalid parameters, E = "<< E <<", nu = "<< v
                       <<", G = "<< theMat->shearModule.getValue() <<".\n";
                invalidMat[theMat->getID()] = true;
              }
            }
          }
          else if (theMat->poissonsRatio.getValue() < 0.0 ||
                   theMat->poissonsRatio.getValue() > 0.5)
          {
            // Reset the invalid Poissons ratio for shell and solid elements
            double E = theMat->youngsModule.getValue();
            double G = theMat->shearModule.getValue();
            double v = G > 0.0 ? 0.5*E/G - 1.0 : -1.0;
            if (v >= 0.0 && v < 0.5)
            {
              nNotes++;
              ListUI <<"\n   * Note: Material "<< theMat->getID()
                     <<" has an invalid Poisson's ratio "
                     << theMat->poissonsRatio.getValue()
                     <<", but is used by shell- and/or solid element.\n"
                     <<"           Resetting to "<< v <<" = 0.5*E/G - 1.\n";
              theMat->poissonsRatio = round(v,10); // use 10 significant digits
            }
            else if (invalidMat.find(theMat->getID()) == invalidMat.end())
            {
              nWarnings++;
              ListUI <<"\n  ** Warning: Material "<< theMat->getID()
                     <<" has invalid parameters, E = "<< E <<", nu = "
                     << theMat->poissonsRatio.getValue()
                     <<", G = "<< G <<".\n";
              invalidMat[theMat->getID()] = true;
            }
          }
        }
        else if (++noMat < 11)
        {
          ok = false;
          ListUI <<"\n *** Error: Material "<< MID <<" referenced by \""
                 << curTyp <<"\" element "<< curElm->getID()
                 <<" through property "<< PID <<" does not exist.\n";
        }
      }
      else if (++noProp < 11)
      {
        ok = false;
        ListUI <<"\n *** Error: Property "<< PID <<" referenced by \""
               << curTyp <<"\" element "<< curElm->getID()
               <<" does not exist.\n";
      }
    }
  }

  if (noProp > 10)
    ListUI <<"\n *** Error: Non-existing property detected for "
           << noProp <<" elements.\n";

  if (noMat > 10)
    ListUI <<"\n *** Error: Non-existing material detected for "
           << noMat <<" elements.\n";

  propMID.clear();
  massCID.clear();
  solidPID.clear();
  shellPID.clear();
  weldGS.clear();
  weld.clear();
  return ok;
}

////////////////////////////////////////////////////////////////////////////////

int FFlNastranReader::resolveBeamAttributes (FFlElementBase* curElm, bool& ok)
{
  FaVec3* eVec[2] = { NULL, NULL };

  int EID = curElm->getID();
  int Ecc = curElm->getAttributeID("PBEAMECCENT");
  if (Ecc > 0)
  {
    // This beam has eccentricities.
    // Check if they need to be transformed to the global coordinate system

    FFlPBEAMECCENT* theEc = GET_ATTRIBUTE(FFlPBEAMECCENT,"PBEAMECCENT",Ecc);
    if (!theEc)
    {
      ok = false;
      std::cerr <<"FFlNastranReader::resolveBeamAttributes: Internal error, "
		<<"EID="<< EID <<" PID="<< Ecc <<", missing eccentricty."
		<< std::endl;
      return 0;
    }

    eVec[0] = &theEc->node1Offset.data();
    eVec[1] = &theEc->node2Offset.data();

    for (int i = 0; i < 2; i++)
    {
      int ID = curElm->getNodeID(i+1);
      IntMap::iterator cit = nodeCDID.find(ID);
      if (cit != nodeCDID.end())
      {
        // This element node has a local coordinate system for displacements
        // in which also the eccentricity vector is defined
#if FFL_DEBUG > 1
        std::cout <<" Transforming eccentricity vector "<< i+1
                  <<" for beam element "<< EID <<", CD = "<< cit->second;
#endif
        ok &= this->transformPoint(*eVec[i],cit->second,true);
#if FFL_DEBUG > 1
        std::cout << std::endl;
#endif
        eVec[i]->round(10); // round to 10 significant digits
      }
    }
  }

  // Resolve property ID and beam orientation
  std::map<int,BEAMOR*>::iterator bit = bOri.find(EID);
  if (bit == bOri.end())
  {
    ok = false;
    std::cerr <<"FFlNastranReader::resolveBeamAttributes: Internal error, "
	      <<"EID="<< EID <<", missing orientation vector."<< std::endl;
    return 0;
  }

  BEAMOR& curOr = *(bit->second);
  BEAMOR* defOr = curOr.isBAR ? barDefault : beamDefault;
  if (defOr)
  {
    // A default entry was specified
    // fill empty entries with default values for this element
    if (curOr.empty[0] && !defOr->empty[0])
    {
      curOr.PID = defOr->PID;
      curOr.empty[0] = false;
    }
    if (curOr.empty[1] && !defOr->empty[1])
    {
      curOr.X[0] = defOr->X[0];
      curOr.empty[1] = false;
    }
    if (curOr.empty[2] && !defOr->empty[2])
    {
      curOr.X[1] = defOr->X[1];
      curOr.empty[2] = false;
    }
    if (curOr.empty[3] && !defOr->empty[3])
    {
      curOr.X[2] = defOr->X[3];
      curOr.empty[3] = false;
    }
    if (curOr.empty[4] && !defOr->empty[4])
    {
      curOr.G0 = defOr->G0;
      curOr.empty[4] = false;
    }
  }

  int PID = curOr.PID;
  ok &= curElm->setAttribute("PBEAMSECTION",PID);

  FaVec3 x, y;
  if (curOr.G0 > 0)
    // Orientation (local y-axis) is given by vector from local node 1 and G0
    ok &= this->getElementAxis(curElm,1,-curOr.G0,y);

  else
  {
    // Orientation (local y-axis) is given component-wise
    IntMap::iterator cit = nodeCDID.find(curElm->getNodeID(1));
    if (cit == nodeCDID.end())
      y = curOr.X;
    else
    {
      // This element node has a local coordinate system for displacements
      // in which also the orientation vector is defined
#if FFL_DEBUG > 1
      std::cout <<" Transforming orientation vector for beam element "<< EID
		<<", CD = "<< cit->second;
#endif
      if (this->transformPoint(curOr.X,cit->second,true))
        y = curOr.X;
      else
        ok = false;
#if FFL_DEBUG > 1
      std::cout << std::endl;
#endif
    }
  }

  // Now find the local Z-axis of the beam element
  if (!this->getElementAxis(curElm,1,2,x))
    ok = false;
  else
  {
    FFlPORIENT* theOr = CREATE_ATTRIBUTE(FFlPORIENT,"PORIENT",EID);
    theOr->directionVector = (x^y).truncate().round(10);
    int newPID = myLink->addUniqueAttribute(theOr);
#if FFL_DEBUG > 1
    if (newPID == EID)
      std::cout <<"Attribute PORIENT, ID = "<< newPID
                <<", Fields: "<< theOr->directionVector << std::endl;
#endif
    ok &= curElm->setAttribute("PORIENT",newPID);
  }

  // Special resolving of beams with neutral axis offset
  FFlPBEAMSECTION* theSec = GET_ATTRIBUTE(FFlPBEAMSECTION,"PBEAMSECTION",PID);
  if (!theSec) return PID; // No neutral axis offset for this beam

  double Sy = theSec->Sy.getValue();
  double Sz = theSec->Sz.getValue();
  if (fabs(Sy)+fabs(Sz) <= 1.0e-16)
    return PID; // Negligible neutral axis offset for this beam

  if (Ecc == 0)
  {
    Ecc = EID; // This beam element does not have eccentricities yet, create one
    FFlPBEAMECCENT* theEc = CREATE_ATTRIBUTE(FFlPBEAMECCENT,"PBEAMECCENT",Ecc);
    eVec[0] = &theEc->node1Offset.data();
    eVec[1] = &theEc->node2Offset.data();
    myLink->addAttribute(theEc);
    ok &= curElm->setAttribute("PBEAMECCENT",Ecc);
  }

  if (fabs(theSec->phi.getValue()) > 1.0e-6)
  {
    // Transform neutral axis offset to local element axes (from principal axes)
    double fi = theSec->phi.getValue()*M_PI/180.0;
    double cf = cos(fi);
    double sf = sin(fi);
    double S1 = cf*Sy - sf*Sz;
    double S2 = cf*Sz + sf*Sy;
    Sy = S1;
    Sz = S2;
  }

  // Update the eccentricity vectors with neutral axis offset
  FaVec3 offset = (Sy*y.normalize() + Sz*(x^y).normalize());
#if FFL_DEBUG > 1
  std::cout <<"Attribute PBEAMECCENT, ID = "<< Ecc <<", offset = "<< offset;
  for (int i = 0; i < 2; i++)
  {
    std::cout <<"\n                       E"<< i+1 <<" = "<< *eVec[i];
    (*eVec[i] -= offset).round(10);
    std::cout <<" --> "<< *eVec[i];
  }
  std::cout << std::endl;
#else
  for (int i = 0; i < 2; i++) (*eVec[i] -= offset).round(10);
#endif

  return PID;
}

////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::resolveWeldElement (FFlElementBase* curElm,
					   int& newEID, const int PID)
{
  int EID = curElm->getID();
  IntMap::iterator git = weldGS.find(EID);
  if (git != weldGS.end())
  {
    // A grid point to be projected onto one or two surface patches is given
    int SHID;
    bool ok = true;
    FFlNode* nS = myLink->getNode(git->second);
    if (!nS)
    {
      ok = false;
      ListUI <<"\n *** Error: Non-existing node "<< git->second;
    }
    else
      for (size_t i = 1; i <= 2; i++)
        if ((SHID = curElm->getNodeID(i)) <= 0)
        {
          // Element node "i" must be created as the projection of node nS
          // onto a given surface patch
          std::vector<int> G(1,myLink->getNewNodeID());
          if (SHID < 0)
          {
            // The patch to project node nS onto is defined by element -SHID
            // The referred element must be a shell element
            FFlElementBase* refElm = myLink->getElement(-SHID);
            if (!refElm)
            {
              ok = false;
              ListUI <<"\n *** Error: Non-existing element "<< -SHID;
            }
            else if (refElm->getCathegory() != FFlTypeInfoSpec::SHELL_ELM)
            {
              ok = false;
              ListUI <<"\n *** Error: Element "<< -SHID <<" is not a shell";
            }
            else
            {
              // Create a property-less WAVGM element where the new element
              // node is the reference node, and all element nodes of the
              // referred shell element are the independent nodes
              int nelnod = refElm->getNodeCount();
              for (int k = 1; k <= nelnod; k++)
                G.push_back(refElm->getNodeID(k));
#if FFL_DEBUG > 2
              std::cout <<"WAVGM element "<< newEID <<", Nodes:";
              for (int node : G) std::cout <<" "<< node;
              std::cout << std::endl;
#endif
              FFlElementBase* theElm = ElementFactory::instance()->create("WAVGM",newEID++);
              if (!theElm)
              {
                ListUI <<"\n *** Error: Failure creating WAVGM element "<< --newEID <<".\n";
                return false;
              }

              theElm->setNodes(G);
              if (!myLink->addElement(theElm))
                return false;
            }
          }
          else if (weld.size() < 2)
          {
	    std::cerr <<"FFlNastranReader::resolveWeldElement: Internal error."
		      << std::endl;
            return false;
          }
          else
          {
            // The patch to project node nS onto is explicitly defined by a
            // series of nodes for which a WAVGM element already has been made
            std::map<int,FFlElementBase*>::iterator wit = weld[i-1].find(EID);
            if (wit == weld[i-1].end())
            {
	      std::cerr <<"FFlNastranReader::resolveWeldElement: Internal error."
			<< std::endl;
              return false;
            }
            wit->second->setID(newEID++);
            wit->second->setNode(1,G.front());
            int nelnod = wit->second->getNodeCount();
            for (int k = 2; k <= nelnod && k <= 5; k++)
              G.push_back(wit->second->getNodeID(k));
          }

          // Get the global position of all patch vertices
          std::vector<FaVec3> vx;
          FFlNode* pchNode;
          for (size_t k = 1; k < G.size(); k++)
            if ((pchNode = myLink->getNode(G[k])))
              vx.push_back(pchNode->getPos());
            else
            {
              ok = false;
              ListUI <<"\n *** Error: Non-existing node "<< G[k];
            }

          // Compute the globalized coordinate system of this patch
          // where the surface normal defines the local Z-axis
          FaMat34 Tpch;
          if (vx.size() == 3)
            Tpch.makeGlobalizedCS(vx[0],vx[1],vx[2]);
          else if (vx.size() == 4)
            Tpch.makeGlobalizedCS(vx[0],vx[1],vx[2],vx[3]);
          else
            ok = false;

          if (ok)
          {
            // Now define the new element node of the weld element as the
            // projection of node nS onto the XY-plane of the surface patch
            FaVec3 X = Tpch.projectOnXY(nS->getPos());
            if (myLink->addNode(new FFlNode(G.front(),X.round(10))))
              curElm->setNode(i,G.front());
            else
              return false;
          }
        }

    if (!ok)
    {
      ListUI <<"\n            referred by Weld element "<< EID <<".\n";
      return false;
    }
  }

  // Find the actual length of the beam element and its cross section diameter
  FaVec3 Xaxis;
  if (!this->getElementAxis(curElm,1,2,Xaxis))
    return false;

  FFlPBEAMSECTION* theSec = GET_ATTRIBUTE(FFlPBEAMSECTION,"PBEAMSECTION",PID);
  if (!theSec)
  {
    ListUI <<"\n *** Error: Non-existing beam property "<< PID
           <<" referred by CWELD element "<< EID <<".\n";
    return false;
  }

  double rLength = Xaxis.length();
  double diameter = 2.0*sqrt(theSec->crossSectionArea.getValue()/M_PI);

  // Check if the beam element is too short or too long
  // such that we need to assign an effective length
  if (rLength < 0.2*diameter)
    rLength = 0.2*diameter;
  else if (rLength > 5.0*diameter)
    rLength = 5.0*diameter;
  else
    return true; // the actual length is OK

#if FFL_DEBUG > 1
    std::cout <<"Effective length for WELD element "<< EID
	      <<" changed from "<< Xaxis.length()
	      <<" to "<< rLength << std::endl;
#endif

  // Create and assign an effective length for this beam element
  FFlPEFFLENGTH* theEff = CREATE_ATTRIBUTE(FFlPEFFLENGTH,"PEFFLENGTH",EID);
  theEff->length = round(rLength,10); // round to 10 significant digits
  int newPID = myLink->addUniqueAttribute(theEff);
  return curElm->setAttribute("PEFFLENGTH",newPID);
}

////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::getElementAxis (FFlElementBase* curElm,
				       const int n1, const int n2,
				       FaVec3& axis)
{
  // Find the vector defined by the global nodes GA and GB
  const int GA = n1 > 0 ? curElm->getNodeID(n1) : -n1;
  const int GB = n2 > 0 ? curElm->getNodeID(n2) : -n2;

  // Here we assumed all nodes already have been transformed to global
  // coordinates, i.e. resolveCoordinates must have been called.
  const FFlNode* nA = myLink->getNode(GA);
  const FFlNode* nB = myLink->getNode(GB);
  if (nA && nB)
  {
    axis = nB->getPos() - nA->getPos();
#if FFL_DEBUG > 2
    std::cout <<"Local axis "<< GA <<"-"<< GB <<" for "<< curElm->getTypeName()
              <<" element "<< curElm->getID() <<": "<< axis << std::endl;
#endif
    return true;
  }
  else if (nB)
    ListUI <<"\n *** Error: Non-existing node "<< GA;
  else if (nA)
    ListUI <<"\n *** Error: Non-existing node "<< GB;
  else
    ListUI <<"\n *** Error: Non-existing nodes "<< GA <<", "<< GB;
  ListUI <<" referred by "<< curElm->getTypeName()
         <<" element "<< curElm->getID() <<".\n";
  return false;
}

////////////////////////////////////////////////////////////////////////////////

void FFlNastranReader::resolveSpringAttributes (FFlElementBase* curElm,
						const double S, const int C,
						bool& ok)
{
  // Set up the local stiffness matrix for the spring
  int i, j, k;
  double K[6][6];
  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++)
      K[i][j] = 0.0;

  bool isTranslation = false;
  bool isRotation = false;
  for (i = C; i > 0; i /= 10)
  {
    j = i%10-1;
    if (j >= 0 && j < 3)
    {
      if (isRotation)
      {
        nWarnings++;
        ListUI <<"\n  ** Warning: Invalid component code "<< C
               <<" for spring element "<< curElm->getID() <<".\n             "
               <<" Only the rotational DOFs will be included.\n";
      }
      else
      {
        isTranslation = true;
        K[  j][  j] =  S;
        K[3+j][3+j] =  S;
        K[  j][3+j] = -S;
        K[3+j][  j] = -S;
      }
    }
    else if (j >= 3 && j < 6)
    {
      if (isTranslation)
      {
        nWarnings++;
        ListUI <<"\n  ** Warning: Invalid component code "<< C
               <<" for spring element "<< curElm->getID() <<".\n             "
               <<" Only the translational DOFs will be included.\n";
      }
      else
      {
        isRotation = true;
        K[j-3][j-3] =  S;
        K[j  ][j  ] =  S;
        K[j-3][j  ] = -S;
        K[j  ][j-3] = -S;
      }
    }
  }

#if FFL_DEBUG > 1
  std::cout <<"\nStiffness matrix for spring element "<< curElm->getID() <<"\n";
  printMatrix6(K);
#endif

  // Define the spring type (either translational or rotational)
  FFlPSPRING* theSpr = GET_ATTRIBUTE(FFlPSPRING,"PSPRING",curElm->getID());
  if (isTranslation)
    theSpr->type = FFlPSPRING::TRANS_SPRING;
  else if (isRotation)
    theSpr->type = FFlPSPRING::ROT_SPRING;
  else
  {
    ok = false;
    ListUI <<"\n *** Error: Invalid component code "<< C
           <<" for spring element "<< curElm->getID() <<".\n";
  }

  // Transform the stiffness matrix to global coordinate system
  for (i = 1; i <= 2; i++)
  {
    IntMap::iterator cit = nodeCDID.find(curElm->getNodeID(i));
    if (cit != nodeCDID.end())
    {
      FFlNode* node = myLink->getNode(curElm->getNodeID(i));
      if (!node)
      {
        ok = false;
        ListUI <<"\n *** Error: Non-existing node "<< curElm->getNodeID(i)
               <<" referred by spring element "<< curElm->getID() <<".\n";
      }
      else
        ok &= this->transformSymmMatrix6(K,node->getPos(),cit->second,i);
    }
  }

  if (!ok) return;

  for (i = k = 0; i < 6; i++)
    for (j = 0; j <= i; j++)
      theSpr->K[k++] = round(K[i][j],10); // round to 10 significant digits

#if FFL_DEBUG > 1
  theSpr->print();
#endif
}

////////////////////////////////////////////////////////////////////////////////

void FFlNastranReader::resolveBushAttributes (FFlElementBase* curElm,
					      const double S, const int CID,
					      bool& ok)
{
  int EID = curElm->getID();
  int Ecc = curElm->getAttributeID("PBUSHECCENT");
  if (Ecc > 0)
  {
    // This bushing element has an explicit eccentricity vector.
    // Check if it needs to be transformed to the global coordinate system
    IntMap::iterator cit = sprPID.find(EID);
    if (cit != sprPID.end())
    {
      // The eccentricity vector is defined in a local coordinate system
      FFlPBUSHECCENT* theEc = GET_ATTRIBUTE(FFlPBUSHECCENT,"PBUSHECCENT",Ecc);
      if (!theEc)
      {
        ok = false;
	std::cerr <<"FFlNastranReader::resolveBushAttributes: Internal error, "
		  <<"EID="<< EID <<" PID="<< Ecc
		  <<", missing eccentricty."<< std::endl;
      }
      else
      {
#if FFL_DEBUG > 1
        std::cout <<" Transforming eccentricity vector for bushing element "
                  << EID <<", CD = "<< cit->second;
#endif
        ok &= this->transformPoint(theEc->offset.data(),cit->second,true);
#if FFL_DEBUG > 1
        std::cout << std::endl;
#endif
        theEc->offset.data().round(10); // round to 10 significant digits
      }
    }
  }
  else if (S > 0.0 && S <= 1.0)
  {
    // The bushing element is located along the line between node 1 and 2
    FaVec3 eVec;
    if (!this->getElementAxis(curElm,1,2,eVec))
      ok = false;
    else if (eVec.sqrLength() > 1.0e-16)
    {
      FFlPBUSHECCENT* myEc = CREATE_ATTRIBUTE(FFlPBUSHECCENT,"PBUSHECCENT",EID);
      myEc->offset = (eVec*S).round(10);
#if FFL_DEBUG > 1
      std::cout <<"Attribute PBUSHECCENT, ID = "<< EID
                <<", Fields: "<< myEc->offset << std::endl;
#endif
      myLink->addAttribute(myEc);
      ok &= curElm->setAttribute("PBUSHECCENT",EID);
    }
  }

  int OID = 0;
  if (CID > 0)
  {
    // A local element coordinate system is explicitly given
    std::map<int,CORD*>::iterator cit = cordSys.find(CID);
    if (cit == cordSys.end())
    {
      ListUI <<"\n *** Error: Coordinate system "<< CID <<" does not exist.\n";
      ok = false;
    }
    else if (cit->second->type != rectangular)
    {
      nWarnings++;
      ListUI <<"\n  ** Warning: A non-rectangular coordinate system "<< CID
             <<" is specified for bushing element "<< EID
             <<".\n              This is not supported."
             <<" The global system is used instead.\n";
    }
    else if (!this->computeTmatrix(CID,*(cit->second)))
      ok = false;
    else if (!myLink->getAttribute("PCOORDSYS",CID))
    {
      FFlPCOORDSYS* mySys = CREATE_ATTRIBUTE(FFlPCOORDSYS,"PCOORDSYS",CID);
      mySys->Origo = cit->second->Origo.round(10);
      mySys->Zaxis = cit->second->Zaxis.round(10);
      mySys->XZpnt = cit->second->XZpnt.round(10);
#if FFL_DEBUG > 1
      mySys->print();
#endif
      myLink->addAttribute(mySys);
    }
  }
  else if (CID < 0 || (OID = curElm->getAttributeID("PORIENT")) > 0)
  {
    // No element coordinate system is explicitly given, not even the global.
    // That is allowed only if the nodes are not coincident.
    FaVec3 x, y;
    if (!this->getElementAxis(curElm,1,2,x))
      ok = false;
    else if (x.sqrLength() < 1.0e-16)
    {
      ok = false;
      ListUI <<"\n *** Error: Bushing element "<< EID <<" has coincident nodes."
             <<"\n            An element coordinate system must then be given."
             <<"\n";
    }
    else if (OID > 0)
    {
      // An orientation vector is given either explicitly or trough a third node
      FFlPORIENT* theOr = GET_ATTRIBUTE(FFlPORIENT,"PORIENT",OID);
      if (CID < 0)
        // Orientation (local y-axis) is given by vector from node 1 and -CID
        ok &= this->getElementAxis(curElm,1,-CID,y);

      else
      {
        // Orientation (local y-axis) is given component-wise
        y = theOr->directionVector.getValue();
        IntMap::iterator cit = nodeCDID.find(curElm->getNodeID(1));
        if (cit != nodeCDID.end())
        {
          // This element node has a local coordinate system for displacements
          // in which also the orientation vector is defined
#if FFL_DEBUG > 1
          std::cout <<" Transforming orientation vector for bushing element "
                    << EID <<", CD = "<< cit->second;
#endif
          ok &= this->transformPoint(y,cit->second,true);
#if FFL_DEBUG > 1
          std::cout << std::endl;
#endif
        }
      }

      // Now find the local Z-axis of the bushing element
      theOr->directionVector = (x^y).truncate().round(10);
#if FFL_DEBUG > 1
      std::cout <<"Attribute PORIENT, ID = "<< EID
                <<", Fields: "<< theOr->directionVector << std::endl;
#endif
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::transformMassMatrix (FFlElementBase* elm,
					    const int CID, FaVec3* X)
{
  int EID = elm->getID();
  int PID = elm->getAttributeID("PMASS");

  FFlPMASS* theMass = NULL;
  if (PID > 0) theMass = GET_ATTRIBUTE(FFlPMASS,"PMASS",PID);
  if (!theMass)
  {
    std::cerr <<"FFlNastranReader::transformMassMatrix: Internal error, "
	      <<"EID="<< EID << ", missing mass property."<< std::endl;
    return false;
  }

  double M[6][6];
  size_t i, j, k;
  for (i = k = 0; i < 6; i++)
    for (j = 0; j <= i; j++)
    {
      if (k >= theMass->M.data().size())
        M[i][j] = 0.0;
      else
        M[i][j] = theMass->M.data()[k++];
      if (j < i) M[j][i] = M[i][j];
    }

#if FFL_DEBUG > 1
  std::cout <<"\nMass matrix for element "<< EID <<"\n";
  printMatrix6(M);
#endif

  FFlNode* n1 = myLink->getNode(elm->getNodeID(1));
  if (!n1)
  {
    ListUI <<"\n *** Error: Non-existing node "<< elm->getNodeID(1)
           <<" referred by mass element "<< EID <<".\n";
    return false;
  }

  if (X)
  {
    if (CID < 0)
      *X -= n1->getPos();
    else if (CID > 0)
      if (!this->transformPoint(*X,CID,true))
	return false;

    if (X->sqrLength() > 0.0)
    {
      FFaAlgebra::eccTransform6(M,*X);
#if FFL_DEBUG > 1
      std::cout <<"\nTransformed matrix:\n";
      printMatrix6(M);
#endif
    }
  }

  if (CID > 0)
    if (!this->transformSymmMatrix6(M,n1->getPos(),CID))
      return false;

  for (i = k = 0; i < 6; i++)
    for (j = 0; j <= i; j++)
      if (k < theMass->M.data().size())
        theMass->M.data()[k++] = round(M[i][j],10); // use 10 significant digits
      else if (M[i][j] != 0.0)
      {
        theMass->M.data().resize(k+1,0.0);
        theMass->M.data()[k++] = round(M[i][j],10); // use 10 significant digits
      }
      else
        k++;

  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::transformSymmMatrix6 (double mat[6][6], const FaVec3& X,
					     const int CID, const int node)
{
  std::map<int,CORD*>::iterator cit = cordSys.find(CID);
  if (cit == cordSys.end())
  {
    ListUI <<"\n *** Error: Coordinate system "<< CID <<" does not exist.\n";
    return false;
  }

  FaMat33 T; // Compute transformation matrix for CID at the global point X
  if (!this->getTmatrixAtPt(CID,*(cit->second),X,T))
    return false;

  // Perform a congruence transformation
  double* M[6] = { mat[0], mat[1], mat[2], mat[3], mat[4], mat[5] };
  if (!FFaAlgebra::congruenceTransform(M,T,2,node))
    return false;

#if FFL_DEBUG > 1
  std::cout <<"\nTransformed matrix:\n";
  printMatrix6(mat);
#endif
  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::transformVec3 (FaVec3& v, const FaVec3& X, const int CID)
{
  std::map<int,CORD*>::iterator cit = cordSys.find(CID);
  if (cit == cordSys.end())
  {
    ListUI <<"\n *** Error: Coordinate system "<< CID <<" does not exist.\n";
    return false;
  }

  FaMat33 T; // Compute transformation matrix for CID at the global point X
  if (!this->getTmatrixAtPt(CID,*(cit->second),X,T))
    return false;

  v = T * v; // Transform vector v to global coordinates

  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::resolveLoads ()
{
  bool ok = true;
  for (const LoadFaceMap::value_type& load : loadFace)
    ok &= this->resolveLoadFace(load.first,load.second);
  for (const std::pair<FFlLoadBase* const,int>& load : loadCID)
    ok &= this->resolveLoadDirection(load.first,load.second);

  loadFace.clear();
  loadCID.clear();
  return ok;
}

////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::resolveLoadFace (const FFlLoadBase* load,
                                        const std::pair<int,int>& nodes) const
{
  int EID, faceNum = 0;
  if (!load->getTarget(EID,faceNum))
    return false;

  FFlElementBase* elm = myLink->getElement(EID);
  if (!elm)
  {
    ListUI <<"\n *** Error: Non-existing element "<< EID
	   <<" referred by pressure load "<< load->getID() <<".\n";
    return false;
  }

  faceNum = elm->getFaceNum(nodes.first,nodes.second);
  if (faceNum < 1)
  {
    ListUI <<"\n *** Error: The nodes "<< nodes.first <<" and "<< nodes.second
	   <<" do not define a face on element "<< EID
	   <<", referred by pressure load "<< load->getID() <<".\n";
    return false;
  }

  const_cast<FFlLoadBase*>(load)->setTarget(EID,faceNum);

  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool FFlNastranReader::resolveLoadDirection (const FFlLoadBase* load, int CID)
{
  int EID, faceNum = 0;
  if (!load->getTarget(EID,faceNum))
    return false;

  const FFlCLOAD* cload = dynamic_cast<const FFlCLOAD*>(load);
  if (cload)
  {
    FFlNode* node = myLink->getNode(EID);
    if (!node)
    {
      ListUI <<"\n *** Error: Non-existing node "<< EID
             <<" referred by concentrated load.\n";
      return false;
    }

    FaVec3& Pval = const_cast<FFlCLOAD*>(cload)->P.data();
    if (!this->transformVec3(Pval,node->getPos(),CID))
      return false;

    Pval.round(10);
    return true;
  }

  return false; // Surface loads not yet implemented
}
