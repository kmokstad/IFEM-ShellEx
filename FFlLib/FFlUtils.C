// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlUtils.H"
#include "FFlLib/FFlLinkHandler.H"
#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlFEParts/FFlPWAVGM.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"


bool FFl::convertMPCsToWAVGM (FFlLinkHandler* part, const FFl::MPCMap& mpcs)
{
  using Doubles       = std::vector<double>;
  using DoublesMap    = std::map<int,Doubles>;
  using DoublesMapMap = std::map<int,DoublesMap>;

  // Create WAVGM elements for multi-point constraints with common slave node
  for (const std::pair<const int,MPC>& mpcGroup : mpcs)
  {
    // Find the element nodes
    std::vector<int> nodes = { mpcGroup.first };
    for (const std::pair<const short int,DepDOFs>& mpc : mpcGroup.second)
      for (const DepDOF& dof : mpc.second)
        if (std::find(nodes.begin(),nodes.end(),dof.node) == nodes.end())
          nodes.push_back(dof.node);

#if FFL_DEBUG > 1
    std::cout <<"\nWAVGM element nodes:";
    for (int n : nodes) std::cout <<" "<< n;
    std::cout << std::endl;
#endif

    size_t nRow = 0;
    DoublesMapMap dofWeights;
    for (const std::pair<const short int,DepDOFs>& mpc : mpcGroup.second)
      if (mpc.first > 0 && mpc.first < 7)
      {
        // Find the weight matrix associated with this slave DOF
        DoublesMap& dofWeight = dofWeights[mpc.first];
        for (const DepDOF& dof : mpc.second)
        {
          Doubles& weights = dofWeight[dof.lDof];
          weights.resize(nodes.size()-1,0.0);
          for (size_t iNod = 1; iNod < nodes.size(); iNod++)
            if (dof.node == nodes[iNod])
            {
              weights[iNod-1] = dof.coeff;
              break;
            }
        }
        nRow += 6; // Assuming all DOFs in the master node is referred
#if FFL_DEBUG > 1
        std::cout <<"Weight matrix associated with slave dof "<< mpc.first;
        for (const std::pair<const int,Doubles>& weights : dofWeight)
        {
          std::cout <<"\n\t"<< weights.first <<":";
          for (double c : weights.second) std::cout <<" "<< c;
        }
        std::cout << std::endl;
#endif
      }

    int refC = 0;
    size_t indx = 1;
    size_t nMst = nodes.size() - 1;
    int indC[6] = { 0, 0, 0, 0, 0, 0 };
    Doubles weights(nRow*nMst,0.0);
    for (const std::pair<const int,DoublesMap>& dofw : dofWeights)
    {
      refC = 10*refC + dofw.first; // Compressed slave DOFs identifier
      indC[dofw.first-1] = indx;   // Index to first weight for this slave DOF
      for (const std::pair<const int,Doubles>& dof : dofw.second)
        for (size_t j = 0; j < dof.second.size(); j++)
          weights[indx+6*j+dof.first-2] = dof.second[j];
      indx += 6*nMst;
    }

    int id = part->getNewElmID();
    FFlElementBase* newElem = ElementFactory::instance()->create("WAVGM",id);
    if (!newElem)
    {
      ListUI <<"\n *** Error: Failure creating WAVGM element "<< id <<".\n";
      return false;
    }

    newElem->setNodes(nodes);
    if (!part->addElement(newElem))
      return false;

    FFlAttributeBase* myAtt = AttributeFactory::instance()->create("PWAVGM",id);
    FFlPWAVGM* newAtt = static_cast<FFlPWAVGM*>(myAtt);
    newAtt->refC = -refC; // Hack: Negative refC means explicit constraints
    newAtt->weightMatrix.data().swap(weights);
    for (int j = 0; j < 6; j++)
      newAtt->indC[j] = indC[j];

    if (part->addUniqueAttributeCS(myAtt))
      newElem->setAttribute(myAtt);
    else
      return false;
  }

  return true;
}
