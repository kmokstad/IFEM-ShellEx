// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlPWAVGM.H"


FFlPWAVGM::FFlPWAVGM(int id) : FFlAttributeBase(id)
{
  this->addField(refC);
  refC = 0;

  for (int i = 0; i < 6; i++)
  {
    this->addField(indC[i]);
    indC[i] = 0;
  }
  this->addField(weightMatrix);
}


FFlPWAVGM::FFlPWAVGM(const FFlPWAVGM& obj) : FFlAttributeBase(obj)
{
  this->addField(refC);
  refC = obj.refC;

  for (int i = 0; i < 6; i++)
  {
    this->addField(indC[i]);
    indC[i] = obj.indC[i];
  }

  this->addField(weightMatrix);
  weightMatrix = obj.weightMatrix;
}


bool FFlPWAVGM::isIdentic(const FFlAttributeBase* otherAttrib) const
{
  const FFlPWAVGM* other = dynamic_cast<const FFlPWAVGM*>(otherAttrib);
  if (!other) return false;

  if (this->refC != other->refC) return false;

  for (int i = 0; i < 6; i++)
    if (this->indC[i] != other->indC[i]) return false;

  if (this->weightMatrix != other->weightMatrix) return false;

  return true;
}


void FFlPWAVGM::init()
{
  FFlPWAVGMTypeInfoSpec::instance()->setTypeName("PWAVGM");
  FFlPWAVGMTypeInfoSpec::instance()->setDescription("Weighted average motion properties");
  FFlPWAVGMTypeInfoSpec::instance()->setCathegory(FFlTypeInfoSpec::OTHER_PROP);

  AttributeFactory::instance()->registerCreator
    (FFlPWAVGMTypeInfoSpec::instance()->getTypeName(),
     FFaDynCB2S(FFlPWAVGM::create,int,FFlAttributeBase*&));
}


FFlAttributeBase* FFlPWAVGM::removeWeights (const std::vector<int>& nodeIndices, size_t nNod)
{
  FFlAttributeBase* newAtt = AttributeFactory::instance()->create("PWAVGM",this->getID());
  FFlPWAVGM* myAtt = static_cast<FFlPWAVGM*>(newAtt);

  // Scale the remaining weights such that they maintain the same sum
  size_t nRow = weightMatrix.getValue().size() / (nNod-1);
  size_t newNod = nNod-1 - nodeIndices.size();
  std::vector<double> newW, remW(nRow,0.0);
  newW.reserve(nRow*newNod);
  size_t i, j, k = 0;
  std::vector<double>::const_iterator wit = weightMatrix.getValue().begin();
  for (i = 1; i < nNod; i++, wit += nRow)
    if (k >= nodeIndices.size() || (int)i < nodeIndices[k])
      newW.insert(newW.end(),wit,wit+nRow);
    else for (++k, j = 0; j < nRow; j++)
      remW[j] += *(wit+j);

  for (i = k = 0; i < newW.size(); i += nRow)
    for (j = 0; j < nRow; j++, k++)
      newW[k] += remW[j] / double(newNod);

  myAtt->weightMatrix.data().swap(newW);

  // Compute new component indices
  int iC;
  for (j = 0; j < 6; j++)
    if ((iC = indC[j].getValue()) > 0)
      myAtt->indC[j] = (iC-1)*newNod/nNod + 1;

  myAtt->refC = refC.getValue();

  return newAtt;
}
