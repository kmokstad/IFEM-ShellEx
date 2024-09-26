// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_GROUP_BASE_H
#define FFL_GROUP_BASE_H

#include "FFlLib/FFlNamedPartBase.H"


class FFlGroupBase : public FFlNamedPartBase
{
protected:
  FFlGroupBase(int id, const std::string& groupName) : FFlNamedPartBase(id)
  {
    this->setName(groupName);
  }

public:
  virtual ~FFlGroupBase() {}

  virtual FFlGroupBase* clone() const  = 0;

  virtual void addElement(int ID, bool sortOnInsert = false) = 0;
  virtual void sortElements(bool removeDuplicates = false) = 0;

  virtual void calculateChecksum(FFaCheckSum* cs, int csMask) const = 0;
};

#endif

