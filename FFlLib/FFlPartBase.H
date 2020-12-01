// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_PART_BASE_H
#define FFL_PART_BASE_H

#include <cstddef>

class FFlTypeInfoSpec;
class FFaCheckSum;


/*!
  \brief Base class for all object types of a FE model.
*/

class FFlPartBase
{
protected:
  FFlPartBase(int ID) : myID(ID) {}
  FFlPartBase(const FFlPartBase&) = default;
  FFlPartBase& operator=(const FFlPartBase&) = delete;

public:
  virtual ~FFlPartBase() {}

  virtual FFlTypeInfoSpec* getTypeInfoSpec() const = 0;

  int getID() const { return myID; }
  void setID(int newID) { myID = newID; }

protected:
  void checksum(FFaCheckSum* cs, int csMask = 0) const;

private:
  int myID;
};


class FFlFEPartBaseLess
{
public:
  bool operator()(const FFlPartBase* lhs, const FFlPartBase* rhs) const
  { return lhs->getID() < rhs->getID(); }

  bool operator()(const FFlPartBase* lhs, const int rhs) const
  { return lhs->getID() < rhs; }

  bool operator()(const int lhs, const FFlPartBase* rhs) const
  { return lhs < rhs->getID(); }
};

#endif
