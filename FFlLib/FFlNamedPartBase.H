// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_NAMED_PART_BASE_H
#define FFL_NAMED_PART_BASE_H

#include <string>

#include "FFlLib/FFlPartBase.H"

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif


class FFlNamedPartBase : public FFlPartBase
{
public:
  enum {FFL_USED                = 0x0001,
	FFL_USED_MASK           = 0x000f,
	FFL_HAS_VIS_ELM         = 0x0010,
	FFL_HAS_VIS_ELM_MASK    = 0x00f0,
	FFL_HAS_HIDDEN_ELM      = 0x0100,
	FFL_HAS_HIDDEN_ELM_MASK = 0x0f00
  };

protected:
  FFlNamedPartBase(int ID) : FFlPartBase(ID), visibilityStatus(0) {}

  FFlNamedPartBase(const FFlNamedPartBase& obj)
    : FFlPartBase(obj), userName(obj.userName), visibilityStatus(0) {}

  FFlNamedPartBase& operator=(const FFlNamedPartBase&) = delete;

public:
  const std::string& getName() const { return userName; }
  void setName(const std::string& groupName);

  void addVisibilityStatus(int val) { visibilityStatus |= val; }
  int  getVisibilityStatus() const  { return visibilityStatus; }
  void resetVisibilityStatus()      { visibilityStatus = 0;    }

protected:
  void checksum(FFaCheckSum* cs, int csMask) const;

private:
  std::string userName;
  int visibilityStatus;
};

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
