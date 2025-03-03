// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_TYPE_INFO_SPEC_H
#define FFL_TYPE_INFO_SPEC_H

#include <string>

#include "FFaLib/FFaPatterns/FFaSingelton.H"


#define FFL_TYPE_INFO(class) \
virtual FFlTypeInfoSpec* getTypeInfoSpec() const \
{ return FFaSingelton<FFlTypeInfoSpec,class>::instance(); }


class FFlTypeInfoSpec
{
public:
  FFlTypeInfoSpec() : myCathegory(OTHER) {}

  const std::string& getTypeName() const { return myTypeName; }
  void setTypeName(const std::string& n) { myTypeName = n; }

  const std::string& getDescription() const { return myDescription; }
  void setDescription(const std::string& d) { myDescription = d; }

  enum Cathegory{SOLID_ELM, SHELL_ELM, BEAM_ELM, CONSTRAINT_ELM,
                 STRC_ELM, OTHER_ELM,
                 MATERIAL_PROP, MASS_PROP, GEOMETRY_PROP,
                 STRC_PROP, VISUAL_PROP, OTHER_PROP,
                 COORD_SYS, USER_DEF_GROUP, NODE,
		 LOAD, OTHER};

  Cathegory getCathegory() const { return myCathegory; }
  void setCathegory(Cathegory c) { myCathegory = c; }

private:
  std::string myTypeName;
  std::string myDescription;
  Cathegory   myCathegory;
};

#endif
