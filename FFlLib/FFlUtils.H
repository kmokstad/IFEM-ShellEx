// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_UTILS_H
#define FFL_UTILS_H

#include <vector>
#include <map>

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE { class FFlLinkHandler; }
using FFlLinkHandler = FF_NAMESPACE::FFlLinkHandler;
#else
class FFlLinkHandler;
#endif


namespace FFl
{
  struct DepDOF
  {
    int node;
    int lDof;
    double coeff;
    DepDOF(int n = 1, int d = 1, double c = 1.0) : node(n), lDof(d), coeff(c) {}
  };

  using DepDOFs = std::vector<DepDOF>;
  using MPC     = std::map<short int,DepDOFs>;
  using MPCMap  = std::map<int,MPC>;

  bool convertMPCsToWAVGM(FFlLinkHandler* part, const MPCMap& mpcs);
}

#endif
