// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_READER_BASE_H
#define FFL_READER_BASE_H

#include <cstddef>

class FFlLinkHandler;
class FFaProfiler;


class FFlReaderBase
{
public:
  FFlReaderBase(FFlLinkHandler* p) : myLink(p), myProfiler(NULL) {}
  virtual ~FFlReaderBase() {}

protected:
  static int searchKeyword(const char* fileName, const char* keyWord,
                           int maxLines = 100);
  static int searchKeyword(const char* fileName, const char** keyWord,
                           int maxLines = 100);

  FFlLinkHandler* myLink;
  FFaProfiler*    myProfiler;
};

#endif
