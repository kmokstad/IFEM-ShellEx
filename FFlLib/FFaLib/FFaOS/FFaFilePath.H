// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFaFilePath.H
  \brief Utilities for file path handling.
*/

#ifndef FFA_FILE_PATH_H
#define FFA_FILE_PATH_H

#include <string>
#include <vector>


namespace FFaFilePath //! Utilities for file path handling
{
  // In place convertion Win32 <-> Unix

  size_t checkName (char* thePath, size_t length);
  std::string& checkName (std::string& thePath);
  std::string& unixStyle (std::string& thePath);
  std::string  unixStyle (const std::string& thePath);

  // Some checking

  bool hasPath        (const std::string& fileName);
  bool isRelativePath (const std::string& pathName);
  bool isExtension    (const std::string& fileName, const std::string& ext);

  // Basic string operations on paths

  std::string getPath      (const std::string& fullPath, bool keepSlash = true);
  std::string getFileName  (const std::string& fullPath);
  std::string getExtension (const std::string& fileName);
  std::string getBaseName  (const std::string& fileName, bool removePath = false);

  std::string appendFileNameToPath (const std::string& path, const std::string& fName);
  std::string& appendToPath   (std::string& path,  const std::string& fName);
  std::string& makeItAbsolute (std::string& fName, const std::string& path);
  std::string& makeItRelative (std::string& fName, const std::string& path);
  std::string& setPath        (std::string& fName, const std::string& path);

  char getPathSeparator ();

  std::string& addExtension (std::string& fileName, const std::string& ext);
  std::string& addExtension (std::string& fileName, const std::vector<std::string>& exts);

  //! \brief Returns the relative path name with respect to an absolute path.
  std::string getRelativeFilename (const std::string& currentDirectory,
                                   const std::string& absoluteFilename);

  //! \brief Checks for illegal characters in file names.
  std::string distillName (const std::string& filePath, bool removePath = false);
}

#endif
