// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFaFilePath.C
  \brief Utilities for file path handling.
*/

#include <algorithm>
#include <cstring>
#include <cctype>

#include "FFaLib/FFaOS/FFaFortran.H"
#include "FFaLib/FFaOS/FFaFilePath.H"

#if defined(win32) || defined(win64)
//! \brief Wrong file separator on Windows.
#define WRONG_SLASH '/'
//! \brief Right file separator on Windows.
#define RIGHT_SLASH '\\'
#else
//! \brief Wrong file separator on UNIX.
#define WRONG_SLASH '\\'
//! \brief Right file separator on UNIX.
#define RIGHT_SLASH '/'
#endif


/*!
  \brief Converts file pathnames from UNIX to Windows syntax and vice versa.
*/

size_t FFaFilePath::checkName (char* path, size_t length)
{
  size_t newLength = length;
#if defined(win32) || defined(win64)
  // Check if we have a cygwin path: //D/models/bar ==> D:\models\bar
  //                                 /cygdrive/D/models/bar ==> D:\models\bar
  size_t lshift = 0;
  if (length > 12 && !strncmp(path,"/cygdrive/",10))
    lshift = 9;
  else if (length > 4 && path[0] == '/' && path[1] == '/')
    lshift = 1;
  if (lshift && isalpha(path[1+lshift]) && path[2+lshift] == '/')
  {
    path[0] = path[1+lshift];
    path[1] = ':';
    path[2] = RIGHT_SLASH;
    for (size_t i = 3+lshift; i < length; i++) path[i-lshift] = path[i];
    newLength = length-lshift; // The modified path is lshift characters shorter
  }
#endif
  for (size_t i = 0; i < newLength; i++)
    if (path[i] == WRONG_SLASH) path[i] = RIGHT_SLASH;

  return newLength;
}


/*!
   \brief C++ wrapper for FFaFilePath::checkName(char*,size_t).
*/

std::string& FFaFilePath::checkName (std::string& thePath)
{
  if (thePath.empty()) return thePath;

  // Touch the string in case it's pointer is still referring to another string
  // immediately after an assignment operator call, etc.
  char c = thePath[0]; thePath[0] = '\0'; thePath[0] = c;

  size_t newLength = checkName((char*)thePath.c_str(),thePath.length());
  if (newLength < thePath.length()) thePath.assign(thePath,0,newLength);

  return thePath;
}


/*!
  \brief Fortran wrapper for FFaFilePath::checkName(char*,size_t).
*/

SUBROUTINE(ffa_checkpath,FFA_CHECKPATH) (char* thePath, const int nchar)
{
  int newLength = FFaFilePath::checkName(thePath,nchar);
  if (newLength < nchar) memset(thePath+newLength,' ',nchar-newLength);
}


/*!
  \brief Modify the given path to using UNIX-style path separators.
*/

std::string& FFaFilePath::unixStyle (std::string& thePath)
{
  for (size_t i = 0; i < thePath.size(); i++)
    if (thePath[i] == '\\')
      thePath[i] = '/';

  return thePath;
}


/*!
  \brief Return a copy of the given path with  UNIX-style path separators.
*/

std::string FFaFilePath::unixStyle (const std::string& thePath)
{
  std::string newPath(thePath);
  return unixStyle(newPath);
}


/*!
  \brief Checks if the provided \a fileName has a path.
*/

bool FFaFilePath::hasPath (const std::string& fileName)
{
  return fileName.find_first_of("/\\") < fileName.size();
}


/*!
  \brief Returns the path part of \a fullPath.
  \details It is equivalent to the UNIX command "dirname",
  except that the trailing '/', after removing the file name part,
  is retained when \a keepTrailingSlash is true.
*/

std::string FFaFilePath::getPath (const std::string& fullPath,
                                  bool keepTrailingSlash)
{
  size_t slashPos = fullPath.find_last_of("/\\");

  // Ignore trailing '/'-character(s) in fullPath
  size_t prevSPos = fullPath.size();
  while (slashPos == prevSPos-1 && slashPos > 0)
  {
    prevSPos = slashPos;
    slashPos = fullPath.find_last_of("/\\",slashPos-1);
  }

  if (slashPos > fullPath.size())
    return std::string("");
  else if (keepTrailingSlash || slashPos == 0)
    return fullPath.substr(0,slashPos+1);
  else
    return fullPath.substr(0,slashPos);
}


/*!
  \brief Returns the file name part of \a fullPath.
*/

std::string FFaFilePath::getFileName (const std::string& fullPath)
{
  size_t slashPos = fullPath.find_last_of("/\\");
  if (slashPos < fullPath.size())
    return fullPath.substr(slashPos+1);
  else
    return fullPath;
}


/*!
  \brief Returns the file name without extension.
  \details Beware that if \a fName has a path,
  that will be included as well, unless \a removePath is true.
*/

std::string FFaFilePath::getBaseName (const std::string& fName, bool removePath)
{
  int slashPos = fName.find_last_of("/\\");
  if (slashPos > (int)fName.size()) slashPos = -1;
  size_t dotPos = fName.find_last_of('.');
  if (dotPos > fName.size() || (int)dotPos < slashPos)
    dotPos = fName.find_last_not_of(' ')+1;

  if (removePath && slashPos > -1)
    return fName.substr(slashPos+1,dotPos-slashPos-1);
  else if (dotPos < fName.size())
    return fName.substr(0,dotPos);
  else
    return fName;
}


/*!
  \brief Fortran wrapper for FFaFilePath::getBaseName(const std::string&,bool).
*/

SUBROUTINE(ffa_getbasename,FFA_GETBASENAME) (const char* fName,
#ifdef _NCHAR_AFTER_CHARARG
                                             const int nFchar,
                                             char* bName,
#else
                                             char* bName,
                                             const int nFchar,
#endif
                                             const int nBchar)
{
  std::string baseName = FFaFilePath::getBaseName(std::string(fName,nFchar),true);
  int newLength = baseName.size();
  if (newLength > nBchar) newLength = nBchar;
  memcpy(bName,baseName.c_str(),newLength);
  if (newLength < nBchar) memset(bName+newLength,' ',nBchar-newLength);
}


/*!
  \brief Returns the extension of \a fName.
*/

std::string FFaFilePath::getExtension (const std::string& fName)
{
  if (fName.empty())
    return std::string("");

  size_t slashPos = fName.find_last_of("/\\");
  if (slashPos > fName.size())
    slashPos = 0;

  size_t dotPos = fName.find_last_of('.');

  if (dotPos < fName.size()-1 && dotPos >= slashPos)
    return fName.substr(dotPos+1);
  else
    return std::string("");
}


/*!
  \brief Checks if \a ext is the extension of \a fName.
*/

bool FFaFilePath::isExtension(const std::string& fName, const std::string& ext)
{
  int i = fName.size();
  int j = ext.size();
  while (i > 1 && j > 0)
    if (fName[--i] != ext[--j])
      return false;

  return j == 0 && fName[i-1] == '.' ? true : false;
}


/*!
  \brief Appends \a fileName to \a path.
  \details If \a fileName begins with one or more "../" and there are
  a corresponding number of directory levels in \a path, the resulting pathname
  is modified such that, e.g., \code /a/bb/../ccc \endcode is reduced to
  \code /a/ccc \endcode
*/

std::string FFaFilePath::appendFileNameToPath(const std::string& path,
                                              const std::string& fileName)
{
  if (path.empty()) return fileName;
  if (fileName.empty()) return path;

  int iStart = fileName[0] == '.' && fileName.find_first_of("/\\") == 1 ? 2 : 0;
  int lPath  = path.size() - 1;
  int iSlash = path.find_last_of("/\\");
  int jSlash = iSlash;
  if (iSlash == lPath)
    jSlash = path.find_last_of("/\\",iSlash-1);
  else
    iSlash = -1;

  while (fileName.substr(iStart,2) == ".." && jSlash < lPath)
    if (path.substr(jSlash+1,2) == "..") break; // TT #3004
    else if (fileName[iStart+2] == '/' || fileName[iStart+2] == '\\')
    {
      iStart += 3;
      iSlash = jSlash;
      if (iSlash == 0) break;
      jSlash = path.find_last_of("/\\",iSlash-1);
    }
    else
      break;

  if (iStart > 0 && iSlash >= 0)
    return path.substr(0,iSlash+1) + fileName.substr(iStart);
  else if (iSlash == lPath)
    return path + fileName.substr(iStart);
  else
    return path + RIGHT_SLASH + fileName.substr(iStart);
}


/*!
  \brief Appends \a fileName to \a path.
*/

std::string& FFaFilePath::appendToPath (std::string& path,
                                        const std::string& fileName)
{
  path = appendFileNameToPath(path,fileName);

  return path;
}


/*!
  \brief Makes \a fileName an absolute pathname by prefixing it with \a absPath.
  \details If \a fileName already is an absolute path, that name is returned.
*/

std::string& FFaFilePath::makeItAbsolute (std::string& fileName,
                                          const std::string& absPath)
{
  if (isRelativePath(fileName))
    fileName = appendFileNameToPath(absPath,fileName);

  return fileName;
}


/*!
  \brief Replaces the current path in \a fileName with \a path.
*/

std::string& FFaFilePath::setPath (std::string& fileName,
                                   const std::string& path)
{
  fileName = appendFileNameToPath(path,getFileName(fileName));

  return fileName;
}


/*!
  \brief Returns the path separator for this system.
*/

char FFaFilePath::getPathSeparator ()
{
  return RIGHT_SLASH;
}


/*!
  \brief Adds \a ext to \a fileName, if not already present.
  \details If \a fileName has another extension, that extension
  is replaced by \a ext.
*/

std::string& FFaFilePath::addExtension (std::string& fName, const std::string& ext)
{
  if (!isExtension(fName,ext))
  {
    int slashPos = fName.find_last_of("/\\");
    if (slashPos > (int)fName.size()) slashPos = -1;
    size_t dotPos = fName.find_last_of('.');

    if (dotPos < fName.size()-1 && (int)dotPos > slashPos)
      fName.replace(dotPos+1,std::string::npos,ext);
    else
      fName.append('.' + ext);
  }

  return fName;
}


/*!
  \brief Checks \a fName for existing extension in \a exts.
  \details If not found, adds the first.
*/

std::string& FFaFilePath::addExtension (std::string& fName,
                                        const std::vector<std::string>& exts)
{
  if (exts.empty()) return fName;

  for (size_t i = 1; i < exts.size(); i++)
    if (isExtension(fName,exts[i])) return fName;

  return addExtension(fName,exts.front());
}


/*!
  \brief Determines if \a path is relative or absolute.
*/

bool FFaFilePath::isRelativePath (const std::string& path)
{
  if (path.empty()) return true;

  int i = 0;
  if (path.size() > 1 && isalpha(path[0]) && path[1] == ':')
  {
    if (path.size() < 3) return false;
    i = 2; // drive, e.g. A:
  }

  return path[i] != '/' && path[i] != '\\';
}


/*!
  Given the absolute current directory and an absolute file name,
  the relative file name is returned. For example,
  if the current directory is \code C:\foo\bar \endcode
  and we provide the filename \code C:\foo\whee\text.txt \endcode
  then this method will return \code ..\whee\text.txt \endcode

  \author Rob Fisher, rfisher@iee.org http://come.to/robfisher
*/

std::string FFaFilePath::getRelativeFilename (const std::string& currentDirectory,
                                              const std::string& absoluteFilename)
{
  //! \cond DO_NOT_DOCUMENT
  // The number of characters at the start of an absolute filename. e.g. in DOS,
  // absolute filenames start with "X:\" so this value should be 3, in UNIX they
  // start with "/" so this value should be 1.
#if defined(win32) || defined(win64)
#define ABSOLUTE_NAME_START 3
#else
#define ABSOLUTE_NAME_START 1
#endif
  // On windows, the string comparison should be case insensitive
#if defined(win32) || defined(win64)
#define EQUAL(a,b) (toupper(a) == toupper(b))
#else
#define EQUAL(a,b) (a == b)
#endif
  //! \endcond

  int cdLen = currentDirectory.size();
  int afLen = absoluteFilename.size();

  // Make sure the names are not too short
  // Bugfix 20/11/2015: If current directory is empty, return the absolute path
  if (cdLen < ABSOLUTE_NAME_START) return absoluteFilename;
  // If both are relative path names, start over with an arbitrary common root
  if (isRelativePath(currentDirectory) && isRelativePath(absoluteFilename))
  {
#if defined(win32) || defined(win64)
    std::string root("A:\\tmp\\");
#else
    std::string root("/tmp/");
#endif
    return getRelativeFilename(root+currentDirectory,root+absoluteFilename);
  }
  if (afLen <= ABSOLUTE_NAME_START) return "";

  // Handle DOS names that are on different drives:
  if (!EQUAL(currentDirectory[0],absoluteFilename[0])) return absoluteFilename;

  // They are on the same drive, find out how much of the current directory
  // is in the absolute filename
  int i = ABSOLUTE_NAME_START;
  while (i < afLen && i < cdLen)
    if (EQUAL(currentDirectory[i],absoluteFilename[i]))
      i++;
    else
      break;

  if (i == cdLen)
  {
    if (absoluteFilename[i-1] == RIGHT_SLASH)
      // the whole current directory name is in the file name,
      // so we just trim off the current directory name to get the
      // current file name.
      return absoluteFilename.substr(i);

    if (absoluteFilename[i] == RIGHT_SLASH)
      // a directory name might have a trailing slash but a relative
      // file name should not have a leading one...
      return absoluteFilename.substr(i+1);
  }

  // The file is not in a child directory of the current directory, so we
  // need to step back the appropriate number of parent directories by
  // using "..\"s.  First find out how many levels deeper we are than the
  // common directory
  int afMarker = i;
  int levels = 1;

  // count the number of directory levels we have to go up to get to the
  // common directory
  while (i < cdLen)
    if (currentDirectory[++i] == RIGHT_SLASH)
      // make sure it's not a trailing slash
      if (currentDirectory[++i] != '\0')
        levels++;

  // move the absolute filename marker back to the start of the directory name
  // that it has stopped in.
  while (afMarker > 0 && absoluteFilename[afMarker-1] != RIGHT_SLASH)
    afMarker--;

  // add the appropriate number of "..\"s.
  std::string relativeFilename;
  for (i = 0; i < levels; i++)
  {
    relativeFilename += "..";
    relativeFilename += RIGHT_SLASH;
  }

  // copy the rest of the filename into the result string
  relativeFilename += absoluteFilename.substr(afMarker);

  return relativeFilename;
}


/*!
  \brief Makes \a fileName a relative pathname with respect to \a absPath.
*/

std::string& FFaFilePath::makeItRelative (std::string& fileName,
                                          const std::string& absPath)
{
  if (!absPath.empty() && !fileName.empty())
    fileName = getRelativeFilename(absPath,fileName);

  return fileName;
}


/*!
  Utility method to get rid of characters we don't want to have in a file name.
  If \a removePath is true, all characters up to and including the last '/' or
  '\\' character are removed, assuming they are the directory path.
*/

std::string FFaFilePath::distillName (const std::string& filePath,
                                      bool removePath)
{
  std::string fName = removePath ? getFileName(filePath) : filePath;

  // Loop over fName to replace all unwanted characters by "_".
  // Maybe we are being a bit rough, but that is always safer.

  for (size_t i = 0; i < fName.size(); i++)
    if (fName[i] == ' ' ||
	fName[i] == ':' ||
	fName[i] == ';' ||
	fName[i] == ',' ||
	fName[i] == '!' ||
	fName[i] == '#' ||
	fName[i] == '$' ||
	fName[i] == '@' ||
	fName[i] == '/' ||
	fName[i] == '\\' ||
	fName[i] == '{' ||
	fName[i] == '}' ||
	fName[i] == '[' ||
	fName[i] == ']' ||
	fName[i] == ')' ||
	fName[i] == '(' ||
	fName[i] == '?' ||
	fName[i] == '+' ||
	fName[i] == '*' ||
	fName[i] == '=' ||
	fName[i] == '&' ||
	fName[i] == '%' ||
	fName[i] == '\'' ||
	fName[i] == '\"' ||
	fName[i] == '|' ||
	fName[i] == '~')
      fName[i] = '_';

  return fName;
}
