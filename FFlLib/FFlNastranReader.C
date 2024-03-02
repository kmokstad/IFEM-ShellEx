// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <cstring>
#include <cctype>

#include "FFlLib/FFlIOAdaptors/FFlNastranReader.H"
#include "FFlLib/FFlIOAdaptors/FFlReaders.H"
#include "FFlLib/FFlLinkHandler.H"
#include "FFlLib/FFlElementBase.H"
#include "FFlLib/FFlGroup.H"

#include "FFaLib/FFaAlgebra/FFaCheckSum.H"
#include "FFaLib/FFaOS/FFaFilePath.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"
#ifdef FFL_TIMER
#include "FFaLib/FFaProfiler/FFaProfiler.H"
#endif
#include "Admin/FedemAdmin.H"

#if _MSC_VER > 1300
#define strncasecmp _strnicmp
#endif

#ifdef FFL_TIMER
#define START_TIMER(function) myProfiler->startTimer(function);
#define STOPP_TIMER(function) myProfiler->stopTimer(function);
#else
#define START_TIMER(function)
#define STOPP_TIMER(function)
#endif

#define MAX_HEADER_LINES 1000

static std::string mainPath;

static bool identFoundSet = false;
static bool procOK        = true; // Set to false if parsing errors detected
static int  startBulk     = 0;    // Used to remember where the bulk data starts

int FFlNastranReader::nWarnings = 0;
int FFlNastranReader::nNotes = 0;

static const char* bulkIdent1 = "BEGIN BULK";
static const char* bulkIdent2 = "GRID";
static const char* bulkIdent3 = "Text Input for Bulk Data";
static const char* bulkIdent4 = "INCLUDE";
static const char* setIdent   = "SET";
static const char* endOfBulk  = "ENDDATA";

//! Convenience function to check for bulk data.
static char isBulkData(const char* line)
{
  if (strncasecmp(line,bulkIdent1,strlen(bulkIdent1)) == 0)
    return 'y'; // Yes, we have bulk data
  else if (strncmp(line,bulkIdent2,strlen(bulkIdent2)) == 0)
    return 'p'; // We should parse this line too (GRID entry)
  else if (strncmp(line+9,bulkIdent3,strlen(bulkIdent3)) == 0)
    return 'y'; // Yes, we have bulk data
  else if (strncmp(line,bulkIdent4,strlen(bulkIdent4)) == 0)
    return 'p'; // We should parse this line too (INCLUDE entry)

  static const char* allKeys[55] = {
    "ASET",  "ASET1", "BAROR", "BEAMOR", "CBAR",
    "CBEAM", "CBUSH", "CHEXA", "CELAS1", "CELAS2",
    "CONROD","CONM1", "CONM2", "CORD1C", "CORD1R",
    "CORD1S","CORD2C","CORD2R","CORD2S", "CPENTA",
    "CQUAD4","CQUAD8","CROD",  "CTETRA", "CTRIA3",
    "CTRIA6","CWELD", "FORCE", "GRDSET", "GRID",
    "MAT1",  "MAT2",  "MAT8",  "MAT9",   "MOMENT",
    "PBAR",  "PBARL", "PBEAM", "PBEAML", "PBUSH",
    "PCOMP", "PELAS", "PLOAD2","PLOAD4", "PROD",
    "PSHELL","PSOLID","PWELD", "QSET1",  "RBAR",
    "RBE2",  "RBE3",  "SPC",   "SPC1",    NULL };

  // Check all legal bulk entry keywords,
  // in case they are not listed in the usual order
  for (const char** keyw = allKeys; *keyw; ++keyw)
  {
    size_t n = strlen(*keyw);
    if (strncmp(line,*keyw,n) == 0)
      if (isspace(line[n]) || line[n] == '*' || line[n] == ',')
        return 's'; // Probably bulk data use this as starting line when parsing
  }
  return 0;
}


FFlNastranReader::FFlNastranReader (FFlLinkHandler* link, const int startHere)
  : FFlReaderBase(link), lineCounter(startHere)
{
#ifdef FFL_TIMER
  myProfiler = new FFaProfiler("NastranReader profiler");
  myProfiler->startTimer("FFlNastranReader");
#endif
#ifdef FFL_DEBUG
  std::cout <<"FFlNastranReader: start of bulk data: "<< startHere << std::endl;
#endif

  gridDefault = NULL;
  barDefault  = NULL;
  beamDefault = NULL;
  attChkSum   = NULL;
  sizeOK      = true;
}


FFlNastranReader::~FFlNastranReader ()
{
  for (std::pair<const int,CORD*>& cs : cordSys) delete cs.second;
  for (std::pair<const int,BEAMOR*>& bor : bOri) delete bor.second;
  for (std::pair<const int,FaVec3*>& mx : massX) delete mx.second;

  delete gridDefault;
  delete barDefault;
  delete beamDefault;
  delete attChkSum;

#ifdef FFL_TIMER
  myProfiler->stopTimer("FFlNastranReader");
  myProfiler->report();
  delete myProfiler;
#endif
}


void FFlNastranReader::init ()
{
  FFlReaders::instance()->registerReader("Nastran Bulk Data","nas",
					 FFaDynCB2S(FFlNastranReader::readerCB,const std::string&,FFlLinkHandler*),
					 FFaDynCB2S(FFlNastranReader::identifierCB,const std::string&,int&),
					 "Nastran Bulk Data reader v2.0",
					 FedemAdmin::getCopyrightString());

  FFlReaders::instance()->addExtension("Nastran Bulk Data","bdf");
  FFlReaders::instance()->setDefaultReader("Nastran Bulk Data");
}


void FFlNastranReader::identifierCB (const std::string& fname, int& positiveId)
{
  if (fname.empty()) return;

  startBulk = 0;
  positiveId = -1;
  identFoundSet = false;
  std::ifstream fs(fname.c_str());
  if (!fs) return;

  char line[256];
  for (int lCounter = 0; !fs.eof() && lCounter < MAX_HEADER_LINES; lCounter++)
  {
    fs.getline(line,255);
    char answer = isBulkData(line);
    if (!startBulk)
    {
      if (answer == 'y')
        startBulk = lCounter+1; // BEGIN BULK detected, start parsing next line
      else if (answer == 'p' || answer == 's')
        startBulk = lCounter; // Parse this line too (valid bulk entry detected)
    }
    if (answer == 'y' || answer == 'p')
    {
      positiveId = lCounter+1;
      return;
    }
    else if (!answer && !strncmp(line,setIdent,strlen(setIdent)))
      identFoundSet = true;
  }
  positiveId = 0;
  if (!identFoundSet && !startBulk)
    return;

  // No Nastran file keyword were found among the first MAX_HEADER_LINES lines.
  // However, a SET keyword was found, so continue the search throughout the
  // whole file to find the starting line of the bulk data section.
  for (int lCounter = MAX_HEADER_LINES; !fs.eof() && !positiveId; lCounter++)
  {
    fs.getline(line,255);
    char answer = isBulkData(line);
    if (!startBulk)
    {
      if (answer == 'y')
        startBulk = lCounter+1; // BEGIN BULK detected, start parsing next line
      else if (answer == 'p' || answer == 's')
        startBulk = lCounter; // Parse this line too (valid bulk entry detected)
    }
    if (answer == 'y' || answer == 'p')
      positiveId = lCounter+1;
  }
}


void FFlNastranReader::readerCB (const std::string& fname, FFlLinkHandler* link)
{
  nWarnings = nNotes = 0;
  FFlNastranReader reader(link,startBulk);
  mainPath = FFaFilePath::getPath(fname);
#ifdef FFL_DEBUG
  std::cout <<"FFlNastranReader: fileName = \""<< fname <<"\"\n"
	    <<"FFlNastranReader: mainPath = \""<< mainPath <<"\""<< std::endl;
#endif
  bool stillOk = reader.read(fname);
  bool setsOk  = true;

  if (identFoundSet)
    setsOk = reader.processSet(fname,startBulk);
  else if (reader.setsArePresent(fname))
    setsOk = reader.processSet(fname,startBulk);

  if (!setsOk)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: Invalid formatting of Nastran SET definitions.\n "
           <<"             Resulting FFlGroups can not be handled correctly.\n";
  }

  if (!reader.resolve(stillOk))
    link->deleteGeometry(); // Parsing failure, delete all FE data
  else if (nWarnings+nNotes > 0)
    ListUI <<"\n  ** Parsing FE data file \""<< fname
           <<"\" succeeded.\n     However, "<< nWarnings
           <<" warning(s) and "<< nNotes <<" note(s) were reported.\n"
           <<"     Please review the message(s) and check the FE data file.\n";
}


bool FFlNastranReader::read (const std::string& fname, bool includedFile)
{
  // If the fname is a relative path, make it into an absolute path, assuming
  // the given fname is relative to the location of the main bulk data file
  std::string fileName(fname);
  if (includedFile)
    if (FFaFilePath::isRelativePath(fname) && !mainPath.empty())
      fileName = FFaFilePath::appendFileNameToPath(mainPath,fname);
  FFaFilePath::checkName(fileName);

  std::ifstream fs(fileName.c_str());
  if (!fs)
  {
    ListUI <<"\n *** Error: Can not open Nastran bulk data file "
           << fileName <<"\n";
    return false;
  }

  // Skip the first 'lineCounter' lines (contain non-bulk Nastran entries)
  char line[256];
  bool dmapIsIncluded = false;
  for (int l = 0; l < lineCounter; l++)
    if (fs.eof())
    {
      ListUI <<"\n *** Error: Premature end-of-file encountered\n";
      return false;
    }
    else
    {
      fs.getline(line,255);
      // Check if this FE part is reduced externally by Nastran
      if (strncasecmp(line,"ASSIGN",6) == 0)
        this->processAssignFile(line);
      else if (strncasecmp(line,"INCLUDE",7) == 0)
        dmapIsIncluded = true;
    }

  size_t numOP2 = myLink->getOP2files().size();
  if (!dmapIsIncluded || numOP2 < 3)
    myLink->clearOP2files();
  else
  {
    nNotes++;
    ListUI <<"\n   * Note: "<< numOP2 <<" OP2-files were detected.\n"
	   <<"           The FE part is assumed to be externally reduced.\n";
  }

#ifdef FFL_DEBUG
  std::cout <<"FFlNastranReader: starting bulk data parsing at line "
	    << lineCounter+1 << std::endl;
#endif
  return this->read(fs);
}


bool FFlNastranReader::read (std::istream& is)
{
  START_TIMER("read")

  ignoredBulk.clear();
  sxErrorBulk.clear();
  BulkEntry entry;
  bool stillOK = procOK = true;
  while (sizeOK && (stillOK = this->getNextEntry(is,entry)))
    if (entry.name == endOfBulk)
      break;
    else if (!entry.cont.empty())
      ucEntries.push_back(entry);
    else if (!this->processThisEntry(entry))
      procOK = false;

  if (!sizeOK || !stillOK)
  {
    ListUI <<" *** Parsing Nastran bulk data aborted due to the above error.\n";
    if (!entry.name.empty())
    {
      ListUI <<"     Entry causing the trouble: \""<< entry.name
             <<"\"\n     Data fields read:";
      for (const std::string& fld : entry.fields) ListUI <<" \""<< fld <<"\"";
      ListUI <<"\n";
    }
  }

  if (!ignoredBulk.empty() || !sxErrorBulk.empty())
    ListUI <<"\n";

  for (const std::pair<const std::string,int>& igb : ignoredBulk)
  {
    nNotes++;
    ListUI <<"   * Note: "<< igb.second
           <<" instances of the unsupported bulk-entry \""<< igb.first
           <<"\" have been ignored.\n";
  }

  for (const std::pair<const std::string,int>& sxb : sxErrorBulk)
  {
    nNotes++;
    ListUI <<"   * Note: Syntax error has been detected in "<< sxb.second
           <<" instances of the bulk-entry \""<< sxb.first <<"\"\n";
  }

#ifdef FFL_DEBUG
  std::cout <<"FFlNastranReader: processed "<< lineCounter <<" lines (done)."
	    << std::endl;
#endif
  STOPP_TIMER("read")
  return sizeOK && stillOK && procOK;
}


bool FFlNastranReader::setsArePresent (const std::string& fname)
{
  std::ifstream fs(fname.c_str());
  if (!fs) return false;

  char line[256];
  for (int lCounter = 0; !fs.eof() && lCounter < MAX_HEADER_LINES; lCounter++)
  {
    fs.getline(line,255);
    if (strncmp(line,setIdent,strlen(setIdent)) == 0)
      return true;
    else if (isBulkData(line))
      return false;
  }
  return false;
}


bool FFlNastranReader::resolve (bool stillOK)
{
  START_TIMER("resolve")

  if (!ucEntries.empty())
  {
    // Try to process the uncompleted entries also
    nWarnings++;
    ListUI <<"\n  ** Warning: There are uncompleted bulk entries\n";
    for (BulkEntry& entry : ucEntries)
    {
      ListUI <<"              "<< entry.name;
      if (entry.fields.empty())
        ListUI <<"  nField = 0";
      else
        ListUI <<" "<< entry.fields.front()
               <<" ...  nField = "<< entry.fields.size();
      ListUI <<", Continuation Field = \""<< entry.cont;
      if (!stillOK)
        ListUI <<"\" (not processed)\n";
      else
      {
        ListUI <<"\"\n";
        stillOK = this->processThisEntry(entry);
      }
    }
  }

#ifdef FFL_DEBUG
  std::cout <<"FFlNastranReader: resolving coordinates ..."<< std::endl;
#endif
  if (stillOK) stillOK = this->resolveCoordinates();
#ifdef FFL_DEBUG
  std::cout <<"FFlNastranReader: resolving attributes ..."<< std::endl;
#endif
  if (stillOK) stillOK = this->resolveAttributes();
#ifdef FFL_DEBUG
  std::cout <<"FFlNastranReader: resolving loads ..."<< std::endl;
#endif
  if (stillOK) stillOK = this->resolveLoads();
#ifdef FFL_DEBUG
  std::cout <<"FFlNastranReader: resolve done."<< std::endl;
#endif

  STOPP_TIMER("resolve")
  return stillOK;
}


bool FFlNastranReader::getNextEntry (std::istream& is, BulkEntry& entry)
{
  START_TIMER("getNextEntry")

  int status = 2;
  std::string field;
  static std::string last_entry;

  // Read the first field containing the name of the bulk-entry
  while (field.empty() && status == 2)
    if (!(status = this->getNextField(is,field)))
    {
      entry.name.erase();
      STOPP_TIMER("getNextEntry")
      return false;
    }
    else if (field == endOfBulk)
    {
      entry.name = field;
      STOPP_TIMER("getNextEntry")
      return true;
    }

  if (field.empty())
    if (last_entry == "ADAPT" || last_entry == "OUTPUT")
    {
      // Ignore continuations of these entries (to avoid errors only)
      is.ignore(BUFSIZ,'\n');
      STOPP_TIMER("getNextEntry")
      return this->getNextEntry(is,entry);
    }

  // Determine the field format
  entry.ffmt = SMALL_FIELD;
  char lChar = field.empty() ? ' ' : field[field.length()-1];
  if (lChar == ',')
  {
    entry.ffmt = FREE_FIELD;
    field.erase(field.end()-1);
    if (!field.empty() && field[field.length()-1] == '*') // Large Free Field
      field.erase(field.end()-1);
  }
  else if (lChar == '*')
  {
    entry.ffmt = LARGE_FIELD;
    field.erase(field.end()-1);
  }

#if FFL_DEBUG > 3
  std::cout <<"FFlNastranReader: entry=\""<< field <<"\" format="<< entry.ffmt
	    << std::endl;
#endif

  // Lambda function with hacks needed to deal with non-standard continuations
  auto&& isContinuation = [](const std::string& field, const std::string& cont,
                             FieldFormat format) -> bool
  {
    // Direct match of continuation field (this is by the book)
    if (field == cont)
      return true;

    // But then...
    // An empty field matches also a continuation field with only a '+'
    if (field.empty())
      return cont == "+" ? true : false;

    // Some files leave out the initial '*' or '+' in the continuation field
    // so see if we get a match by adding a leading '*'/'+' to the field
    bool nonStandard = false;
    if (field == std::string(format == LARGE_FIELD ? "*" : "+") + cont)
      nonStandard = true;

    // Some large-field formatted files (e.g., from Strand7) might
    // have continuation fields with a leading '+' instead of '*',
    // but with '*' in the first field of the continuation line
    else if (format == LARGE_FIELD && field[0] == '*' && cont[0] == '+')
      nonStandard = field.substr(1) == cont.substr(1);

    if (nonStandard && nWarnings < 1)
    {
      nWarnings++;
      ListUI <<"\n  ** Warning: Bulk-data file may have inconsistent "
             <<"continuation fields,\n              Assuming leading +\n";
    }

    return nonStandard;
  };

  std::vector<BulkEntry>::iterator beit;
  if (field == "INCLUDE")
    entry.ffmt = FREE_FIELD; // Read include filename in free_field format

  else for (beit = ucEntries.begin(); beit != ucEntries.end(); ++beit)

    // Check if we have found the continuation of any of the uncompleted entries
    if (isContinuation(field,beit->cont,beit->ffmt))
    {
      // Yes, the current line contains the continuation of entry *beit
      // Now continue to read data fields, if any, for this uncompleted entry
#if FFL_DEBUG > 3
      std::cout <<"FFlNastranReader: continuing on entry "<< beit->name
                <<" "<< beit->fields.front() <<" ..."<< std::endl;
#endif

      beit->cont.erase();
      if (status != 2)
        if (!this->getFields(is,*beit))
        {
          entry.name = beit->name;
          entry.fields = beit->fields;
          STOPP_TIMER("getNextEntry")
          return false;
        }

      if (beit->cont.empty())
      {
        // This entry is now completed, process it and delete it from the list
        if (!this->processThisEntry(*beit)) procOK = false;
        ucEntries.erase(beit);
      }

      // Start over again reading a new entry
      STOPP_TIMER("getNextEntry")
      return this->getNextEntry(is,entry);
    }

  // We found a new bulk entry, now read the data fields, if any
  last_entry = field;
  entry.name = field;
  entry.cont.erase();
  entry.fields.clear();
  bool ok = status == 2 ? true : this->getFields(is,entry);
  STOPP_TIMER("getNextEntry")
  return ok;
}


bool FFlNastranReader::getFields (std::istream& is, BulkEntry& entry)
{
  std::string field;
  int status = 1;
  int nFields = (entry.ffmt == LARGE_FIELD ? 5 : 9);

  // Read through the data line
  for (int i = 1; i < nFields; i++)
  {
    status = this->getNextField(is,field,entry.ffmt);
    if (status == 0) return false;
    if (status == 1 || !field.empty()) entry.fields.push_back(field);
    if (status == 2)
    {
      // We have reached end-of-line before reading the expected number of bytes
      // If the first character on the next line is a space, '+' or '*',
      // that line is still considered as a continuation of the present line
      char c;
      if (is.get(c) && !is.eof())
      {
        is.putback(c);
        if (strchr(" +*",c) || (c == ',' && entry.ffmt == FREE_FIELD))
        {
          // Yes, the next line is a continuation, but first fill up the
          // the remaining fields of the current line, if any, with ""
          if (field.empty()) entry.fields.push_back("");
          for (int j = i+1; j < nFields; j++) entry.fields.push_back("");
          entry.cont = "+";
        }
      }
      return true;
    }
  }

  // Special treatment of non-standard Free Field entries with 9 fields or more.
  // Continue to read until a continuation marker or EOL is detected.
  if (entry.ffmt == FREE_FIELD)
  {
    while ((status = this->getNextField(is,field,entry.ffmt)) == 1)
      entry.fields.push_back(field);
    if (status == 2)
    {
      // We have reached end-of-line before reading the expected number of bytes
      // If the first character on the next line is a space, ',', '+' or '*',
      // that line is still considered as a continuation of the present line
      char c;
      if (is.get(c) && !is.eof())
      {
	is.putback(c);
	if (!field.empty() && !strchr(" ,+*",c))
	{
	  // Not a continuation field, but a regular data field ends this line
	  entry.fields.push_back(field);
	  field = "";
	}
      }
    }
    else
      return false;
  }

  // Check for continuation field
  else if (!this->getNextField(is,field,CONT_FIELD))
    return false;

  // Some large field entries may use just an asterix as a continuation marker,
  // but this is also used to identify large field lines, we replace it with +
  if (entry.ffmt == LARGE_FIELD && field == "*")
    entry.cont = "+";

  else if (!field.empty())
    entry.cont = field;

  // Free field entries may continue also if they end with just a ','
  else if (entry.ffmt == FREE_FIELD)
  {
    char c;
    if (is.get(c) && !is.eof())
    {
      // Check that the next line actually is a continuation line
      if (c == ' ' || c == ',') entry.cont = "+";
      is.putback(c);
    }
  }

  // Small field entries may continue even without a continuation marker
  else if (entry.ffmt == SMALL_FIELD)
  {
    char c;
    if (is.get(c) && !is.eof())
    {
      // Check that the next line actually is a continuation line
      if (c == '+') entry.cont = "+";
      is.putback(c);
    }
  }

#if FFL_DEBUG > 4
  if (!entry.cont.empty())
    std::cout <<"FFlNastranReader: continuation field=\""<< entry.cont
	      <<"\""<< std::endl;
#endif

  return true;
}


int FFlNastranReader::getNextField (std::istream& is, std::string& field,
				    const FieldFormat size)

{
  // Return value: 0 - error, 1 - ok, 2 - ok, but end-of-line detected

  START_TIMER("getNextField")

  int nChar, retval = 1;
  size_t blank = 0;
  bool goOn = true;
  char c;

  field.erase();
  for (nChar = 1; goOn; nChar++)
    if (!is.get(c) || is.eof())
    {
      if (size == UNDEFINED)
      {
	// End-of-file encountered while searching for the next entry
	field = endOfBulk;
      }
      else
      {
	// End-of-file encountered while searching for a data field
	ListUI <<"\n *** Error: Premature end-of-file encountered."
	       <<" Nastran bulk-data is corrupt.\n";
 	retval = 0;
      }
      break;
    }
    else if (c == '\n' || c == '\r')
    {
      // End-of-line encountered
      lineCounter++;
      goOn = false;
      retval = 2;
#if FFL_DEBUG > 4
      std::cout <<"c=EOL\n";
#endif
      if (c == '\r')
      {
	// In some files a carriage-return is always followed by a line-feed
	if (is.get(c) && !is.eof())
	  if (c != '\n') is.putback(c);
      }
      if (size == CONT_FIELD && field.empty())
      {
	// We were looking for a continuation marker but just found newline
	// If the first character on the next line is a space, '+' or '*',
	// that line is still considered as a continuation of the present line
	if (is.get(c) && !is.eof())
	{
	  if (strchr(" +*",c)) field = "+";
	  is.putback(c);
	}
      }
    }
    else if (c == '$')
    {
      // This is a comment line, ignore it
      char line[BUFSIZ];
      is.getline(line,BUFSIZ,'\n');
      if (strlen(line) > 1)
      {
	// Store the comments for extraction of attribute names, etc.
	lastComment.first = lineCounter;
	lastComment.second += c;
	lastComment.second += line;
	lastComment.second += '\n';
      }
      lineCounter++;
      goOn = false;
      retval = 2;
#if FFL_DEBUG > 4
      std::cout <<"c='$' "<< nChar <<" comment line (ignored)\n";
#endif
    }
    else if (size == FREE_FIELD)
    {
      // Free field format, read everything until next ',' or ' ' or '\t'
      if (c != ',' && c != ' ' && c != '\t')
	field += c;
      else if (c == ',' || !field.empty())
	goOn = false;
#if FFL_DEBUG > 4
      std::cout <<"c='"<< c <<"' "<< nChar <<" free field\n";
#endif
    }
    else if (c == '\t')
    {
      // Tab-character encountered, currently disallowed
      ListUI <<"\n *** Error: Tabulators are not allowed in a bulk data file.\n"
             <<"            Replace them by space characters and try again.\n"
             <<"            Line: "<< lineCounter+1 <<"\n";
      if (field.length() > 0)
	ListUI <<"            Field: \""<< field
	       << std::string(blank,' ') << c <<"\"\n"
	       << std::string(20+field.length()+blank,' ') <<"^\n";
      retval = 0;
      break;
    }
    else if (size > 0)
    {
      // Small or large field format, read exactly size characters
      if (nChar == size)
      {
	goOn = false;
	if (c != ' ') field += c;
#if FFL_DEBUG > 4
	std::cout <<"c='"<< c <<"' "<< nChar <<" field finished\n";
#endif
      }
      else if (c != ' ')
      {
	if (blank)
	{
	  ListUI <<"\n *** Error: Embedded blanks are not allowed.\n"
		 <<"            Line: "<< lineCounter+1 <<"\n"
		 <<"            Field: \""<< field
		 << std::string(blank,' ') << c <<"\"\n"
		 << std::string(20+field.length(),' ')
		 << std::string(blank,'^') <<"\n";
	  retval = 0;
	  break;
	}
 	field += c;
#if FFL_DEBUG > 4
	std::cout <<"c='"<< c <<"' "<< nChar <<" field"<< size <<"\n";
#endif
      }
      else if (field.length() > 0)
	blank++;
    }
    else
    {
      // No field format yet, we are reading a new entry-name field
      if (nChar == SMALL_FIELD && size != CONT_FIELD) goOn = false;
      if (c == ',') goOn = false;
      if (c != ' ') field += c;
#if FFL_DEBUG > 4
      std::cout <<"c='"<< c <<"' "<< nChar
		<< (size==CONT_FIELD ? " continuation field\n":" field name\n");
#endif
      if (nChar == (int)strlen(endOfBulk))
	if (field == endOfBulk) break;
    }

#if FFL_DEBUG > 3
  std::cout <<"FFlNastranReader: field=\""<< field
	    <<"\" retval="<< retval << std::endl;
#endif
#ifdef FFL_DEBUG
  if (lineCounter%1000 == 0 && retval == 2)
    std::cout <<"FFlNastranReader: processed "<< lineCounter <<" lines"
	      << std::endl;
#endif
  STOPP_TIMER("getNextField")
  return retval;
}


bool FFlNastranReader::processThisEntry (BulkEntry& entry)
{
  START_TIMER("processThisEntry")
#if FFL_DEBUG > 2
  std::cout <<"Entry: \""<< entry.name <<"\"\nField:";
  for (const std::string& fld : entry.fields) std::cout <<" \""<< fld <<"\"";
  std::cout << std::endl;
#endif
  bool ok = this->processThisEntry (entry.name,entry.fields);
  if (sizeOK) entry.fields.clear();
  lastComment = std::make_pair(0,"");
  STOPP_TIMER("processThisEntry")
  return ok;
}


bool FFlNastranReader::processSet (const std::string& fname, const int startBulk)
{
  START_TIMER("processSet")

  std::ifstream fs(fname.c_str());
  if (!fs) return false;

#ifdef FFL_DEBUG
  std::cout <<"FFlNastranReader: starting SET parsing within lines 1-"
	    << startBulk-1 << std::endl;
#endif

  char line[256];
  int nError = 0;
  FFlGroup* aGroup = NULL;
  for (int lCounter = 0; !fs.eof() && lCounter < startBulk; lCounter++)
  {
    fs.getline(line,255);
    // Check for comment line possibly containing a group name
    if (line[0] == '$' && strlen(line) > 8)
      lastComment = std::make_pair(lCounter+1,line);
    // Check for the "SET" keyword
    else if (strncmp(line,setIdent,strlen(setIdent)) == 0)
    {
      // Delete trailing whitespaces, if any
      int startLin = lCounter+1;
      int lastChar = strlen(line);
      while (lastChar > 0 && isspace(line[--lastChar]));

      std::string setLine(line,lastChar+1);
      while (setLine[setLine.length()-1] == ',' && !fs.eof())
      {
        // The SET definition is continued over several lines,
        // read the entire SET definition into one long string
        fs.getline(line,255);
        // Delete trailing whitespaces, if any
        lastChar = strlen(line);
        while (lastChar > 0 && isspace(line[--lastChar]));
        setLine = setLine + std::string(line,lastChar+1);
        lCounter++;
      }

      if (lastComment.first > 0)
        if (lastComment.second.substr(0,18) == "$*  Group (nodes):")
          continue; // Ignore node groups from NX

      // Process the set definition and create an associated FFlGroup
      aGroup = this->processThisSet(setLine,startLin,lCounter+1);
      if (aGroup)
      {
        aGroup->sortElements(true);
        if (!myLink->addGroup(aGroup))
          aGroup = NULL;
        // Check if the group is named
        else if (lastComment.first > 0 && lastComment.first < startLin)
          if (extractNameFromComment(lastComment.second))
          {
            lastComment.first = 0;
            aGroup->setName(lastComment.second);
          }
      }
      else
        nError++;
    }
    else if (isBulkData(line))
      break;
    else
      aGroup = NULL;

    // Check if the last read comment contains a name of the last created group
    if (aGroup && lastComment.first > 0)
    {
      lastComment.first = 0;
      if (extractNameFromComment(lastComment.second))
        aGroup->setName(lastComment.second);
    }
  }

#ifdef FFL_DEBUG
  std::cout <<"FFlNastranReader: SET parsing done."<< std::endl;
#endif
  STOPP_TIMER("processSet")
  return nError == 0;
}


FFlGroup* FFlNastranReader::processThisSet (std::string& setLine,
#ifdef FFL_DEBUG
					    const int startL, const int stopL
#else
					    const int, const int
#endif
){
  START_TIMER("processThisSet")

  // Create an FFlGroup with the corresponding SET Id
  char* endPtr;
  int setId = strtol(setLine.c_str()+3,&endPtr,10);
  FFlGroup* aGroup = new FFlGroup(setId,"Nastran SET");
#ifdef FFL_DEBUG
  std::cout <<"FFlNastranReader: Processing SET "<< setId;
  if (startL == stopL)
    std::cout <<" line "<< startL << std::endl;
  else if (startL == stopL-1)
    std::cout <<" lines "<< startL <<", "<< stopL << std::endl;
  else
    std::cout <<" lines "<< startL <<" - "<< stopL << std::endl;
#endif

  int oldNotes = nNotes;

  // Process the SET definition
  int lastNumberAdded = 0;
  const char* p = setLine.c_str() + setLine.find_first_of('=',1) + 1;
  char* charPtr = strtok((char*)p," ,");
  while (charPtr)
  {
    if (*charPtr == 'A') // "ALL" keyword
    {
      ElementsCIter eit;
      for (eit = myLink->elementsBegin(); eit != myLink->elementsEnd(); ++eit)
        aGroup->addElement((*eit)->getID());
      STOPP_TIMER("processThisSet")
      return aGroup;
    }
    else if (*charPtr == 'T') // "THRU" keyword
    {
      charPtr  = strtok(NULL," ,");
      int low  = lastNumberAdded+1;
      int high = charPtr ? strtol(charPtr,&endPtr,10) : 0;
      if (high == 0)
      {
        ListUI <<"\n *** Syntax Error in Nastran SET definition: \""
               << setLine <<" ...\" (ignored)\n";
        delete aGroup;
        STOPP_TIMER("processThisSet")
        return NULL;
      }
      charPtr = strtok(NULL," ,");

      int excludedFromSet = -999;
      if (charPtr)
        if (*charPtr == 'E') // "EXCEPT" keyword
        {
          charPtr = strtok(NULL," ,");
          excludedFromSet = charPtr ? strtol(charPtr,&endPtr,10) : 0;
          if (excludedFromSet == 0)
          {
            ListUI <<"\n *** Syntax Error in Nastran SET definition: \""
                   << setLine <<" ...\" (ignored)\n";
            delete aGroup;
            STOPP_TIMER("processThisSet")
            return NULL;
          }
          charPtr = strtok(NULL," ,");
        }

      for (int j = low; j <= high; j++)
        if (j == excludedFromSet)
        {
          if (!charPtr)
            excludedFromSet = -999;
          else
          {
            excludedFromSet = strtol(charPtr,&endPtr,10);
            if (excludedFromSet == 0)
            {
              ListUI <<"\n *** Syntax Error in Nastran SET definition: \""
                     << setLine <<" ...\" (ignored)\n";
              delete aGroup;
              STOPP_TIMER("processThisSet")
              return NULL;
            }
            charPtr = strtok(NULL," ,");
          }
        }
        else if (myLink->getElement(j))
          aGroup->addElement(j);
        else if (nNotes++ < oldNotes+10)
          ListUI <<"\n   * Note: Ignoring non-existing element "<< j
                 <<" in Nastran SET "<< setId <<"\n";

      if (excludedFromSet > 0)
      {
        lastNumberAdded = excludedFromSet;
        if (myLink->getElement(excludedFromSet))
          aGroup->addElement(excludedFromSet);
        else if (nNotes++ < oldNotes+10)
          ListUI <<"\n   * Note: Ignoring non-existing element "
                 << excludedFromSet <<" in Nastran SET "<< setId <<"\n";
      }
    }
    else
    {
      int currNumb = strtol(charPtr,&endPtr,10);
      if (currNumb == 0)
      {
        ListUI <<"\n *** Syntax Error in Nastran SET definition: \""
               << setLine <<" ...\" (ignored)\n";
        delete aGroup;
        STOPP_TIMER("processThisSet")
        return NULL;
      }
      charPtr = strtok(NULL," ,");
      lastNumberAdded = currNumb;
      if (myLink->getElement(currNumb))
        aGroup->addElement(currNumb);
      else if (nNotes++ < oldNotes+10)
        ListUI <<"\n   * Note: Ignoring non-existing element "<< currNumb
               <<" in Nastran SET "<< setId <<"\n";
    }
  }

  if (nNotes > oldNotes+10)
  {
    nWarnings++;
    ListUI <<"\n  ** Warning: "<< nNotes-oldNotes
	   <<" non-existing elements were detected for Nastran SET "<< setId
	   <<".\n              Only the 10 first are reported."
	   <<"\n              Please verify that the model is consistent.\n";
  }

  STOPP_TIMER("processThisSet")
  return aGroup;
}


void FFlNastranReader::processAssignFile (const std::string& line)
{
  size_t i = line.find("output2");
  if (i == std::string::npos) i = line.find("OUTPUT2");
  if (i == std::string::npos) return;

  size_t j = line.find('\'',i+7);
  if (j == std::string::npos) return;
  size_t k = line.find('\'',j+1);
  if (k == std::string::npos || k < j+2) return;

  std::string op2file = line.substr(j+1,k-j-1);
  FFaFilePath::makeItAbsolute(op2file,mainPath);
  myLink->addOP2file(FFaFilePath::checkName(op2file));

  nNotes++;
  ListUI <<"\n   * Note: OP2 file detected: "<< op2file <<"\n";
}


bool FFlNastranReader::extractNameFromComment (std::string& commentLine,
					       bool first)
{
#if FFL_DEBUG > 1
  std::cout <<"FFlNastranReader: Processing comment\n"<< commentLine;
#endif

  // Check for I-DEAS syntax
  size_t pos = first ? commentLine.find("name: ") : commentLine.rfind("name: ");
  if (pos < commentLine.size())
  {
    // Erase everything before "name: "
    commentLine.erase(0,pos+6);
    pos = commentLine.find_first_of("\r\n");
    if (pos < commentLine.size()) commentLine.erase(pos);
#if FFL_DEBUG > 1
    if (!commentLine.empty())
      std::cout <<"\tFound name: "<< commentLine << std::endl;
#endif
    return !commentLine.empty();
  }

  // Check for NX syntax
  pos = first ? commentLine.find("$*  NX ") : commentLine.rfind("$*  NX ");
  if (pos < commentLine.size())
  {
    size_t pos2 = commentLine.substr(pos).find(": ");
    if (pos2 < commentLine.size()-pos)
    {
      // Erase everything before ": "
      commentLine.erase(0,pos+pos2+2);
      pos = commentLine.find_first_of("\r\n");
      if (pos < commentLine.size()) commentLine.erase(pos);
#if FFL_DEBUG > 1
      if (!commentLine.empty())
	std::cout <<"\tFound name: "<< commentLine << std::endl;
#endif
      return !commentLine.empty();
    }
  }

  // Check for new NX syntax for element groups
  pos = first ? std::string::npos : commentLine.rfind("$*  Group (elements):");
  if (pos < commentLine.size())
  {
    size_t pos2 = commentLine.substr(pos).rfind(": ");
    if (pos2 < commentLine.size()-pos)
    {
      // Erase everything before the last ": "
      commentLine.erase(0,pos+pos2+2);
      pos = commentLine.find_first_of("\r\n");
      if (pos < commentLine.size()) commentLine.erase(pos);
#if FFL_DEBUG > 1
      if (!commentLine.empty())
        std::cout <<"\tFound name: "<< commentLine << std::endl;
#endif
      return !commentLine.empty();
    }
  }

  // Check for HyperMesh syntax
  pos = first ? commentLine.find("$HMNAME") : commentLine.rfind("$HMNAME");
  if (pos == std::string::npos)
    pos = first ? commentLine.find("$HMSET") : commentLine.rfind("$HMSET");
  if (pos < commentLine.size())
  {
    // Erase everything outside the last ""-pair
    commentLine.erase(0,commentLine.find('"',pos)+1);
    commentLine.erase(commentLine.find('"'));
#if FFL_DEBUG > 1
    if (!commentLine.empty())
      std::cout <<"\tFound name: "<< commentLine << std::endl;
#endif
    return !commentLine.empty();
  }

  return false;
}
