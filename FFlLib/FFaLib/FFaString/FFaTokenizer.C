// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFaTokenizer.C
  \brief One-level token hierarchies.
*/

#include <iterator>

#include "FFaLib/FFaString/FFaTokenizer.H"


/*!
  \brief Base class for generic token input.
*/

class FFaTokenInput
{
protected:
  //! \brief Default constructor.
  FFaTokenInput() {}

public:
  //! \brief Empty destructor.
  virtual ~FFaTokenInput() {}

  //! \brief Initializes the instance, returning the first character to process.
  virtual void init(int&) {}
  //! \brief Checks for end-of-file.
  virtual bool eof() const = 0;
  //! \brief Returns the next character to process.
  virtual int get() = 0;
};


/*!
  \brief Class for token input from file.
*/

class FFaFileData : public FFaTokenInput
{
  FILE* fd; //!< The file to read from

public:
  //! \brief The constructor initializes the file pointer.
  FFaFileData(FILE* f) : fd(f) {}
  //! \brief Empty destructor.
  virtual ~FFaFileData() {}

  //! \brief Checks for end-of-file.
  virtual bool eof() const { return feof(fd); }
  //! \brief Returns the next character to process.
  virtual int get() { return getc(fd); }
};


/*!
  \brief Class for token input from an input stream.
*/

class FFaStreamData : public FFaTokenInput
{
  std::istream& is; //!< The input stream to read from

public:
  //! \brief The constructor initializes the input stream reference.
  FFaStreamData(std::istream& s) : is(s) {}
  //! \brief Empty destructor.
  virtual ~FFaStreamData() {}

  //! \brief Checks for end-of-file.
  virtual bool eof() const { return is.eof(); }
  //! \brief Returns the next character to process.
  virtual int get() { return is.get(); }
};


/*!
  \brief Class for token input from a string range.
*/

class FFaIteratorData : public FFaTokenInput
{
  std::string::const_iterator it;  //!< Running iterator
  std::string::const_iterator jt;  //!< Iterator pointing to next character
  std::string::const_iterator end; //!< End iterator

public:
  //! \brief The constructor initializes the iterators.
  FFaIteratorData(std::string::const_iterator b,
                  std::string::const_iterator e) : it(b), jt(b), end(e) {}
  //! \brief Empty destructor.
  virtual ~FFaIteratorData() {}

  //! \brief Initializes the instance, returning the first character to process.
  virtual void init(int& c) { if (it != end) { c = this->get(); } }
  //! \brief Checks for end-of-file.
  virtual bool eof() const { return jt == end; }
  //! \brief Returns the next character to process.
  virtual int get() { jt = it++; return *jt; }
};


void FFaTokenizer::createTokens (FILE* tokenFile)
{
  FFaFileData pFile(tokenFile);
  this->createTokens(pFile);
}


void FFaTokenizer::createTokens (std::istream& tokenStream)
{
  FFaStreamData pStream(tokenStream);
  this->createTokens(pStream);
}


std::string::const_iterator
FFaTokenizer::createTokens (std::string::const_iterator tokenBegin,
			    std::string::const_iterator tokenEnd)
{
  FFaIteratorData pData(tokenBegin,tokenEnd);
  return tokenBegin + this->createTokens(pData);
}


int FFaTokenizer::createTokens (FFaTokenInput& tokenData)
{
  this->clear();
  this->reserve(10);

  std::string token;
  token.reserve(512);

  int  c = eb;
  int  counter = 0, subEntryCounter = 0;
  bool readingText = false;

  for (tokenData.init(c); !tokenData.eof(); c = tokenData.get(), counter++)
  {
    // Check for sub-entry begin and end, because we only want
    // to split "this" entry into tokens, and not the sub-entries
    if (!readingText && (c == eb || c == '<' || c == '[' || c == '{'))
      subEntryCounter++;

    // Check for quoted text and toggle if we shall
    // ignore other "keywords" than the text quote
    if (c == '"')
      readingText = !readingText;

    // Insert the char into the current token string,
    // unless equal to the token separator of "this" entry, or begin/end.
    // Skip if it is a quote, and we don't want it (we always want it
    // in a sub-entry, which we are not supposed to touch here).
    if ((!iAmStrippingQuotes || subEntryCounter != 1 || c != '"') && c > 0)
    {
      if (token.capacity() < token.size()+1)
	token.reserve(token.capacity()*3);

      if (readingText)
	token += (char)c;
      else if (subEntryCounter > 1 || (c != ts && c != eb && c != ee))
	if (!isspace(c)) token += (char)c; // Strip whitespace
    }

    // Check if end of field in "this" entry, emit the field text if it is
    if (!readingText && subEntryCounter == 1 && (c == ts || c == ee))
    {
      if (this->capacity() < this->size()+1)
	this->reserve(this->capacity()*3);
      this->push_back(std::string());
      token.swap(this->back());
    }

    // If hit a sub-entry end, decrement the sub-entry counter
    if (!readingText && (c == ee || c == '>' || c == ']' || c == '}'))
      subEntryCounter--;

    // Continue with next char in input if not finished
    if (subEntryCounter == 0) break;
  }

  return counter;
}
