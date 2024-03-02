// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFaTokenizer.H
  \brief One-level token hierarchies.
*/

#ifndef FFA_TOKENIZER_H
#define FFA_TOKENIZER_H

#include <vector>
#include <string>
#include <cstdio>
#include <iostream>

class FFaTokenInput;


/*!
  \brief Class creating a one-level token hierarchy.
  \details The class is essentially a list of strings, with some added
  functionality for parsing the string elements (tokens) based on the specified
  start-, stop- and separator characters.
*/

class FFaTokenizer : public std::vector<std::string>
{
public:
  //! \brief The constructor defines the start, stop and separator characters.
  //! \param[in] entryBegin Character indicating the start of the token list
  //! \param[in] entryEnd Character indicating the end of the token list
  //! \param[in] separator Character separating each token
  FFaTokenizer(char entryBegin, char entryEnd, char separator = ',')
    : eb(entryBegin), ee(entryEnd), ts(separator), iAmStrippingQuotes(true) {}

  //! \brief Constructor parsing the tokens from a file.
  //! \param[in] tokenFile File to parse the tokens from
  //! \param[in] entryBegin Character indicating the start of the token list
  //! \param[in] entryEnd Character indicating the end of the token list
  //! \param[in] separator Character separating each token
  //! \param[in] sq If \e true, strip quotes ("-characters) from the tokens
  FFaTokenizer(FILE* tokenFile, char entryBegin, char entryEnd,
               char separator = ',', bool sq = true)
    : eb(entryBegin), ee(entryEnd), ts(separator), iAmStrippingQuotes(sq)
  {
    this->createTokens(tokenFile);
  }

  //! \brief Constructor parsing the tokens from an input stream.
  //! \param[in] tokenStream Input stream to parse the tokens from
  //! \param[in] entryBegin Character indicating the start of the token list
  //! \param[in] entryEnd Character indicating the end of the token list
  //! \param[in] separator Character separating each token
  //! \param[in] sq If \e true, strip quotes ("-characters) from the tokens
  FFaTokenizer(std::istream& tokenStream, char entryBegin, char entryEnd,
               char separator = ',', bool sq = true)
    : eb(entryBegin), ee(entryEnd), ts(separator), iAmStrippingQuotes(sq)
  {
    this->createTokens(tokenStream);
  }

  //! \brief Constructor parsing the tokens from a string.
  //! \param[in] tokenString String to parse the tokens from
  //! \param[in] entryBegin Character indicating the start of the token list
  //! \param[in] entryEnd Character indicating the end of the token list
  //! \param[in] separator Character separating each token
  //! \param[in] sq If \e true, strip quotes ("-characters) from the tokens
  FFaTokenizer(const std::string& tokenString, char entryBegin, char entryEnd,
               char separator = ',', bool sq = true)
    : eb(entryBegin), ee(entryEnd), ts(separator), iAmStrippingQuotes(sq)
  {
    this->createTokens(tokenString.begin(), tokenString.end());
  }

  //! \brief Empty destructor.
  virtual ~FFaTokenizer() {}

  //! \brief Creates tokens from a string range.
  std::string::const_iterator createTokens(std::string::const_iterator tBegin,
                                           std::string::const_iterator tEnd);

protected:
  //! \brief Creates tokens from a file.
  void createTokens(FILE* tokenFile);
  //! \brief Creates tokens from an input stream.
  void createTokens(std::istream& tokenStream);
  //! \brief Creates tokens from generic input.
  int createTokens(FFaTokenInput& tokenData);

private:
  char eb; //!< Start character
  char ee; //!< End character
  char ts; //!< Token separator

  bool iAmStrippingQuotes; //!< If \e true, strip "-characters from tokens
};

#endif
