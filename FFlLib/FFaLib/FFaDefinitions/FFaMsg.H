// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFA_MSG_H
#define FFA_MSG_H

#include <string>
#include <vector>
#include <stack>
#include <map>

// The symbol ERROR might have already been defined if, e.g., windows.h
// is included before this header file. To avoid compiler problems due to name
// conflict with the FFaDialogType enum value below, it must be undefined here.
// TODO: Rename FFaMsg::ERROR to something else?
#ifdef ERROR
#undef ERROR
#endif

//! \brief Output stream for feedback messages to the user.
//! \details This symbol works like an output stream.
//! Use this instead of std::cout and std::cerr when the message is meaningful
//! to the user too (and not only for the developers).
#define ListUI FFaMsg::getMessager()

class FaVec3;


/*!
  \brief Interface class for various user feedback.
*/

class FFaMsg
{
public:
  //! \brief Default constructor.
  FFaMsg() {}
  //! \brief Empty destructor.
  virtual ~FFaMsg() {}

  //! \brief Sets the real messager object to be used.
  static void setMessager(FFaMsg* messager = NULL);
  //! \brief Returns the current messager object.
  static FFaMsg& getMessager();
  //! \brief Returns \e true if we are using the default messager.
  static bool usingDefault() { return ourCurrentMessager == NULL; }

  //! \brief Enums defining various dialog box types.
  enum FFaDialogType
  {
    DISMISS_WARNING,
    DISMISS_ERROR,
    DISMISS_INFO,
    WARNING,
    ERROR,
    OK,
    OK_CANCEL,
    YES_NO,
    YES_NO_CANCEL,
    _ALL_,
    OK_ALL_CANCEL,
    YES_ALL_NO,
    YES_ALL_NO_ALL,
    YES_ALL_NO_CANCEL,
    YES_ALL_NO_ALL_CANCEL,
    YES_NO_ALL,
    YES_NO_ALL_CANCEL,
    GENERIC
  };

  // User feedback methods. These will be routed
  // to the messager provided, or to the default one.

  static void resetToAllAnswer(const FFaDialogType type = _ALL_);

  static int dialog(const std::string& message,
		    const FFaDialogType type = DISMISS_WARNING,
		    const char** genericButtons = (const char**)NULL);

  static int dialog(const std::string& message,
		    const FFaDialogType type,
		    const std::vector<std::string>& buttonTexts);

  static int dialog(int& returnedSelectionIdx, const std::string& message,
		    const FFaDialogType type,
		    const std::vector<std::string>& buttonTexts,
		    const std::vector<std::string>& selectionList);

  static void list(const std::string& str, bool onScreen = false);
  static void tip (const std::string& str);

  static void setStatus   (const std::string& statusText);
  static void changeStatus(const std::string& statusText);
  static void pushStatus  (const std::string& statusText);
  static void popStatus   ();

  static void enableSubSteps (int steps);
  static void setSubStep     (int step);
  static void disableSubSteps();

  static void displayTime(int hour, int min, int sec);
  static void clearTime  ();

  static void setSubTask(const std::string& taskText);

  static void enableProgress (int nSteps);
  static void setProgress    (int progressStep);
  static void disableProgress();

  // Stream-like operators for easy printing of messages to output list
  friend FFaMsg& operator<< (FFaMsg&, const std::string&);
  friend FFaMsg& operator<< (FFaMsg&, const char*);
  friend FFaMsg& operator<< (FFaMsg&, const FaVec3&);
  friend FFaMsg& operator<< (FFaMsg&, double);
  friend FFaMsg& operator<< (FFaMsg&, char);
  friend FFaMsg& operator<< (FFaMsg&, int);
  friend FFaMsg& operator<< (FFaMsg&, long long int);
  friend FFaMsg& operator<< (FFaMsg&, size_t);

  //! \brief Opens a list file containing a copy of the output list window.
  virtual bool openListFile(const char* = NULL, bool = false) { return true; }
  //! \brief Closes the list file.
  virtual void closeListFile() {}

protected:

  // Virtual methods reimplemented to make real GUI meaning

  virtual int dialogVt(const std::string& message,
		       const FFaDialogType type,
		       const char** genericButtons);
  virtual int dialogVt(const std::string& message,
		       const FFaDialogType type,
		       const std::vector<std::string>& buttonTexts);
  virtual int dialogVt(int& returnedSelectionIdx,
		       const std::string& message,
		       const FFaDialogType type,
		       const std::vector<std::string>& buttonTexts,
		       const std::vector<std::string>& selectionList);

  virtual void listVt(const std::string& str, bool onScreen = false);
  virtual void tipVt (const std::string& str);

  virtual void setStatusVt   (const std::string&) {}
  virtual void changeStatusVt(const std::string&) {}
  virtual void pushStatusVt  (const std::string&) {}
  virtual void popStatusVt   () {}

  static void showStatus();

  virtual void enableSubStepsVt (int) {}
  virtual void setSubStepVt     (int) {}
  virtual void disableSubStepsVt() {}

  virtual void displayTimeVt(int, int, int) {}
  virtual void clearTimeVt  () {}

  virtual void setSubTaskVt(const std::string&) {}

  virtual void enableProgressVt (int) {}
  virtual void setProgressVt    (int) {}
  virtual void disableProgressVt() {}

  static std::stack<std::string> ourStatuses;

private:
  static std::map<FFaDialogType,int> ourToAllAnswer;

  static FFaMsg* ourCurrentMessager;
};

#endif
