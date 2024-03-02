// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include <cstdio>

#include "FFaLib/FFaDefinitions/FFaMsg.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"

std::stack<std::string>             FFaMsg::ourStatuses;
std::map<FFaMsg::FFaDialogType,int> FFaMsg::ourToAllAnswer;
FFaMsg*                             FFaMsg::ourCurrentMessager = NULL;


void FFaMsg::setMessager(FFaMsg* messager)
{
  if (ourCurrentMessager)
    delete ourCurrentMessager;
  ourCurrentMessager = messager;
}

FFaMsg& FFaMsg::getMessager()
{
  if (ourCurrentMessager)
    return *ourCurrentMessager;

  static FFaMsg defaultMessager;
  return defaultMessager;
}


/*!
  \brief Shows a dialog according to the dialog type.
  \return For the OK_CANCEL dialog: 1 for "Ok" and 0 for "Cancel".
  \return For the YES_NO dialog: 1 for "Yes" and 0 for "No".
  \return For the YES_NO_CANCEL dialog: 1, 0 and 2, respectively.
  \return For a GENERIC dialog: The index of the pushed button (maximum 3).
  \return For all other dialogs: Always 0.
*/

int FFaMsg::dialog(const std::string& statusText, const FFaDialogType type,
		   const char** genericButtons)
{
  if (type < _ALL_ || type == GENERIC)
    return getMessager().dialogVt(statusText,type,genericButtons);

  // We have an "... to all" dialog, check if the answer is aready given
  std::map<FFaDialogType,int>::iterator ait = ourToAllAnswer.find(type);
  if (ait != ourToAllAnswer.end())
    return ait->second; // return the previous answer

  // Set up a new dialog with "to all" options.
  // Note that "Cancel" always means "Cancel all".

  int button = -1;
  std::vector<std::string> buttonTexts;
  switch (type)
    {
    case OK_ALL_CANCEL:
      buttonTexts.push_back("Ok");
      buttonTexts.push_back("Ok to all");
      buttonTexts.push_back("Cancel");
      button = getMessager().dialogVt(statusText,type,buttonTexts);
      if (button >= 1)
	button = ourToAllAnswer[type] = button-1;
      break;

    case YES_ALL_NO:
    case YES_ALL_NO_ALL:
    case YES_ALL_NO_CANCEL:
      buttonTexts.push_back("Yes");
      buttonTexts.push_back("Yes to all");
      buttonTexts.push_back("No");
      if (type == YES_ALL_NO_ALL)
	buttonTexts.push_back("No to all");
      else if (type == YES_ALL_NO_CANCEL)
	buttonTexts.push_back("Cancel");
      button = getMessager().dialogVt(statusText,type,buttonTexts);
      if (button == 1)
	button = ourToAllAnswer[type] = 0;
      else if (button == 2)
	button = 1;
      else if (button == 3)
	button = ourToAllAnswer[type] = type == YES_ALL_NO_ALL ? 1 : 2;
      break;

    case YES_ALL_NO_ALL_CANCEL:
      buttonTexts.push_back("Yes");
      buttonTexts.push_back("Yes to all");
      buttonTexts.push_back("No");
      buttonTexts.push_back("No to all");
      buttonTexts.push_back("Cancel");
      button = getMessager().dialogVt(statusText,type,buttonTexts);
      if (button == 1)
	button = ourToAllAnswer[type] = 0;
      else if (button == 2)
	button = 1;
      else if (button >= 3)
	button = ourToAllAnswer[type] = button-2;
      break;

    case YES_NO_ALL:
    case YES_NO_ALL_CANCEL:
      buttonTexts.push_back("Yes");
      buttonTexts.push_back("No");
      buttonTexts.push_back("No to all");
      if (type == YES_NO_ALL_CANCEL)
	buttonTexts.push_back("Cancel");
      button = getMessager().dialogVt(statusText,type,buttonTexts);
      if (button >= 2)
	button = ourToAllAnswer[type] = button-1;
      break;

    default:
      break;
    }

  if (button == 0 || button == 1)
  {
    // Swap the Yes/No reply values (such that Yes/Ok=1 and No=0).
    // This is to be consistent with the Yes/No dialog.
    button = 1-button;
    ait = ourToAllAnswer.find(type);
    if (ait != ourToAllAnswer.end())
      ait->second = 1 - ait->second;
  }

  return button;
}


void FFaMsg::resetToAllAnswer(const FFaDialogType type)
{
  if (type == _ALL_)
    ourToAllAnswer.clear();
  else
    ourToAllAnswer.erase(type);
}


int FFaMsg::dialog(const std::string& message,
		   const FFaDialogType type,
		   const std::vector<std::string>& buttonTexts)
{
  return getMessager().dialogVt(message, type, buttonTexts);
}

int FFaMsg::dialog(int& returnedSelectionIdx,
		   const std::string& message,
		   const FFaDialogType type,
		   const std::vector<std::string>& buttonTexts,
		   const std::vector<std::string>& selectionList)
{
  return getMessager().dialogVt(returnedSelectionIdx, message,
				type, buttonTexts, selectionList);
}


void FFaMsg::list(const std::string& statusText, bool onScreen)
{
  getMessager().listVt(statusText,onScreen);
}


void FFaMsg::tip(const std::string& statusText)
{
  getMessager().tipVt(statusText);
}


void FFaMsg::setStatus(const std::string& statusText)
{
  getMessager().setStatusVt(statusText);
}


void FFaMsg::changeStatus(const std::string& statusText)
{
  getMessager().changeStatusVt(statusText);
}


void FFaMsg::pushStatus(const std::string& statusText)
{
  getMessager().pushStatusVt(statusText);
}

void FFaMsg::popStatus()
{
  getMessager().popStatusVt();
}


void FFaMsg::enableSubSteps(int steps)
{
  getMessager().enableSubStepsVt(steps);
}

void FFaMsg::setSubStep(int step)
{
  getMessager().setSubStepVt(step);
}

void FFaMsg::disableSubSteps()
{
  getMessager().disableSubStepsVt();
}


void FFaMsg::displayTime(int hour, int min, int sec)
{
  getMessager().displayTimeVt(hour,min,sec);
}

void FFaMsg::clearTime()
{
  getMessager().clearTimeVt();
}

void FFaMsg::setSubTask(const std::string& taskText)
{
  getMessager().setSubTaskVt(taskText);
}


void FFaMsg::enableProgress(int nSteps)
{
  getMessager().enableProgressVt(nSteps);
}

void FFaMsg::setProgress(int progressStep)
{
  getMessager().setProgressVt(progressStep);
}

void FFaMsg::disableProgress()
{
  getMessager().disableProgressVt();
}


FFaMsg& operator<< (FFaMsg& messager, const std::string& str)
{
  messager.listVt(str);
  return messager;
}

FFaMsg& operator<< (FFaMsg& messager, const char* str)
{
  if (str) messager.listVt(str);
  return messager;
}

FFaMsg& operator<< (FFaMsg& messager, const FaVec3& vec)
{
  char str[64];
  snprintf(str,64,"%g %g %g",vec.x(),vec.y(),vec.z());
  messager.listVt(str);
  return messager;
}

FFaMsg& operator<< (FFaMsg& messager, double d)
{
  char str[16];
  snprintf(str,16,"%g",d);
  messager.listVt(str);
  return messager;
}

FFaMsg& operator<< (FFaMsg& messager, char c)
{
  char str[2] = { c, '\0' };
  messager.listVt(str);
  return messager;
}

FFaMsg& operator<< (FFaMsg& messager, int i)
{
  messager.listVt(std::to_string(i));
  return messager;
}

FFaMsg& operator<< (FFaMsg& messager, long long int i)
{
  messager.listVt(std::to_string(i));
  return messager;
}

FFaMsg& operator<< (FFaMsg& messager, size_t i)
{
  messager.listVt(std::to_string(i));
  return messager;
}


// Default implementations of the virtual methods.
// These are used when the program is run in terminal mode, without GUI.

int FFaMsg::dialogVt(const std::string& msg, const FFaDialogType, const char**)
{
  std::cout << msg << std::endl;
  return -1;
}


int FFaMsg::dialogVt(const std::string& message, const FFaDialogType,
		     const std::vector<std::string>& buttonTexts)
{
  std::cout << message << std::endl;
  for (size_t i = 0; i < buttonTexts.size(); i++)
    std::cout <<"["<< buttonTexts[i] <<"] = " << i <<" "<< std::flush;

  int answer = -1;
  std::cin >> answer;
  return answer;
}


int FFaMsg::dialogVt(int& returnedSelectionIdx,
		     const std::string& message,
		     const FFaDialogType,
		     const std::vector<std::string>& buttonTexts,
		     const std::vector<std::string>& selectionList)
{
  std::cout << message << std::endl;
  for (size_t i = 0; i < selectionList.size(); i++)
    std::cout <<"["<< selectionList[i] <<"] = "<< i <<" "<< std::flush;

  returnedSelectionIdx = -1;
  if (selectionList.size())
    std::cin >> returnedSelectionIdx;

  for (size_t j = 0; j < buttonTexts.size(); j++)
    std::cout <<"["<< buttonTexts[j] <<"] = "<< j <<" "<< std::flush;

  int answer = -1;
  std::cin >> answer;
  return answer;
}


void FFaMsg::listVt(const std::string& str, bool)
{
  std::cout << str << std::flush;
}


void FFaMsg::tipVt(const std::string& str)
{
  std::cout <<"Tip: "<< str << std::endl;
}


void FFaMsg::showStatus()
{
  if (ourStatuses.empty()) return;

  std::cout <<"Status "<< ourStatuses.size() <<": ";
  for (size_t i = 1; i < ourStatuses.size(); i++) std::cout <<" ";
  std::cout << ourStatuses.top() << std::endl;
}
