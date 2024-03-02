// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFA_GENERIC_FACTORY_H
#define FFA_GENERIC_FACTORY_H

#include <string>
#include <vector>
#include <map>

#include "FFaLib/FFaPatterns/FFaSingelton.H"
#include "FFaLib/FFaDynCalls/FFaDynCB.H"


// TODO: Create doxygen documentation
// TODO: Create copy constructor creator
// TODO: Change class T to be a specific type, not a mix of pointer and type.
//       This will affect all users - beware...

template <class T, class Key = std::string, class ArgType = int>
class FFaGenericFactory : public FFaSingelton< FFaGenericFactory<T,Key,ArgType> >
{
  typedef std::map<Key, FFaDynCB2<ArgType,T*&> > CreatorType;
  typedef typename CreatorType::const_iterator   CreatorTypeCIter;

public:
  FFaGenericFactory() {}
  virtual ~FFaGenericFactory() {}

  bool registerCreator(const Key& key, const FFaDynCB2<ArgType,T*&>& creator)
  {
    return myCreatorMap.insert(std::make_pair(key,creator)).second;
  }

  T* create(const Key& key, ArgType id)
  {
    T* retVal = 0;
    typename CreatorType::iterator creator = myCreatorMap.find(key);
    if (creator != myCreatorMap.end())
      creator->second.invoke(id,retVal);
    return retVal;
  }

  void getKeys(std::vector<Key>& keys) const
  {
    keys.clear();
    keys.reserve(myCreatorMap.size());
    for (CreatorTypeCIter it = myCreatorMap.begin(); it != myCreatorMap.end(); ++it)
      keys.push_back(it->first);
  }

  void clear() { myCreatorMap.clear(); this->removeInstance(); }

private:
  CreatorType myCreatorMap;
};

#endif
