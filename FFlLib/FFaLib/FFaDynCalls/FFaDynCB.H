// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFaDynCB.H
  \brief Macros for definition of dynamic callback objects.
  \details Enables callbacks to point at methods in specific objects,
  or to standard static methods or functions.
*/

#ifndef FFA_DYNCB_H
#define FFA_DYNCB_H

///////////////////////////////////////////////////////////////////////////////
//
// User interface:
//    Use these classes only through the following macros and methods:
//
//    To create a DynCB:
//

#define FFaDynCB0M(ownerType, ownerPt, methodName) \
          FFaDynCB0(new FFaMethodCB0< ownerType > (ownerPt, &ownerType::methodName))
#define FFaDynCB1M(ownerType, ownerPt, methodName, P1Type) \
          FFaDynCB1 < P1Type > (new FFaMethodCB1< ownerType , P1Type > (ownerPt, &ownerType::methodName))
#define FFaDynCB2M(ownerType, ownerPt, methodName, P1Type, P2Type) \
          FFaDynCB2 < P1Type, P2Type > (new FFaMethodCB2< ownerType , P1Type , P2Type > (ownerPt, &ownerType::methodName))
#define FFaDynCB3M(ownerType, ownerPt, methodName, P1Type, P2Type, P3Type)  \
          FFaDynCB3 < P1Type, P2Type, P3Type > (new FFaMethodCB3< ownerType , P1Type , P2Type , P3Type > (ownerPt, &ownerType::methodName))
#define FFaDynCB4M(ownerType, ownerPt, methodName, P1Type, P2Type, P3Type, P4Type)  \
          FFaDynCB4 < P1Type, P2Type, P3Type, P4Type > (new FFaMethodCB4< ownerType , P1Type , P2Type , P3Type, P4Type > (ownerPt, &ownerType::methodName))

#define FFaDynCB0S(function) \
          FFaDynCB0(new FFaStaticCB0 (function))
#define FFaDynCB1S(function, P1Type) \
          FFaDynCB1 < P1Type > (new FFaStaticCB1 < P1Type > (function))
#define FFaDynCB2S(function, P1Type, P2Type) \
          FFaDynCB2 < P1Type, P2Type > (new FFaStaticCB2 < P1Type , P2Type > (function))
#define FFaDynCB3S(function, P1Type, P2Type, P3Type) \
          FFaDynCB3 < P1Type, P2Type, P3Type > (new FFaStaticCB3 < P1Type , P2Type , P3Type > (function))
#define FFaDynCB4S(function, P1Type, P2Type, P3Type, P4Type) \
          FFaDynCB4 < P1Type, P2Type, P3Type, P4Type  > (new FFaStaticCB4 < P1Type , P2Type , P3Type, P4Type > (function))

//   To invoke the callback:
//      aDynCB->invoke(arg1,arg2);
//
//   To copy a call back:
//      FFaDynCB<N>* aCopyDynCB = aDynCB->copy();
//
//
// End user interface
///////////////////////////////////////////////////////////////////////////////


// one type :

#define FFADYNCB_H_TEMPLATE1            template< class P1Type >
#define FFADYNCB_H_ADDED_TEMPLATE_ARGS1 , class P1Type
#define FFADYNCB_H_TYPE_LIST1           < P1Type >
#define FFADYNCB_H_ADDED_TYPE_LIST1     , P1Type
#define FFADYNCB_H_ARGUMENT_LIST1       P1Type p1
#define FFADYNCB_H_PARAM_LIST1          p1


// Two types :

#define FFADYNCB_H_TEMPLATE2            template< class P1Type, class P2Type >
#define FFADYNCB_H_ADDED_TEMPLATE_ARGS2 , class P1Type, class P2Type
#define FFADYNCB_H_TYPE_LIST2           < P1Type, P2Type >
#define FFADYNCB_H_ADDED_TYPE_LIST2     , P1Type, P2Type
#define FFADYNCB_H_ARGUMENT_LIST2       P1Type p1, P2Type p2
#define FFADYNCB_H_PARAM_LIST2          p1, p2

// Three types :

#define FFADYNCB_H_TEMPLATE3            template< class P1Type, class P2Type , class P3Type >
#define FFADYNCB_H_ADDED_TEMPLATE_ARGS3 , class P1Type, class P2Type, class P3Type
#define FFADYNCB_H_TYPE_LIST3           < P1Type, P2Type , P3Type >
#define FFADYNCB_H_ADDED_TYPE_LIST3     , P1Type, P2Type, P3Type
#define FFADYNCB_H_ARGUMENT_LIST3       P1Type p1, P2Type p2, P3Type p3
#define FFADYNCB_H_PARAM_LIST3          p1, p2, p3

// Four types :

#define FFADYNCB_H_TEMPLATE4            template< class P1Type, class P2Type , class P3Type, class P4Type >
#define FFADYNCB_H_ADDED_TEMPLATE_ARGS4 , class P1Type, class P2Type, class P3Type, class P4Type
#define FFADYNCB_H_TYPE_LIST4           < P1Type, P2Type , P3Type, P4Type >
#define FFADYNCB_H_ADDED_TYPE_LIST4     , P1Type, P2Type, P3Type, P4Type
#define FFADYNCB_H_ARGUMENT_LIST4       P1Type p1, P2Type p2, P3Type p3, P4Type p4
#define FFADYNCB_H_PARAM_LIST4          p1, p2, p3, p4



//
// Macro that defines the different dynCBN templates
// One macro call for each number of arguments we want to support :
//

#define MAKE_BASE_DYNCB(nr, templateDef, argumentList)	\
							\
templateDef class FFaBaseDynCB##nr			\
{							\
public:							\
  FFaBaseDynCB##nr () {}				\
  virtual ~FFaBaseDynCB##nr () {}			\
							\
  virtual void invoke (argumentList) = 0;		\
  virtual FFaBaseDynCB##nr * copy () const = 0;		\
};


MAKE_BASE_DYNCB(0, ,)
MAKE_BASE_DYNCB(1, FFADYNCB_H_TEMPLATE1, FFADYNCB_H_ARGUMENT_LIST1)
MAKE_BASE_DYNCB(2, FFADYNCB_H_TEMPLATE2, FFADYNCB_H_ARGUMENT_LIST2)
MAKE_BASE_DYNCB(3, FFADYNCB_H_TEMPLATE3, FFADYNCB_H_ARGUMENT_LIST3)
MAKE_BASE_DYNCB(4, FFADYNCB_H_TEMPLATE4, FFADYNCB_H_ARGUMENT_LIST4)


#define MAKE_METHODCB(nr, addedTemplateArgs, argumentList, parameterList, typeList, addedTypeList)	\
													\
template<class T addedTemplateArgs> class FFaMethodCB##nr : public FFaBaseDynCB##nr typeList		\
{													\
  typedef void (T::*methodType) (argumentList);								\
													\
public:													\
  FFaMethodCB##nr (T* anObj, methodType aMethod) : myObject(anObj), myMethodPt(aMethod) {}		\
  virtual ~FFaMethodCB##nr () {}									\
													\
  virtual void invoke (argumentList) { if (myObject) (myObject->*this->myMethodPt) (parameterList); }	\
													\
  virtual FFaBaseDynCB##nr typeList * copy () const							\
  {													\
    return new FFaMethodCB##nr <T addedTypeList> (myObject, myMethodPt);				\
  }													\
													\
private:												\
  T* myObject;												\
  methodType myMethodPt;										\
};

MAKE_METHODCB(0, ,,,,)
MAKE_METHODCB(1, FFADYNCB_H_ADDED_TEMPLATE_ARGS1, FFADYNCB_H_ARGUMENT_LIST1, FFADYNCB_H_PARAM_LIST1, FFADYNCB_H_TYPE_LIST1, FFADYNCB_H_ADDED_TYPE_LIST1)
MAKE_METHODCB(2, FFADYNCB_H_ADDED_TEMPLATE_ARGS2, FFADYNCB_H_ARGUMENT_LIST2, FFADYNCB_H_PARAM_LIST2, FFADYNCB_H_TYPE_LIST2, FFADYNCB_H_ADDED_TYPE_LIST2)
MAKE_METHODCB(3, FFADYNCB_H_ADDED_TEMPLATE_ARGS3, FFADYNCB_H_ARGUMENT_LIST3, FFADYNCB_H_PARAM_LIST3, FFADYNCB_H_TYPE_LIST3, FFADYNCB_H_ADDED_TYPE_LIST3)
MAKE_METHODCB(4, FFADYNCB_H_ADDED_TEMPLATE_ARGS4, FFADYNCB_H_ARGUMENT_LIST4, FFADYNCB_H_PARAM_LIST4, FFADYNCB_H_TYPE_LIST4, FFADYNCB_H_ADDED_TYPE_LIST4)


#define MAKE_STATICCB(nr, templateDef, argumentList, parameterList, typeList)	\
										\
templateDef class FFaStaticCB##nr : public FFaBaseDynCB##nr typeList		\
{										\
  typedef void (*functionType)(argumentList);					\
										\
public:										\
  FFaStaticCB##nr (functionType aFunction) : myFunctionPt(aFunction) {}		\
  virtual ~FFaStaticCB##nr () {}						\
										\
  virtual void invoke (argumentList) { myFunctionPt(parameterList); }		\
										\
  virtual FFaBaseDynCB##nr typeList * copy () const				\
  {										\
    return new FFaStaticCB##nr typeList (myFunctionPt);				\
  }										\
										\
private:									\
  functionType myFunctionPt;							\
};

MAKE_STATICCB(0, ,,,)
MAKE_STATICCB(1, FFADYNCB_H_TEMPLATE1, FFADYNCB_H_ARGUMENT_LIST1, FFADYNCB_H_PARAM_LIST1, FFADYNCB_H_TYPE_LIST1)
MAKE_STATICCB(2, FFADYNCB_H_TEMPLATE2, FFADYNCB_H_ARGUMENT_LIST2, FFADYNCB_H_PARAM_LIST2, FFADYNCB_H_TYPE_LIST2)
MAKE_STATICCB(3, FFADYNCB_H_TEMPLATE3, FFADYNCB_H_ARGUMENT_LIST3, FFADYNCB_H_PARAM_LIST3, FFADYNCB_H_TYPE_LIST3)
MAKE_STATICCB(4, FFADYNCB_H_TEMPLATE4, FFADYNCB_H_ARGUMENT_LIST4, FFADYNCB_H_PARAM_LIST4, FFADYNCB_H_TYPE_LIST4)


#define MAKE_DYNCB(nr, templateDef, argumentList, parameterList, typeList)			\
												\
templateDef class FFaDynCB##nr									\
{												\
public:												\
  FFaDynCB##nr (FFaBaseDynCB##nr typeList *aBaseDynCB = 0) : myBaseDynCBPtr(aBaseDynCB) {}	\
  FFaDynCB##nr (const FFaDynCB##nr typeList &rhs) 						\
  {												\
    myBaseDynCBPtr = rhs.myBaseDynCBPtr ? rhs.myBaseDynCBPtr->copy() : 0;			\
  }												\
  ~FFaDynCB##nr () { if (myBaseDynCBPtr) delete myBaseDynCBPtr; }				\
												\
  void invoke(argumentList) { if (myBaseDynCBPtr) myBaseDynCBPtr->invoke(parameterList); }	\
												\
  FFaDynCB##nr &operator= (const FFaDynCB##nr typeList &rhs)					\
  {												\
    if (this != &rhs)										\
    {												\
      if (myBaseDynCBPtr) delete myBaseDynCBPtr;						\
      myBaseDynCBPtr = rhs.myBaseDynCBPtr ? rhs.myBaseDynCBPtr->copy() : 0;			\
    }												\
    return *this;										\
  }												\
												\
  bool empty() const { return myBaseDynCBPtr == 0; }						\
  void erase() { if (myBaseDynCBPtr) delete myBaseDynCBPtr; myBaseDynCBPtr = 0; }		\
  FFaDynCB##nr* copy() const { return new FFaDynCB##nr (*this); }				\
												\
private:											\
  FFaBaseDynCB##nr typeList *myBaseDynCBPtr;							\
};

MAKE_DYNCB(0, ,,,)
MAKE_DYNCB(1, FFADYNCB_H_TEMPLATE1, FFADYNCB_H_ARGUMENT_LIST1, FFADYNCB_H_PARAM_LIST1, FFADYNCB_H_TYPE_LIST1)
MAKE_DYNCB(2, FFADYNCB_H_TEMPLATE2, FFADYNCB_H_ARGUMENT_LIST2, FFADYNCB_H_PARAM_LIST2, FFADYNCB_H_TYPE_LIST2)
MAKE_DYNCB(3, FFADYNCB_H_TEMPLATE3, FFADYNCB_H_ARGUMENT_LIST3, FFADYNCB_H_PARAM_LIST3, FFADYNCB_H_TYPE_LIST3)
MAKE_DYNCB(4, FFADYNCB_H_TEMPLATE4, FFADYNCB_H_ARGUMENT_LIST4, FFADYNCB_H_PARAM_LIST4, FFADYNCB_H_TYPE_LIST4)

#endif
