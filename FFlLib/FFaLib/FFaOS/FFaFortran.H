/* SPDX-FileCopyrightText: 2023 SAP SE
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * This file is part of FEDEM - https://openfedem.org
 */
/*!
  \file FFaFortran.H
  \brief Macros to aid cross-language linkage between C/C++ and Fortran.
*/

#ifndef FFA_FORTRAN_H
#define FFA_FORTRAN_H

/*!
  \brief Integer type corresponding to the Fortran integer kind kindmodule::i8.
  \sa kindModule.f90
*/
#if defined(win32) || defined(win64)
#if defined(FT_HAS_INT4_ONLY)
#define f90_int8 __int32
#else
#define f90_int8 __int64
#endif
#elif defined(irix)
#define f90_int8 long long
#elif defined(hpux)
#define f90_int8 off64_t
#else
#define f90_int8 long
#endif

/*!
  \brief Integer type corresponding to the default Fortran integer kind.
*/
#if defined(FT_USE_INT8)
#define f90_int f90_int8
#else
#define f90_int int
#endif

/*!
  \brief Integer type used in type casting to/from C/C++ pointers and Fortran.
  \sa pointerkindmodule in FFrExtractorInterface.f90
*/
#if defined(irix64) || defined(hpux64) || defined(aix64) || defined(win64) || defined(linux64)
#define ptr_int f90_int8
#else
#define ptr_int int
#endif

/*!
  \brief Macro for portable handling of name conventions in mixed C/C++/Fortran.
  \details Since there is no way(?) to translate a lowercase string
  into an uppercase string through preprocessor macros, we must specify both
  the lower- and uppercase version of the name.
*/
#if defined(win32) || defined(win64)
#define F90_NAME(name,NAME) NAME
#elif defined(hpux) || defined(aix) || defined(aix64)
#define F90_NAME(name,NAME) name
#else
#define F90_NAME(name,NAME) name##_
#endif

/*! \brief Defines the C++ linkage type on Windows. */
#if _MSC_VER > 1300
#define _CPFX _cdecl
#elif defined(win32)
#define _CPFX _stdcall
#define _NCHAR_AFTER_CHARARG
#else
#define _CPFX
#endif

/*! \brief Macro for portable handling of linkage in mixed C/C++/F90 calls. */
#ifdef __cplusplus
#define F90_FUNCTION(type,name,NAME) extern "C" type _CPFX F90_NAME(name,NAME)
#else
#define F90_FUNCTION(type,name,NAME) type _CPFX F90_NAME(name,NAME)
#endif


/*! \brief Macro for declaration for Fortran subroutines in C/C++. */
#define       SUBROUTINE(name,NAME) F90_FUNCTION(void,name,NAME)
/*! \brief Macro for declaration for Fortran pointer functions in C/C++. */
#define     PTR_FUNCTION(name,NAME) F90_FUNCTION(void*,name,NAME)
/*! \brief Macro for declaration for Fortran integer functions in C/C++. */
#define INTEGER_FUNCTION(name,NAME) F90_FUNCTION(f90_int,name,NAME)
/*! \brief Macro for declaration for Fortran real functions in C/C++. */
#define    REAL_FUNCTION(name,NAME) F90_FUNCTION(float,name,NAME)
/*! \brief Macro for declaration for Fortran double functions in C/C++. */
#define  DOUBLE_FUNCTION(name,NAME) F90_FUNCTION(double,name,NAME)

#endif
