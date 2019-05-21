/**
 * @file m1d_defs.h
 *  Declarations that are included in all files within the m1d namespace
 */

/*
 * Copywrite 2013 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef M1D_DEFS_H
#define M1D_DEFS_H

/*
 *   Include this file at the top so that Epetra's configuration file is always called before the
 *   contents of our config.h file
 */
#include "Epetra_ConfigDefs.h"

/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for each package and need to
 * be undef'd here to avoid warnings when this file is included from another package.
 * KL 11/25/02
 */
#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif
#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif
#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif
#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif
#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif
#ifdef HAVE_MPI
#undef HAVE_MPI
#endif

#include "configCantrilbat.h"
#include "m1d_exception.h"

#include "zuzax/base/ct_defs.h"
#include "zuzax/base/stringUtils.h"


//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{

//! Macro to describe how to access the data pointer underlying stl vectors directly
/*!
 *  Note we can now use the data() member function, so this is now decremented
 */
#ifndef DATA_PTR
#define DATA_PTR(x)  &((x)[0])
#endif

//! Macro to check pointers before deleting their contents.
/*!
 *   This macro will set the pointer to zero, indicating that the contents are now inaccessible
 */
#define safeDelete(ptr)  if (ptr) { delete (ptr) ; ptr = nullptr; }

//! Utility for setting doubles
const double M1D_DOUBLE_NOTSET = -1.234567E300;

//! Index returned by functions to indicate "no position"
const size_t npos = static_cast<size_t>(-1);

}
//----------------------------------------------------------------------------------------------------------------------------------

#endif
