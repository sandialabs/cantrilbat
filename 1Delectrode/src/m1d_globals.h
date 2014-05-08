/**
 * @file m1d_globals.h
 * This file contains definitions for globally defined functions and data
 * that are not part of any class.
 */

/*
 *   $Id: m1d_globals.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef _M1D_GLOBALS_H
#define _M1D_GLOBALS_H

class Epetra_Comm;

#include <string>
#include <iostream>
#include <fstream>

/*
 * This file contains globally defined functions and data
 * that are not
 * part of any class, but are still part of the m1d namespace.
 *
 *
 *
 */
namespace m1d
{

//!  Global Epetra comm pointer. There is never more than one of these in this application.
extern Epetra_Comm *Comm_ptr;

class Appl;
class ProblemStatement;

extern int s_printLvl_TimingInformation;

extern int s_printLvl_DebugTables;

extern void readEnvironmentalVariables();


//! Print input options and exit
/*!
 *  If this is true, then we print all of the entire input options, as much as possible,
 *  and then we exit before starting the calculation
 */
extern bool PrintInputFormat;

//! Problem statement pointer. This must be defined in the calling program
/*!
 *   This is allocated in m1d_ProblemStatement.cpp
 */
extern ProblemStatement* PSinput_ptr;

//! Return a pointer to the single application object
Appl*
app();

//! Set an error condition in the application class without
//! throwing an exception.
/*!
 * This routine adds an error message to the end of the stack
 * of errors that m1d accumulates in the Appl class.
 * \ingroup errorhandling
 *
 * @param r  Location of the error
 * @param msg  String containing a description of the error
 */
void
setError(std::string r, std::string msg);

//! Prints all of the error messages to stream f.
/*!
 * Write out to ostream, f, all of the saved error messages.
 * Cantera saves a stack of exceptions that it
 * has caught in the Application class. This routine writes
 * out all of the error messages to ostream f, and then
 * clears them from internal storage.
 * @ingroup errorhandling
 */
void
showErrors(std::ostream& f);

//! Discard the last error message
/*!
 *  saves a stack of exceptions that it
 * has caught in the Application class. This routine eliminates
 * the last exception to be added to that stack.
 *
 * @ingroup errorhandling
 */
void
popError();

std::string
intToString(const int p);
}

#endif
