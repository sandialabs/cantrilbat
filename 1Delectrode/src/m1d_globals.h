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
#include <vector>
#include <ostream>

/*
 * This file contains globally defined functions and data
 * that are not
 * part of any class, but are still part of the m1d namespace.
 *
 *
 *
 */
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{

//!  Global Epetra comm pointer. There is never more than one of these in this application.
extern Epetra_Comm* Comm_ptr;

//! Application Class. There is nover more than one of these in the application
class Appl;

class ProblemStatement;

//==================================================================================================================================
//! Level of timing information
extern int s_printLvl_TimingInformation;

//! External int specifing the Print level for the DebugTables
extern int s_printLvl_DebugTables;
//==================================================================================================================================
//! Function to read in environmental variables that are set in the user environment
/*!
 *  Environmental variables recognized:
 *        M1D_PRINT_TIMINGINFORMATION
 *        M1D_PRINT_DEBUGTABLES
 */
extern void readEnvironmentalVariables();

//==================================================================================================================================
//! Global Print input options and exit
/*!
 *  If this is true, then we print all of the entire input options, as much as possible,
 *  and then we exit before starting the calculation
 */
extern bool PrintInputFormat;

//==================================================================================================================================
//! Global Problem statement pointer. This must be defined in the calling program
/*!
 *   This is allocated in m1d_ProblemStatement.cpp
 */
extern ProblemStatement* PSinput_ptr;

//==================================================================================================================================
//! Return a pointer to the single application object
/*!
 *  @return                                      Returns a pointer to the single application object
 */
Appl* app();

//==================================================================================================================================
//! Set an error condition in the application class without throwing an exception.
/*!
 *  This routine adds an error message to the end of the stack
 *  of errors that m1d accumulates in the Appl class.
 *  \ingroup errorhandling
 *
 *  @param[in]               r                   Location of the error
 *  @param[in]               msg                 String containing a description of the error
 */
void setError(std::string r, std::string msg);

//==================================================================================================================================
//! Prints all of the error messages to stream f.
/*!
 *  Write out to ostream, f, all of the saved error messages.  Zuzax saves a stack of exceptions that it
 *  has caught in the Application class. This routine writes out all of the error messages to ostream f, and then
 *  clears them from internal storage.
 *
 *  @param[in]               f                   Reference to an ostream to print the errors
 *  @ingroup errorhandling
 */
void showErrors(std::ostream& f);

//==================================================================================================================================
//! Discard the last error message
/*!
 *  Saves a stack of exceptions that it has caught in the Application class. This routine eliminates
 * the last exception to be added to that stack.
 *
 * @ingroup errorhandling
 */
void popError();

//==================================================================================================================================
//! Global function to change an int into a string
/*!
 *   @param[in]              p                   Integer to be converted
 *
 *  @return                                      Returns a string representation
 */
std::string intToString(const int p);

//==================================================================================================================================
//! Check error between two doubles
/*!
 *  @param[in]               d1                  The first double
 *  @param[in]               d2                  The second double
 *  @param[in]               rtol                The relative tolerance. Defaults to 1.0E-4
 *  @param[in]               atol                the absolute tolerance. Defaults to 1.0E-13
 *  @return                                      Returns true if the doubles agree and false otherwise
 */
bool checkDblAgree(double d1, double d2, double rtol = 1.0E-4, double atol = 1.0E-13);

//==================================================================================================================================
//! Write out a tecplot vector to an ascii file
/*!
 *  Write out a vector of doubles to the output file
 *
 *  @param[in]               ofp        FILE pointer, Must already be open
 *  @param[in]               vals       vector of doubles
 *  @param[in]               numD       Number of decimal places to be written out. Default is 13
 *  @param[in]               numPerLine Number of values per line. The default is 10.
 */
void fwriteTecplotVector(FILE* ofp, const std::vector<double>& vals, int numD = 13, int numPerLine = 10);

//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
