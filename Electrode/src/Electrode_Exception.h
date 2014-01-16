/**
 * @file m1d_exceptions.h
 *  Declarations for error handling within m1d
 */

/*
 *  $Id: m1d_exception.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef ELECTRODE_EXCEPTION_H
#define ELECTRODE_EXCEPTION_H

#include <string>
#include <exception>
#include <vector>

#include "cantera/base/ctexceptions.h" 

namespace Cantera 
{

//! base error class for m1d package inherits from the exception stl
/*!
 *
 */
class Electrode_Error : public Cantera::CanteraError {
public:

  //! Normal Constructor for the m1d_Error base class
  /*!
   * This class doesn't have any storage associated with it. In its
   * constructor, a call to the Application class is made to store
   * the strings associated with the generated error condition.
   *
   * @param procedure String name for the function within which the error was
   *             generated.
   * @param msg  Descriptive string describing the type of error message.
   */
  Electrode_Error(const std::string &procedure, const std::string &msg);

  //! Destructor for base class does nothing
  virtual ~Electrode_Error() throw ();

protected:

  //! Empty base constructor is made protected so that it may be used only by
  //! inherited classes.
  /*!
   *  We want to discourage throwing an error containing no information.
   */
  Electrode_Error();
};

//! Provides a line number
#define XSTR_TRACE_LINE(s) STR_TRACE_LINE(s)

//! Provides a line number
#define STR_TRACE_LINE(s) #s

//! Provides a std::string variable containing the file and line number
/*!
 *   This is a std:string containing the file name and the line number
 */
#define STR_TRACE   (std::string(__FILE__) +  ":" + XSTR_TRACE_LINE(__LINE__))

#ifdef AssertTrace
#undef AssertTrace
#endif

#ifdef AssertThrow
#undef AssertThrow
#endif

#ifdef AssertThrowMsg
#undef AssertThrowMsg
#endif

#ifdef NDEBUG
#  define AssertTrace(expr)                        ((void) (0))
#  define AssertThrow(expr, procedure)             ((void) (0))
#  define AssertThrowMsg(expr,procedure, message)  ((void) (0))
#else

//! Assertion must be true or an error is thrown
/*!
 * Assertion must be true or else a ElectrodeError is thrown. A diagnostic string containing the
 * file and line number,  indicating where the error
 * occured is added to the thrown object.
 *
 * @param expr  Boolean expression that must be true
 *
 * @ingroup errorhandling
 */
#  define AssertTrace(expr)  ((expr) ? (void) 0 : \
			      throw Cantera::Electrode_Error(STR_TRACE, std::string("failed assert: ") + #expr))

//!  Assertion must be true or an error is thrown
/*!
 * Assertion must be true or else a CanteraError is thrown. A diagnostic string indicating where the error
 * occurred is added to the thrown object.
 *
 * @param expr       Boolean expression that must be true
 * @param procedure  Character string or std:string expression indicating the procedure where the assertion failed
 * @ingroup errorhandling
 */
#  define AssertThrow(expr, procedure)   ((expr) ? (void) 0 :\
					  throw Cantera::Electrode_Error(procedure, std::string("failed assert: ") + #expr))

//!  Assertion must be true or an error is thrown
/*!
 * Assertion must be true or else a CanteraError is thrown. A
 * diagnostic string indicating where the error occurred is added
 * to the thrown object.
 *
 * @param expr       Boolean expression that must be true
 * @param procedure  Character string or std:string expression indicating
 *                   the procedure where the assertion failed
 * @param message  Character string or std:string expression contaiing
 *    a descriptive message is added to the thrown error condition.
 *
 * @ingroup errorhandling
 */

# define AssertThrowMsg(expr, procedure, message) \
             ((expr) ? (void) 0 : throw Cantera::Electrode_Error(procedure + std::string(": at failed assert: \"") +\
                                                       std::string(#expr) + std::string("\""), message) )

#endif
        // ==================================================================================================
      } // End of m1d namespace
      // ==================================================================================================
#endif

