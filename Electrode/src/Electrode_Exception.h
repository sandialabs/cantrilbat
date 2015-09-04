/**
 * @file Electrode_Exception.h
 *  Declarations for error handling within the Electrode Object
 *  (see \ref errorHandling and class Electrode_Error
 */


#ifndef ELECTRODE_EXCEPTION_H
#define ELECTRODE_EXCEPTION_H

#include "cantera/base/ctexceptions.h" 

#include <string>
#include <vector>
//----------------------------------------------------------------------------------------------------------------------------------
namespace Cantera 
{
class Electrode;
}

namespace Cantera
{
//==================================================================================================================================
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
   * @param[in] procedure String name for the function within which the error was
   *             generated.
   * @param[in] msg  Descriptive string describing the type of error message.
   */
  Electrode_Error(const std::string &procedure, const std::string &msg);

  //! Destructor for base class does nothing
  virtual ~Electrode_Error() throw ();

protected:

  //! Empty base constructor is made protected so that it may be used only by inherited classes.
  /*!
   *  We want to discourage throwing an error containing no information.
   */
  Electrode_Error();

};

//==================================================================================================================================

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
 * Assertion must be true or else a Electrode_Error is thrown. A diagnostic string containing the
 * file and line number, indicating where the error occurred is added to the thrown object.
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
}
//----------------------------------------------------------------------------------------------------------------------------------
namespace esmodel
{
//==================================================================================================================================

//! Write a warning to the logfile when within an Electrode object
/*!
 *   Writes a warning to the logfile but the code keeps going.
 *
 *   @param[in]             e                     Reference to the electrode object where the warning originated
 *   @param[in]             procedure             String identifying the procedure where the warning occurred.
 *   @param[in]             msg                   String with the message
 */
void Electrode_Warning(const Cantera::Electrode& e,  const std::string &procedure, const std::string &msg);

//! Write a warning to the logfile when not within an Electrode object
/*!
 *   Writes a warning to the logfile but the code keeps going.
 *
 *   @param[in]             procedure             String identifying the procedure where the warning occurred.
 *   @param[in]             msg                   String with the message
 */
void ESModel_Warning(const std::string &procedure, const std::string &msg);

//==================================================================================================================================
//! Assert two numbers are equal up to a number of digits and to an absolute tolerance
/*!
 *   @param[in]      a1                  First float
 *   @param[in]      a2                  Second float
 *   @param[in]      atol                Absolute tolerance - defaults to 1.0E-200, must be positive
 *   @param[in]      digits              Number of digits - defaults to 4
 *
 *   @return                             Returns true if floats are equal up to the specification
 */
extern bool doubleEqual(double a1, double a2, double atol = 1.0E-300, int digits = 6);

//! Assert two numbers are equal up to a number of digits 
/*!
 *   @param[in]      a1                  First float
 *   @param[in]      a2                  Second float
 *   @param[in]      digits              Number of digits - defaults to 4
 *
 *   @return                             Returns true if floats are equal up to the specification
 */
extern bool doubleEqualNoAtol(double a1, double a2, int digits = 6);

//! Assert two vectors of doubles are equal up to a number of digits and to an absolute tolerance
/*!
 *   @param[in]      a1                  first float vector
 *   @param[in]      a2                  second float vector
 *   @param[in]      atol                Absolute tolerance - defaults to 1.0E-200, most be positive
 *   @param[in]      digits              Number of digits - defaults to 4
 *
 *   @return                             Returns true if floats are equal up to the specification
 */
extern bool doubleVectorEqual(const std::vector<double>& a1, const std::vector<double>& a2, double atol = 1.0E-300, int digits = 6);

//! Assert two vectors of doubles are equal up to a number of digits
/*!
 *   This routine asserts that the vector of doubles are all equal up to a number of digits.
 *   @param[in]      a1                  first float vector
 *   @param[in]      a2                  second float vector
 *   @param[in]      digits              Number of digits - defaults to 4
 *
 *   @return                             Returns true if floats are equal up to the specification
 */
extern bool doubleVectorEqualNoAtol(const std::vector<double>& a1, const std::vector<double>& a2, int digits = 6);

//==================================================================================================================================
} // End of esmodel namespace
//----------------------------------------------------------------------------------------------------------------------------------
#endif

