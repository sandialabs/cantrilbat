/**
 * @file Electrode_Exception.h
 *  Declarations for error handling within the Electrode Object
 *  (see \ref errorHandling and class Electrode_Error
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef ELECTRODE_EXCEPTION_H
#define ELECTRODE_EXCEPTION_H

#include "cantera/base/ctexceptions.h" 
#include "Electrode_defs.h"

#include <string>
#include <vector>
//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif 
{
class Electrode;
}
//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
//! base error class for m1d package inherits from the exception stl
/*!
 *  This error class is built on top of the Zuzax error procedures.
 */
class Electrode_Error : public ZZCantera::ZuzaxError {
public:

  //! Normal Constructor for the m1d_Error base class
  /*!
   * This class doesn't have any storage associated with it. In its constructor, 
   * a call to the Zuzax::Application class is made to store the strings associated with the generated error condition.
   *
   * @param[in]              procedure           String name for the function within which the error was  generated.
   * @param[in]              msg                 Descriptive string describing the type of error message.
   */
  Electrode_Error(const std::string &procedure, const std::string &msg);

  //! printf-like constructor for the Electrode_Error base class
  /*!
   * The message is formatted according to the standard c printing routines.
   *
   * @param[in]              procedure           String name for the function within which the error was generated.
   * @param[in]              fmt                 printf-like format string
   * 
   *  Add parameters for fmt string according to the printf, fprintf man pages
   */
  Electrode_Error(const std::string& procedure, const char* fmt, ...);


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

//==================================================================================================================================
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
//==================================================================================================================================
//! Assertion must be true or an error is thrown
/*!
 * Assertion must be true or else a Electrode_Error is thrown. A diagnostic string containing the
 * file and line number, indicating where the error occurred is added to the thrown object.
 *
 * @param[in]                expr                Boolean expression that must be true
 *
 * @ingroup errorhandling
 */
#  define AssertTrace(expr)  ((expr) ? (void) 0 : \
			      throw ZZCantera::Electrode_Error(STR_TRACE, std::string("failed assert: ") + #expr))

//==================================================================================================================================
//!  Assertion must be true or an error is thrown
/*!
 * Assertion must be true or else a CanteraError is thrown. A diagnostic string indicating where the error
 * occurred is added to the thrown object.
 *
 * @param[in]                expr                Boolean expression that must be true
 * @param[in]                procedure           Character string or std:string expression indicating the procedure where
 *                                               the assertion failed
 * @ingroup errorhandling
 */
#  define AssertThrow(expr, procedure)   ((expr) ? (void) 0 :\
					  throw ZZCantera::Electrode_Error(procedure, std::string("failed assert: ") + #expr))

//==================================================================================================================================
//!  Assertion must be true or an error is thrown
/*!
 *  Assertion must be true or else a CanteraError is thrown. A
 *  diagnostic string indicating where the error occurred is added to the thrown object.
 *
 *  @param[in]               expr                Boolean expression that must be true
 *  @param[in]               procedure           Character string or std:string expression indicating
 *                                               the procedure where the assertion failed
 *  @param message  Character string or std:string expression contaiing
 *    a descriptive message is added to the thrown error condition.
 *
 *  @ingroup errorhandling
 */

# define AssertThrowMsg(expr, procedure, message) \
             ((expr) ? (void) 0 : throw ZZCantera::Electrode_Error(procedure + std::string(": at failed assert: \"") +\
                                                       std::string(#expr) + std::string("\""), message) )

//==================================================================================================================================
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
void Electrode_Warning(const ZZCantera::Electrode& e,  const std::string &procedure, const std::string &msg);

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
bool doubleEqual(double a1, double a2, double atol = 1.0E-300, int digits = 6);

//! Assert two numbers are equal up to a number of digits 
/*!
 *   @param[in]      a1                  First float
 *   @param[in]      a2                  Second float
 *   @param[in]      digits              Number of digits - defaults to 4
 *
 *   @return                             Returns true if floats are equal up to the specification
 */
bool doubleEqualNoAtol(double a1, double a2, int digits = 6);

//! Assert two vectors of doubles are equal up to a number of digits and to an absolute tolerance
/*!
 *   @param[in]      a1                  first float vector
 *   @param[in]      a2                  second float vector
 *   @param[in]      atol                Absolute tolerance - defaults to 1.0E-200, most be positive
 *   @param[in]      digits              Number of digits - defaults to 4
 *
 *   @return                             Returns true if floats are equal up to the specification
 */
bool doubleVectorEqual(const std::vector<double>& a1, const std::vector<double>& a2, double atol = 1.0E-300, int digits = 6);

//! Assert two vectors of doubles are equal up to a number of digits
/*!
 *   This routine asserts that the vector of doubles are all equal up to a number of digits.
 *   @param[in]      a1                  first float vector
 *   @param[in]      a2                  second float vector
 *   @param[in]      digits              Number of digits - defaults to 4
 *
 *   @return                             Returns true if floats are equal up to the specification
 */
bool doubleVectorEqualNoAtol(const std::vector<double>& a1, const std::vector<double>& a2, int digits = 6);

//! Create the L0 relative norm of the difference between two vectors
/*!
 * This creates a relative norm that is scaled by the current value of rtol, and that can always be compared to 1
 * for satisfaction of the error criteria. Assumes all vectors have equal lengths.
 *
 *     L0_i    =     || Actual_i - Pred_i || / ( rtol * MAX(||Actual_i||, || Pred_i ||)
 *
 *   if ||Abstol_i|| > rtol * MAX(||Actual_i||, || Pred_i || 
 * 
 *    L0_i    =     || Actual_i - Pred_i || / ( MAX(||Abstol_i|| )
 *
 *  @param[in]               v1                  First vector to compare against. Lengt
 *  @param[in]               v2                  Second vector to compare against. Length >= num
 *  @param[in]               atolVec             Absolute tolerance vector. Length >= num
 *  @param[in]               rtol                Relative tolerance
 * 
 *   @return                                     returns the maximum value of Lo_i for all i.
 */
double l0norm(const std::vector<double>& v1, const std::vector<double>& v2, const std::vector<double>& atolVec,
              const double rtol);

//! Return a relative difference between two doubles given an absolute tolerance
/*!
 *  Returns the dimensionless value representative of the number of digits of agreement between two number above an 
 *  absolute tolerance level.
 *
 *   rel_diff = fabs( a - b ) / MAX( || a || , || b || , atol , 1.0E-300 )
 *
 *  @param[in]               a                   First double
 *  @param[in]               b                   second double
 *  @param[in]               atol                absolute tolerance
 *
 *  @return                                      Returns the relative normalized difference
 */
double relv(double a, double b, double atol);

//==================================================================================================================================
} // End of esmodel namespace
//----------------------------------------------------------------------------------------------------------------------------------
#endif

