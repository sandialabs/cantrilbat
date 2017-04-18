/**
 * @file m1d_exception.h  Declarations for error handling within m1d
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef M1D_EXCEPTION_H
#define M1D_EXCEPTION_H

#include <string>
#include <exception>

#include <cstdarg>
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! Base error class for m1d package. Inherits from the std exception class
/*!
 *  This causes an immediate error exit from the program.
 *
 *  @todo:   This class needs to handle mp program exits better.
 */
class m1d_Error : public std::exception
{
public:

    //! Normal Constructor for the m1d_Error base class
    /*!
     *  This class doesn't have any storage associated with it. In its constructor, a call to the Application class is made to store
     *  the strings associated with the generated error condition.
     *
     * @param[in]              procedure           String name for the function within which the error was generated.
     * @param[in]              msg                 Descriptive string describing the type of error message.
     */
    m1d_Error(const std::string& procedure, const std::string& msg);

    //! printf-like constructor for the m1d_Error base class
    /*!
     * The message is formatted according to the standard c printing routines.
     *
     * @param[in]              procedure           String name for the function within which the error was generated.
     * @param[in]              fmt                 printf-like format string
     *
     *  Add parameters for fmt string according to the printf, fprintf man pages
     */
    m1d_Error(const std::string& procedure, const char* fmt, ...);

    //! Destructor for base class does nothing
    virtual ~m1d_Error() throw ()
    {
    }

protected:

    //! Message associated with the error
    std::string msg_;

    //! Empty base constructor is made protected so that it may be used only by inherited classes.
    /*!
     *  We want to discourage throwing an error containing no information.
     */
    m1d_Error()
    {
    }
};

//==================================================================================================================================
//! Assert two number are equal up to a number of digits
/*!
 *  @param[in] a1          First double
 *  @param[in] a2          Second double
 *  @param[in] atol        Absolute tolerance. Number below this value are not considered to be different.
 *
 *  @param[in] digits      Number of digits of accuracy to be considered. Defaults to 13
 *
 *  @return                Returns true for equality, and false for inequality
 */
bool doubleEqual(double a1, double a2, double atol = 1.0E-13, int digits = 13);

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

//==================================================================================================================================
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
 * Assertion must be true or else a CanteraError is thrown. A diagnostic string containing the
 * file and line number,  indicating where the error
 * occured is added to the thrown object.
 *
 * @param expr  Boolean expression that must be true
 *
 * @ingroup errorhandling
 */
#  define AssertTrace(expr)  ((expr) ? (void) 0 : \
			      throw m1d::m1d_Error(STR_TRACE, std::string("failed assert: ") + #expr))

//==================================================================================================================================
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
					  throw m1d::m1d_Error(procedure, std::string("failed assert: ") + #expr))

//==================================================================================================================================
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
             ((expr) ? (void) 0 : throw m1d::m1d_Error(procedure + std::string(": at failed assert: \"") +\
                                                       std::string(#expr) + std::string("\""), message) )

#endif
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

