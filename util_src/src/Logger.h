/**
 * @file Logger.h
 * Header for Base class for 'loggers' that write text messages to log files
 */
#ifndef MDPUTIL_LOGGER_H
#define MDPUTIL_LOGGER_H

#include <iostream>
#include <cstdlib>

namespace mdpUtil
{

//! Base class for 'loggers' that write text messages to log files.
/*!
 *      This class is used to direct log messages to application- or
 *      environment-specific output. The default is to simply print
 *      the messages to the standard output stream or standard error
 *      stream, but classes may be derived from Logger that implement
 *      other output options. This is important when Cantera is used
 *      in applications that do not display the standard output, such
 *      as MATLAB.  The Cantera MATLAB interface derives a class from
 *      Logger that implements these methods with MATLAB-specific
 *      procedures, insuring that the messages will be passed through
 *      to the user. It would also be possible to derive a class that
 *      displayed the messages in a pop-up window, or redirected them
 *      to a file, etc.
 *
 *      To install a logger, call function setLogger (global.h / misc.cpp).
 *
 *      @ingroup textlogs
 */
class Logger
{
public:

    //! Constructor - empty
    Logger();

    //! Destructor - empty
    virtual ~Logger();

    //! Write a log message.
    /*!
     *  The default behavior is to write to
     *  the standard output. Note that no end-of-line character is
     *  appended to the message, and so if one is desired it must
     *  be included in the string.
     *
     * @param msg      String message to be written to cout
     */
    virtual void write(const std::string& msg);

    //! Write an end of line character and flush output.
    /*!
     *  Some systems treat endl and \n differently. The endl
     *  statement causes a flushing of stdout to the screen.
     */
    virtual void writeendl();

    //! Set the location of output
    /*!
     *   @param[in] outputLocation  Location for the output
     */
    virtual void setOutputStream(std::ostream& outputLocation);

    //! Write an error message and quit.
    /*!
     *  The default behavior is
     *  to write to the standard error stream, and then call
     *  exit(). Note that no end-of-line character is appended to
     *  the message, and so if one is desired it must be included
     *  in the string. Note that this default behavior will
     *  terminate the application Cantera is invoked from (MATLAB,
     *  Excel, etc.) If this is not desired, then derive a class
     *  and reimplement this method.
     *
     * @param msg    Error message to be written to cerr.
     */
    virtual void error(const std::string& msg);

    //! Name of the output stream
    std::ostream* out_;
    
};

}
#endif