/**
 * @file BI_InputError.h
 *   Declarations for exception handling
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef BI_INPUTERROR_H
#define BI_INPUTERROR_H

#include <string>

namespace BEInput
{
//==================================================================================================================================
//! Base class for BlockInput Errors derives from the std::exception class
/*!
 *  This is the base class for all exceptions thrown from the %BlockInput
 *  Package. It derives from std::exception so a general catch can
 *  be used to print out fatal errors from the package.
 */
class BI_InputError : public std::exception
{
public:

    //! Constructor
    /*!
     * @param[in] procedure Procedure where it happened
     * @param[in] msg       Message to be printed out
     */
    explicit BI_InputError(std::string procedure, std::string msg);

    //! Returns the error message string
    /*!
     *   @return           Returns the error message
     */
    std::string errorMessage() const;

    //! virtual destructor
    virtual ~BI_InputError() throw();

    //! Add to the storred message
    /*!
     *  @param   msg      String to add to the storred messaged
     */
    void append(std::string msg);

    //! Function that is inherited from std::exception
    /*!
     *  Guarrantteed not to throw
     *
     *   @return   Returns a null terminated cstring explaining what the error is
     */
    virtual const char* what() const throw();
protected:

    //! String that gets printed out
    std::string m_procedure_msg;
};
//==================================================================================================================================
//! An Unknown Keyline was encountered
class BI_UnknownKeyLine : public BI_InputError
{
public:

    //! Constructor
    /*!
     * @param procedure Procedure where it happened
     * @param lineTok   Value of the LineToken that is not recognized
     */
    BI_UnknownKeyLine(std::string procedure, std::string lineTok);
};
//==================================================================================================================================
//! Unknown BlockEntry name
class BI_UnknownSubBlock : public BI_InputError
{
public:

    //! Constructor
    /*!
     * @param procedure Procedure where it happened
     * @param lineTok   Value of the LineToken that is not recognized
     */
    BI_UnknownSubBlock(std::string procedure, std::string lineTok);
};
//==================================================================================================================================
//! Unknown list Entry was encountered
class BI_UnknownListEntry : public BI_InputError
{
public:

    //! Constructor
    /*!
     * @param procedure Procedure where it happened
     * @param lineTok   Value of the LineToken for the LineEntry
     *                  causing the problem
     */
    BI_UnknownListEntry(std::string procedure, std::string lineTok);
};
//==================================================================================================================================
//! Logic error concerning end of block mismatches
class BI_EndBlockMismatch : public BI_InputError
{
public:

    //! Constructor
    /*!
     * @param procedure Procedure where it happened
     * @param lineTok   Value of the LineToken for the block that is in error
     * @param expected  Expected LineToken
     */
    BI_EndBlockMismatch(std::string procedure, std::string lineTok, std::string expected);
};
//==================================================================================================================================
//! Mismatch concerning the required number of line entries
class BI_MissingRequired : public BI_InputError
{
public:

    //! Constructor
    /*!
     * @param procedure Procedure where it happened
     * @param lineTok   Value of the LineToken causing the error
     * @param numR   Number required
     * @param numP   number found
     */
    BI_MissingRequired(std::string procedure, std::string lineTok, int numR, int numP);
};
//==================================================================================================================================
//! Mismatch concerning the number of vector entries
class BI_MissingRequiredVec : public BI_InputError
{
public:

    //! Constructor
    /*!
     * @param procedure Procedure where it happened
     * @param lineTok   Value of the LineToken causing the error
     * @param numVneeded Number of items needed
     * @param numVfound  Number of items found
     */
    BI_MissingRequiredVec(std::string procedure, std::string lineTok, int numVneeded, int numVfound);
};
//==================================================================================================================================
//! Set the print level
/*!
 * @param ilvl print lvl 1 -> printing
 *                       0 -> no printing
 */
void BI_SetPrintLevel(int ilvl);
//==================================================================================================================================
}
#endif
