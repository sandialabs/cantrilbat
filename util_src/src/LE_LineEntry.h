/**
 * @file LE_LineEntry.h
 *   Declarations for a LineEntry object
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */
/*
 *
 */
#ifndef LE_LINEENTRY_H
#define LE_LINEENTRY_H

#include "tok_input_util.h"
#include "BaseEntry.h"
#include "BI_Dependency.h"

namespace BEInput
{
//==================================================================================================================================
//! A LineEntry object is a single-line entry that can be matched
//! and then acted upon
/*!
 * This class is the base class for handling these single
 * line entries. It declares interfaces for specifying the
 * entries and for processing the entries once a match is
 * made in the input deck.
 *
 * The general format of a %LineEntry in the input deck is:
 *
 *    \f[
 *      KeyName = Arguments
 *    \f]
 *
 * Both the <I>KeyName</I> and the <I>Arguments</I> variables
 * may be consist of multiple
 * tokens. Each is tokenized and compared in a case insensitive format
 * at the start of their processing and before any matching
 * is carried out.
 *
 * We have made this base class a marginally pure virtual class.
 * You can't make an instance of this class; you must
 * make an instance of a derivative of this class.
 *
 * @ingroup blockentryModule
 */
class LineEntry : public BaseEntry
{

public:

    //! Constructor
    /*!
     *  This is the constructor for the %LineEntry object. We enter the
     *  number of required lines of this object and the name of the line
     *  entry here.
     *
     * @param keyName  KeyName of the object
     * @param numRL    Number of required lines.
     */
    LineEntry(const char* keyName, int numRL = 0);

    //! Copy Constructor
    /*!
     * @param b parameter to be copied
     */
    LineEntry(const LineEntry& b);

    //!  Assignment operator
    /*!
     * @param[in] b object to be copied
     *
     *   @return    Returns a reference to the current object
     */
    LineEntry& operator=(const LineEntry& b);

    //! Destructor
    virtual ~LineEntry();

    //! Duplicator for a BaseEntry object
    /*!
     *         This just calls the duplMyselfAsLineEntry command.
     *
     *   @return            Returns a pointer to the duplicate object as a BaseEntry object
     */
    virtual BaseEntry* duplMyselfAsBaseEntry() const;

    //! Duplicator function. This is a pure virtual function here
    /*!
     *   This function duplicates the entry and returns a pointer to a LineEntry. This is a pure virtual function
     *   making direct instances of this class not possible.
     *
     *   @return            Returns a pointer to the duplicate object as a LineEntry object
     */
    virtual LineEntry* duplMyselfAsLineEntry() const = 0;

    //! Process this LineEntry
    /*!
     * This function is called when it has been determined that
     * the current KeyName from the input file matches this
     * object's KeyName.
     *
     * This function then processes the entry.
     * Processing involves checking for the satisfaction
     * of runtime dependencies.
     * Then, the current value is assigned from the token arguments of the LineEntry
     * and the external pointer to the data is assigned the current value.
     *
     *  What's in this fucntion is just a stub routine. The main work is done by
     *  the member function <TT>process_LineEntry()</TT> functions
     *  of derived classes of the LineEntry class.
     *  In this stub routine, we increment the NumProcessedLines variable and
     *   we do optional printing, i.e., everything that is common to
     *   al LineEntry objects.
     *   This stub routine may be called be every LineEntry child class before it
     *   does its own specific processing.
     *
     * @param lineArgTok  Pointer to the token containing the
     *                    arguments to the lineEntry. This is
     *                    everything after the "=" sign.
     */
    virtual void process_LineEntry(const TK_TOKEN* lineArgTok);

    //!  This routine will print out a processed line
    /*!
     *  The default behavior is to print the original line with a "====>"
     *  prefix to indicate that action has been taken on it.
     *
     * @param lineArgTok  Pointer to the token containing the
     *                    arguments to the lineEntry. This is
     *                    everything after the "=" sign.
     */
    virtual void print_ProcessedLine(const TK_TOKEN* lineArgTok) const;

    //! Print out API information about this keyline
    /*!
     * This routine will print out to stdout information about
     * the keyline. This command is used to document the interface.
     *
     * @param ilvl Level of the indentation to use in printing
     */
    virtual void print_usage(int ilvl) const = 0;

protected:

    //! Static function that prints an indent
    /*!
     *
     * @param ilvl Level of the indentation to use in printing
     */
    static void print_indent(int ilvl);

    //! Adjust base address to store the value
    /*!
     * @param addrAdjustment                 Offset in raw bytes to adjust the internal value of the Address that will receive
     *                                       the int.
     */
    virtual void adjustAddress(LONG_PTR addrAdjustment) = 0;

    //! Check for requirements being met at the end of input
    /*!
     * @param[in]   throwSpecificError      If true then you should throw a specific error condition, if you have one. If not,
     *                                      an generic error condition will be thrown on return.
     *
     *  @return                             Returns a boolean indicating whether the specific requirements have been met.
     */
    virtual bool checkRequirements(bool throwSpecificError = false);

    //! class BlockEntry is a friend
    /*!
     * Declare BlockEntry a friend. That way we can used protected functions
     * adjustAddress() and checkRequirements() within BlockEntry functions.
     */
    friend class BlockEntry;
};
}
#endif
