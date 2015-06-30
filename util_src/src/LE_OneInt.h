/**
 * @file LE_OneInt.h
 *  Header for the LineEntry of a single integer
 *  (see \ref blockentryModule and class
 *  \link BEInput::LE_OneInt LE_OneInt\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef LEONEINT_H
#define LEONEINT_H

#include "LE_LineEntry.h"

/* TEMPLATE FOR THE OTHER CARDS */

namespace BEInput
{

//!  Keyline entry of a single int value.
/*!
 *   This sets up the <b>LineEntry</b> special case for the entry of a single
 *   int input. It declares interfaces for specifying the
 *   entry and for processing the entry once a match is
 *   made in the input deck. The general form of the input is
 *
 *   \f[
 *     KeyName = \mathrm{int}
 *   \f]
 *
 *  The <I>KeyName</I> variable  may be consist of multiple
 *  tokens. It is tokenized and compared in a case insensitive format
 *  at the start of the processing and before any matching
 *  is carried out.
 *  \f$ \mathrm{int} \f$ is a single integer. The code will check that
 *  the rhs consists of one and only one well-formed integer.
 *  Additionally, the following tokens are understood and accepted in
 *  place of \f$ \mathrm{int} \f$ : <tt>INT_MAX</tt>,  <tt>INT_MIN</tt>,
 *  and <tt>Default</tt>. If <tt>INT_MAX</tt> or  <tt>INT_MIN</tt> is
 *  encountered, the value is set to the respective value as defined
 *  in the file limits.h. If the value <tt>Default</tt> is encountered,
 *  the #CurrValue is set to the default value, #DefaultVal,
 *  specified for the object.
 *
 *  If the default value specified for the object
 *  is set to <tt>NO_DEFAULT_INT</tt>, which is the initial value for
 *  #DefaultVal, then an error is thrown when <tt>default</tt>
 *  is encountered on the rhs.
 *
 *  When a match is found, the rhs of the card is processed into an integer
 *  using the function #process_LineEntry().
 *  The integer is then assigned to the #CurrValue of this object.
 *  And, it is optionally written out to an external address of an
 *  int, #AddrVal,
 *  that was supplied during the construction of the object.
 *  The object may then later be queried, using either the
 *  #currentTypedValue() or #currentValueAsVoidP()
 *  functions for the value of #CurrValue
 *  during subsequent processing.
 *
 *  A maximum and minimum allowed int value may be set as a requirement.
 *
 * <H2> Example of Usage </H2>
 *
 *   The example below sets up a required keyline entry, putting the value in
 *   the global address globInput.numIts.
 *
 *  @code
 *      struct globInput {
 *         int numIts;
 *      } gI;
 *
 *      BlockEntry *besmd = new BlockEntry("Iteration Options");
 *      int reqd = 1;
 *      LE_OneInt *i2 = new LE_OneInt("Number of Iterations",
 *                                     &gI.numIts,
 *                                     reqd, "numIts");
 *      i2->set_default(50);
 *      i2->set_limits(100000, 2);
 *      besmd->addLineEntry(i2);
 *  @endcode
 *
 *  An example of the input deck entry for this one LineEntry follows.
 *  Note, the line may be enclosed in any number of nested blocks,
 *  as long as the last nested block is named "Iteration Options".
 *  All white space differences and capitalization differences
 *  are ignored.
 *
 *  @code
 *       Start Block Iteration Options
 *         Number of Iterations = 150
 *       End   Block Iteration Options
 *  @endcode
 *
 *  It is an error for the Block "Iteration Options" to not
 *  have exactly one "Number of Iterations" %LineEntry's.
 *
 *  Examples of executing the code and extracting values follow.
 *
 *  The input file containing the entry may be processed via a
 *  command similar to the following
 *
 *  @code
 *     FILE *ifp = fopen("inputFile.txt", "r");
 *     besmd->read_block(ifp);
 *  @endcode
 *
 *  After running the above code,
 *  the value of <tt>gI.numIts</tt> is 150.
 *
 *  Also, we may alternatively query the <tt>besmd</tt> structure to find the
 *  the value in two ways:
 *
 *  @code
 *   LineEntry *le = besmd->searchLineEntry("Number of Iterations");
 *   int value = *(static_cast<const int *>(le->currentValueAsVoidP()));
 *
 *   LineEntry *le = besmd->searchLineEntry("Number of Iterations");
 *   LE_OneInt *le_int = dynamic_cast<LE_OneInt *>(le);
 *   int value = le_int->currentTypedValue();
 *  @endcode
 *
 *  Using this alternative, we could have set the external address used in
 *  the constructor of <B>%LE_OneInt</B> to 0.
 *
 *  <H2> Setting Dependencies on this card </H2>
 *
 *  Dependencies that this card may have on other cards in the input deck
 *  may be specified by calls to the member function #declareDependency().
 *  Dependencies that may be set on this card are listed below:
 *
 *       -  <b>BIDRT_PT_ERROR</b>
 *                 The target %BaseEntry must have occurred before this
 *                 card in the input deck.
 *       -  <b>BIDRT_ANTITHETICAL_ERROR</b>
 *                 The target %BaseEntry must not occur at all
 *                 in the input deck if this card occurs
 *       -  <b>BIDRT_ZERONUMTIMESREQUIRED</b>
 *                 If the dependency check on the target %BaseEntry
 *                 is satisfied, then the number of required
 *                 times this card is needed is turned to zero.
 *       -  <b>BIDRT_ONENUMTR</b>
 *                 If the dependency check on the target %BaseEntry
 *                 is satisfied, then the number of required
 *                 times this card is needed is turned to one.
 *       -  <b>BIDRT_USETOPROCESS </b>
 *                 The target %BaseEntry must have been already been
 *                 processed, and the default value for this card
 *                 is set to the value of the target %BaseEntry.
 *       -  <b>BIDRT_RTINTMM_ERROR</b>
 *                 If the #CurrValue of this entry is
 *                 in a specified range of values, then a
 *                 dependency check is made against the target
 *                 %BaseEntry.
 *
 *  This card may service dependency requests from other cards
 *  using the base service request, #ansDepCheck(), and the service request
 *  that returns an int value, #ansDepCheckOneInt().
 *  The int value returned is the
 *  #CurrValue of the %LE_OneInt object.
 *
 * @ingroup blockentryModule
 */
class LE_OneInt : public LineEntry
{
public:

    //! Main constructor for an int keyline entry.
    /*!
     *   This sets up the line entry special case.
     *   We make sure to call the base class constructor here to do
     *   much of the initialization.
     *
     *  When the keyline in the input file of the following form is
     *  found:
     *
     *  "KeyName" = [int]
     *
     * The int value at address, addrVal, is assigned the value
     * read in from the input file.
     *
     * @param keyName   C character string setting up the name
     *                  of the keyline to match
     * @param addrVal   Address of the int, external to the
     *                  object, which will get assigned the value of
     *                  the expression. (default 0)
     * @param numRL     Number of Required lines in the input file.
     *                  A fault is triggered if this number is nonzero
     *                  and the keyline isn't found in the input file.
     * @param varName   Variable name that is defined by this command.
     *                  This is only used for IO purposes.
     */
    LE_OneInt(const char* keyName, int* addrVal, int numRL = 0,
              const char* varName = 0);

    //! Copy constructor
    /*!
     * @param right Object to be copied
     */
    LE_OneInt(const LE_OneInt& right);

    //! Copy assignment operator
    /*!
     * @param right Object to be copied
     *
     * @return                      Returns a reference to the current object
     */
    LE_OneInt& operator=(const LE_OneInt& right);

    //! Duplicator function
    /*!
     * This function duplicates the entry and returns a pointer to a LineEntry
     *
     * @return                     Returns a pointer to the duplicated object
     */
    virtual LineEntry* duplMyselfAsLineEntry() const;

    //! Process this line Entry, which assigns a boolean variable
    /*!
     * This function is called when it has been determined that
     * the current KeyName from the input file matches this
     * object's KeyName.
     *
     * This function then processes the entry.
     * Processing involves checking for the satisfaction
     * of runtime dependencies.
     * Then, the CurrentValue is assigned and the external pointer
     * to int is assigned the CurrentValue.
     *
     * The NumTimesProcessed field is also incremented.
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
    virtual void print_usage(int ilvl) const;

    //! Adjust base address to store the value
    /*!
     * @param addrAdjustment Offset in raw bytes to adjust the
     *              internal value of the Address that will receive
     *              the int.
     */
    virtual void adjustAddress(LONG_PTR addrAdjustment);

    //! Return the current value as an int
    /*!
     * This is a nonvirtual function since the return type is specific to this child.
     *
     * @return                  Returns the current value as an int
     */
    int currentTypedValue() const;

    //! Return the current value as a const pointer to void
    /*!
     *  The calling function must know how to interpret the
     *  pointer to void. This is not very far fetched, since
     *  it should know that the value must be an int.
     *
     *  @return Returns a pointer to current int value as a pointer to void.
     */
    virtual const void* currentValueAsVoidP() const;

    //! Set the default value of this card
    /*!
     * @param defValue int value that will be the default.
     */
    void set_default(int defValue);

    //!  Set the maximum and minimum value
    /*!
     *  This puts bounds on the permissible entry values.
     *  If the bounds are violated, an Error is thrown.
     *
     * @param maxValue Maximum permissible value of the int
     * @param minValue Minimum permissible value of the int
     */
    void set_limits(int maxValue, int minValue);

    //! Sets the string that will be used as the name of
    //! the variable.
    /*!
     *  @param ps  string name of the variable.
     */
    void set_PrintString(const char* ps);

    //! Answer a query as to whether it can service a request of a
    //! certain type.
    /*!
     *  This call checks to whether requests can be fulfilled by this
     *  LineEntry object. Basically, this object can return wheether it
     *  has been called enough times, and it can return the integer
     *  value that was last processed.
     *
     * @param BIDSR_value Dependency service request type
     *      -> right now this is figured out from the previous two ints.
     *         It is not part of the interface.
     *
     * @return Returns true if the object can validly answer
     *         the query.
     */
    virtual bool DepCanService(BIDSR_TYPE BIDSR_value) const;

    //! Answer a query whether the object has been called before
    //! and what the current_value is
    /*!
     * This call returns a bool indicating whether it has been called before.
     * And if it has been, it returns the last value processed.
     *
     * @param returnInt  Int to be returned. This is the current value
     *
     * @return true if this object has been processed before.
     *              False otherwise
     */
    virtual bool ansDepCheckOneInt(int& returnInt) const;

private:

    //! The address of the int that gets updated by the
    //! processing of the %LineEntry
    /*!
     *  If this is zero, then no external update is carried out.
     */
    int* AddrVal;

    //! Max permissible value of an entry
    /*!
     * The default is set to INT_MAX
     */
    int MaxVal;

    //! Minimum permissible value of an entry
    /*!
     * The default is set to -INT_MAX
     */
    int MinVal;

    //! default value
    /*!
     * The default is set to NO_DEFAULT_INT
     */
    int DefaultVal;

    //! Current value
    /*!
     * This is initially set to NO_DEFAULT_INT
     */
    int CurrValue;

    //! Character string to be printed
    //! This represents the variable name
    char PrintString[MAX_INPUT_STR_LN+1];
};

}

#endif
