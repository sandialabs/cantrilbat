/**
 * @file LE_OneBoolInt.h
 *  Header for the LineEntry of a single bool
 *  (see \ref blockentryModule and class 
 *     \link BEInput::LE_OneBoolInt LE_OneBoolInt\endlink).
 */
/* $Author: hkmoffa $
 * $Revision: 5 $
 * $Date: 2012-02-23 14:34:18 -0700 (Thu, 23 Feb 2012) $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef LEONEBOOLINT_H
#define LEONEBOOLINT_H

#include "LE_LineEntry.h"

namespace BEInput {

  //!  Keyline Entry of a single bool value, storred in an int
  /*!
   *   This sets up the <b>LineEntry</b> special case for the entry of a single
   *   bool input. It declares interfaces for specifying the
   *   entry and for processing the entry once a match is
   *   made in the input deck. The general form of the input is
   *
   *   \f[
   *     KeyName = \mathrm{bool}
   *   \f]
   *
   *  The <I>KeyName</I> variable  may be consist of multiple
   *  tokens. It is tokenized and compared in a case insensitive format
   *  at the start of the processing and before any matching
   *  is carried out.
   *  \f$ \mathrm{bool} \f$ is a single bool. The code will check that
   *  the rhs consists of one and only one well-formed bool.
   *
   *   The following tokens are accepted for true:
   *   TRUE, yes, T, Y
   *
   *   The following tokens are accepted for false:
   *   no, false, n, f
   *
   *  Additionally, the following tokens are understood and accepted in
   *  place of \f$ \mathrm{bool} \f$ :  <tt>Default</tt>.
   *  If the value <tt>Default</tt> is encountered,
   *  the #CurrValue is set to the default value, #DefaultVal,
   *  specified for the object.
   *
   *  If the default value specified for the object 
   *  is set to <tt>NO_DEFAULT_INT</tt>, which is the initial value for
   *  #DefaultVal, then an error is thrown when <tt>default</tt> 
   *  is encountered on the rhs.
   *
   *  When a match is found, the rhs of the card is processed into a bool
   *  using the function #process_LineEntry().
   *  The bool is then assigned to the #CurrValue of this object.
   *  And, it is optionally written out to an external address of an 
   *  int, #AddrVal,
   *  that was supplied during the construction of the object.
   *  The object may then later be queried, using either the
   *  #currentTypedValue() or #currentValueAsVoidP()
   *  functions for the value of #CurrValue
   *  during subsequent processing.
   *
   * <H2> Example of Usage </H2>
   *
   *   The example below sets up a required keyline entry, putting the value in
   *   the global address globInput.doExtraWork
   *
   * @code
   *      struct globInput {
   *         int doExtraWork;
   *      } gI;
   *
   *      BlockEntry *besmd= new BlockEntry("Iteration Options");
   *      int reqd = 1;
   *      LE_OneBoolInt *i2 = new LE_OneBoolInt("Do Extra Work",
   *                                            &gI.doExtraWork,
   *                                            reqd, "doExtraWork");
   *      i2->set_default(true);
   *      besmd->addLineEntry(i2);
   * @endcode
   *
   *  An example of the input deck entry for this one LineEntry follows.
   *  Note, the line may be enclosed in any number of nested blocks,
   *  as long as the last nested block is named "Iteration Options".
   *  All white space differences and capitalization differences
   *  are ignored.
   *
   *   @code
   *       Start Block Iteration Options
   *          Do Extra Work = false
   *       End   Block Iteration Options
   *  @endcode
   *
   *  After running the above code, 
   *  the value of <tt>gI.doExtraWork</tt> is false.
   *
   *  Also, we may alternatively query the <tt>besmd</tt> structure to find the
   *  the value in two ways:
   *
   *  @code
   *   LineEntry *le = besmd->searchLineEntry("Do Extra Work");
   *   int value = *(static_cast<const int *>(le->currentValueAsVoidP()));
   * 
   *   LineEntry *le = besmd->searchLineEntry("Do Extra Work");
   *   LE_OneBool *le_bool = dynamic_cast<LE_OneBool *>(le);
   *   int value = le_bool->currentTypedValue();
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
   *  that
   *  returns an int value, #ansDepCheckOneInt(). The int value returned is the
   *  #CurrValue of the %LE_OneBool object, translated to an int.
   *
   * @ingroup blockentryModule
   *
   */
  class LE_OneBoolInt : public LineEntry {
  public:

    //! Main constructor for a boolean-as-int keyline entry.
    /*!
     *   This sets up the line entry special case.
     *   We make sure to call the base class constructor here to do
     *   much of the initialization. 
     *
     *  When the keyline in the input file of the following form is
     *  found:
     *
     *  "KeyName" = [boolean]
     *
     * The boolean value at address, addrVal, is assigned the value
     * read in from the input file.
     *
     * @param keyName  C Character string setting up the name
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
    LE_OneBoolInt(const char *keyName, int *addrVal, int numRL = 0,
		  const char *varName = 0);

    //! Copy constructor
    LE_OneBoolInt(const LE_OneBoolInt&);

    //! Copy assignment operator
    LE_OneBoolInt& operator=(const LE_OneBoolInt&);

    //! Duplicator function
    /*!
     * This function duplicates the entry and returns a pointer
     * to a LineEntry
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
     * to boolean is assigned the CurrentValue.
     * 
     * The NumTimesProcessed field is also incremented.
     *
     * @param lineArgTok  Pointer to the token containing the
     *                    arguments to the lineEntry. This is 
     *                    everything after the "=" sign.
     */
    void process_LineEntry(const TK_TOKEN *lineArgTok);

    //!  This routine will print out a processed line
    /*!
     *  The default behavior is to print the original line with a "====>"
     *  prefix to indicate that action has been taken on it.
     *
     * @param lineArgTok  Pointer to the token containing the
     *                    arguments to the lineEntry. This is 
     *                    everything after the "=" sign.
     */
    void print_ProcessedLine(const TK_TOKEN *lineArgTok) const;

    //! Print out API information about this keyline
    /*!
     * This routine will print out to stdout information about
     * the keyline. This command is used to document the interface.
     *
     * @param ilvl Level of the indentation to use in printing
     */
    void print_usage(int ilvl) const;

    //! Adjust base address to store the value
    /*!
     * @param addrAdjustment Offset in raw bytes to adjust the
     *              internal value of the Address that will receive
     *              the boolean.
     */
    void adjustAddress(LONG_PTR addrAdjustment);

    //! Return the current value as a int
    /*!
     * This is a nonvirtual function since the return type
     * is specific to this child.
     */
    int currentTypedValue() const;

    //! Return the current value as a const pointer to void
    /*!
     *  The calling function must know how to interpret the
     *  pointer to void. This is not very far fetched, since
     *  it should know that the value must be a boolean.
     *
     *  @return Returns a pointer to current boolean value as a
     *          pointer to void.  
     */
    virtual const void * currentValueAsVoidP() const;

    //! Set the default value of this card
    /*!
     * @param defValue int value that will be the default.
     */
    void set_default(int defValue);

    //! Sets the string that will be used as the name of 
    //! the variable.
    /*!
     *  @param ps  string name of the variable.
     */
    void set_PrintString(const char *ps);

  private:
  
    //!The address of the boolean that gets updated by this command
    /*!
     *  If this is zero, then no external update is carried out.
     */
    int *AddrVal;

    //! Default value
    int DefaultVal;

    //! Current value
    /*!
     *  Initially this gets set to the default. Then when
     *  the keyline is called, this gets set to the value of 
     *  the argument of the keyline.
     */
    int CurrValue; 

    //! PrintString -> name of the variable that this keyLine update,
    //! as it will appear on descriptive printouts.
    /*!
     * This defaults to the name of the keyline
     */
    char PrintString[MAX_INPUT_STR_LN+1];
  };

}
#endif
