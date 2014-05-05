/**
 * @file LE_MultiCStr.h
 *  Header for the LineEntry of multiple C strings.
 *  (see \ref blockentryModule and class 
 *  \link BEInput::LE_MultiCStr LE_MultiCStr\endlink).
 */
/*
 * $Author: hkmoffa $
 * $Revision: 5 $
 * $Date: 2012-02-23 14:34:18 -0700 (Thu, 23 Feb 2012) $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef LEMULTICSTR_H
#define LEMULTICSTR_H

#include "LE_LineEntry.h"

namespace BEInput {


  //!  This class stores character strings from 1 or more 
  //!  keylines. 
  /*!
   *   This sets up the <b>LineEntry</b> special case for the entry of 
   *   Multiple lines of text input.
   *   It declares interfaces for specifying the
   *   entry and for processing the entry once a match is
   *   made in the input deck. The general form of the input is
  
   *    \f[
   *               Keyname = \mathrm{Character String 1} 
   *               Keyname = \mathrm{Character String 2} 
   *    \f]
   *
   *  An arbitrary number of keylines may be entered.
   *  The resulting storage is in the form of a malloced array of character
   *  pointers. Each pointer is individually malloced. The array
   *  has an extra pointer entry that is zero.
   *  An example follows.
   *
   *   address = (char **) mdp_allo_ptr_1(3);
   *   char *address[0] = mdp_alloc_string("Character String1");
   *   char *address[1] = mdp_alloc_string("Character String2");
   *   char *address[2] = 0
   *
   *  The <I>KeyName</I> variable  may be consist of multiple
   *  tokens. It is tokenized and compared in a case insensitive format
   *  at the start of the processing and before any matching
   *  is carried out.
   *  \f$ \mathrm{int} \f$ is a single integer. The code will check that
   *  the rhs consists of well formed string.
   *
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
   *         char** Title;
   *      } gI;
   *
   *      BlockEntry *besmd = new BlockEntry("Run Description");
   *      int reqd = 0;
   *      LE_MultiCStr *t2 = new LE_MultiCStr("Title", &gI.Title,
   *                                     reqd, "Title");
   *      t2->set_default(50);
   *      t2->set_limits(100000, 2);
   *      besmd->addLineEntry(t2);
   *  @endcode
   *
   *  An example of the input deck entry for this one LineEntry follows.
   *  Note, the line may be enclosed in any number of nested blocks,
   *  as long as the last nested block is named "Run Description.
   *  All white space differences and capitalization differences
   *  are ignored.
   *
   *  @code
   *       Start Block Run Description
   *         Title = Analysis of sooting flame
   *         Title =     base case, parameter study strain rate = 2.0
   *         Title =     ethylene flame
   *       End   Block Run Description
   *  @endcode
   *
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
   *  the value of <tt>gI.Title </tt> is
   *   @verbatim
   *        Title[0] =  "Analysis of sooting flame"
   *        Title[1] =  "base case, parameter study strain rate = 2.0"
   *        Title[2] =  " ethylene flame"
   *        Title[3] =  \0
   *   @endverbatim
   *
   *  Also, we may alternatively query the <tt>besmd</tt> structure to find the
   *  the value in two ways:
   *
   *  @code
   *   LineEntry *le = besmd->searchLineEntry("Title");
   *   const char** value = *(static_cast<const char **>(le->currentValueAsVoidP()));
   * 
   *   LineEntry *le = besmd->searchLineEntry("Title");
   *   LE_OneMultiCStr *le_mc = dynamic_cast<LE_OneMultiCStr *>(le);
   *   const char ** value = le_mc->currentTypedValue();
   *  @endcode
   *
   *  Using this alternative, we could have set the external address used in
   *  the constructor of <B>%LE_MultiCStr</B> to 0.
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
   *  using the base service request, #ansDepCheck(). 
   *
   * @ingroup blockentryModule
   */
  class LE_MultiCStr : public LineEntry {
  public:

    //!  Constructor sets up the line entry.
    /*!
     *   This sets up the line entry special case. When the keystring
     *   is matched, it sets the value of the string to the contents
     *   of the rhs.
     *  
     *   The maximum number of tokens allowed in the string is equal
     *   to the limit specified in the tok_input_util.h file.
     *   The minimum number of tokens is equal to one. In other words we
     *   require something in the string as a default.
     *
     *   @param keyName Input keystring to be matched.
     *   @param addrVal  Address of the C ** string variable that 
     *                   will be used to store the result
     *   @param maxval   Maximum number of tokens allowed in the string
     *                   The default is MAXTOKENS
     *   @param minval   Minimum number of tokens allowed in the string
     *                   The default is 1
     *   @param numRL    Number or instances of the keystring that
     *                   must be present for the input deck to be 
     *                   complete. The default is zero.
     *   @param varName  Name of the variable to be used for describing
     *                   output errors.
     *
     *  Example
     * ---------
     *
     * Code Segment:
     *
     *    char *spName;
     *    LE_OneCStr *sN = new LE_OneSCtr("Name of the Species", &spName, 
     *                                     1, 1, 1, "SpeciesName");
     *    bl->addLineEntry(sN);
     *
     *    This sets up a line entry in the block entry, bl, that
     *    matches the keystring "Name of the Species". When found it
     *    sets the string entry to the value on the rhs of the equals
     *    sign. One and only one token is required on the rhs. Anything
     *    else is an error. And, this is a required entry. Meaning that
     *    if it is not found, and checking is turned on, then an
     *    exception is thrown.
     * 
     * Input Entry: 
     *            Name of the species = H2O(l)
     * Result:
     *            spname = "H2O(l)"
     *
     *  Input Entry
     *            Name of the Species = H2O (l)
     * Result
     *               -> exception thrown
     *
     */
    LE_MultiCStr(const char *keyName, char ***addrVal, int maxval = MAXTOKENS,
		 int minval = 1, int numRL = 0, const char *varName = 0);

    //! Destructor
    virtual ~LE_MultiCStr();

    //! Copy constructor
    /*!
     * @param right Object to be copied
     */
    LE_MultiCStr(const LE_MultiCStr& right);

    //! Copy assignment operator
    /*!
     * @param right Object to be copied
     */
    LE_MultiCStr& operator=(const LE_MultiCStr& right);

    //! Duplicator function
    /*!
     * This function duplicates the entry and returns a pointer
     * to a LineEntry
     */
    virtual  LineEntry* duplMyselfAsLineEntry() const;

    //! Process this line Entry, which assigns a std::string variable
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
    void process_LineEntry(const TK_TOKEN *lineArgTok);

    //!  This routine will print out a processed line
    /*!
     *  The default behavior is to print the original line with a "=>"
     *  prefix to indicate that action has been taken on it.
     *  Then, we print out a further line stating we have
     *  set a variable
     *   ====> PrintString = currValue
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
     *              the string.
     */
    virtual void adjustAddress(LONG_PTR addrAdjustment);

    //! Return the current Multilined C string input
    /*!
     * This is a nonvirtual function since the return type
     * is specific to this child.
     *   This returns a const representation of the
     *   value. Note, the member is presumed to be owned
     *   by this object.
     */
    const char **currentTypedValue() const;

    //! Return the current value as a const pointer to void
    /*!
     *  The calling function must know how to interpret the
     *  pointer to void. This is not very far fetched, since
     *  it should know that the value must be an int.
     *
     *  @return Returns a pointer to current CurrValue as a
     *          pointer to void. The CurrValue should be
     *          cast to a const char **. 
     */  
    virtual const void * currentValueAsVoidP() const;
    
    //! Set the default value of this card
    /*!
     * @param defValue C String value that will be the default.
     */
    void set_default(const char *defValue);
  
    //!   This sets a limit on the maximum and minimum number of tokens
    //!   allowed in a valid string.
    /*!
     * @param maxV  Maximum number of tokens. Default is a lot
     * @param minV  minimum number of tokens. Default is zero
     */
    void set_limits(int maxV, int minV);

    //! Sets the string that will be used as the name of 
    //! the variable.
    /*!
     *  @param ps  string name of the variable.
     */
    void set_PrintString(const char *ps);

  private:

    //! Following is the handle of the address of the
    //!  vector of pointers to c character strings that gets updated
    //!  by this command.
    /*!
     *    Note, we assume that *AddrVal is always a malloced value
     *    unless it is set to zero.
     */
    char ***AddrVal;

    //! Max permissible number of permissible tokens
    /*!
     * The default is set to INT_MAX
     */
    int MaxVal;

    //! Min permissible number of permissible tokens
    /*!
     * The default is set to 0
     */
    int MinVal;

    //! default value
    /*!
     * The default is set to the empty string
     */
    char *DefaultVal;

    //! Current value
    /*!
     * This is initially set to the empty string
     */
    char **CurrValue;

    //! Last current value
    /*!
     *  Points to the last C string input
     */
    char *LastCurrValue;

    //! Character string representing the variable name to be printed
    char PrintString[MAX_INPUT_STR_LN+1];
  };

}

#endif
