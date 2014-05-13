/**
 * @file BE_MoleComp.h
 *   Header for the BlockEntry of a set of doubles that fills up a vector
 *   of mole fractions
 *  (see \ref blockentryModule and class 
 *  \link BEInput::BE_MoleComp BE_MoleComp\endlink).
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef BE_MOLECOMP_H
#define BE_MOLECOMP_H

#include "BE_StrDbl.h"

namespace BEInput {


  //!  BlockEntry of a vector of mole fractions (or mass fractions)
  //!  into a vector of doubles, where normalization of the sum of mole
  //!  fractions to one occurs after the end block statement.
  /*!
   *   This sets up the <b>BlockEntry</b> special case for the entry of
   *   multiple
   *   double inputs into a common vector.  It declares interfaces 
   *   for specifying the
   *   entry and for processing the entry once a match is
   *   made in the input deck. The general form of the input is
   *
   *   \f[
   *    \begin{array}[t]{l}
   *      \mathrm{Start\ block\ BlockName} \\
   *         \mathrm{\ \ \ KeyName1\ =\ }\mathrm{double} \\
   *         \mathrm{\ \ \ KeyName2\ =\ }\mathrm{double} \\
   *         \mathrm{\ \ \ \ \ \ \ \ \ ...  }  \\
   *       \mathrm{End\ block\ BlockName} 
   *     \end{array}
   *   \f]
   *
   *  The <I>BlockName</I> variable may be consist of multiple
   *  tokens. It is tokenized and compared in a case insensitive format
   *  at the start of the processing and before any matching
   *  is carried out.
   *   \f$ \mathrm{KeyName1} \f$  \f$ \mathrm{KeyName1} \f$ are character
   *  Keynames generated from a permissible list of length \f$ ListLength\f$.
   *  They may also consist of multiple tokens, and they are 
   *  compared on a case insensitive format.
   *
   *  \f$ \mathrm{double} \f$ represents a single double.
   *  The code will check that
   *  the rhs consists of single tokens.
   *  The token must be a well formed double value.
   *  The following tokens are understood and accepted in
   *  place of \f$ \mathrm{double} \f$ : <tt>DBL_MAX</tt>,  <tt>DBL_MIN</tt>,
   *  and <tt>Default</tt>. If <tt>DBL_MAX</tt> or  <tt>DBL_MIN</tt> is
   *  encountered, the value is set to the respective value as defined
   *  in the file math.h. If the value <tt>Default</tt> is encountered,
   *  the #CurrValue is set to the default value, #DefaultVal,
   *  specified for the object.
   *
   *  The default value specified for the object 
   *  is set to zero, which is the initial value for
   *  #DefaultVal. Therefore, any entry not specified will have
   *  a zero mole fraction assigned to it.
   *
   *  At the wrapup stage of the block processing. All of the
   *  entries within the block are added up. The final results
   *  are normalized to a value of 1.0. If the entries fall 
   *  outside of this range sum < 0.98 or sun > 1.02, an error is thrown.
   *
   *  This routine inherits from BE_StrDbl object, and most 
   *  of the functionality is implemented there.
   *
   *  These doubles are assigned to the m_CurrentVecValues[iList]
   *  entry and #CurrValue of this object.
   *  And, it is optionally written out to an external address, as
   *  (*HndlDblVecAddr)[iList],
   *  that was supplied during the construction of the object.
   *  The object may then later be queried, using either the
   *  #currentTypedValue() or #currentValueAsVoidP()
   *  functions for the value of #CurrValue and m_CurrentVecValues[i]
   *  during subsequent processing.
   *
   *  A maximum and minimum allowed double value is automatically
   *  set to 1.0 and 0.0 respectively.
   *
   * <H2> Example of Usage </H2>
   *
   *
   *  Example:
   *
   *   The example below sets up a required keyline entry, putting the value in
   *   the global address globInput.molarVolume
   *
   * @code
   *      struct globInput {
   *         int nSpecies;
   *         double *moleFractions;
   *         globInput() : nSpecies(3), molarVolume(0) {};
   *      } gI;
   *
   *      BlockEntry *besmd = new BlockEntry("Fluid Properties");
   *      int reqd = 1;
   *      int subReqd = 1;
   *      int listLength = 3;
   *      char *charList[3] = {"H2", "O2", "N2"};
   *      BE_MoleComp *d2 = new BE_MoleComp("Gas Mole Fractions",
   *                                     &gI.moleFractions, reqd,
   *                                     charList, listLength, true,
   *                                     "moleFractions", besmd);
   * @endcode
   *
   *  An example of the input deck entry for this BlockEntry follows.
   *  Note, the line may be enclosed in any number of nested blocks,
   *  as long as the last nested block is named "Fluid Properties".
   *  All white space differences and capitalization differences
   *  are ignored.
   *
   *  @code
   *       Start Block Fluid Properties
   *         Start Block Gas Mole Fractions
   *            H2     = 0.001
   *            O2     = 0.23
   *            N2     = 0.77
   *         End  Block Gas Mole Fractions
   *       End Block   Fluid Properties
   *  @endcode
   *
   *  It is an error for the Block "Fluid Properties" to not
   *  have exactly one "Fluid Properties" sub block. And, within
   *  that block it is required that every possible LineEntry
   *  generated actually be present. 
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
   *  
   *  After running the above code, 
   *  the value of <tt>gI.molarVolume[]</tt> is [0.1, 0.001, 0.002].
   * 
   *  Also, we may alternatively query the <tt>besmd</tt> structure to find the
   *  the value in two ways:
   *
   *  @code
   *   BlockEntry *le = besmd->searchBlockEntry("Gas Mole Fractions");
   *   const double *moleFractions = 
   *         *(static_cast<const double *>(be->currentValueAsVoidP()));
   * 
   *   BlockEntry *be = besmd->searchBlockEntry("Gas Mole Fractions");
   *   BE_MoleComp *be_dbl = dynamic_cast<BE_MoleComp *>(be);
   *   const double * moleFractions = be_dbl->currentTypedValue();
   *  @endcode
   *
   *  Using this alternative, we could have set the external address used in
   *  the constructor of <B>%BE_MoleComp</B> to 0.
   *
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
  class BE_MoleComp : public BE_StrDbl {
  public:

    //! Main constructor for a double keyline entry.
    /*!
     *   This sets up the line entry special case for mole fraction
     *   entry.
     *   We make sure to call the base class constructor here to do
     *   much of the initialization. 
     *
     *  When the keyline in the input file of the following form is
     *  found:
     *
     *  "charListName = [double]"
     *
     * The double value at address, addrVal, is assigned the value
     * read in from the input file.
     *
     * @param blockName   C character string setting up the name
     *                  of the block to match
     * @param hndlAddr   Address of the vector of doubles, external to the 
     *                  object, which will get assigned the value of 
     *                  the expressions. (default 0)
     * @param numTimesRequired Number of Required blocks in the input file.
     *                  A fault is triggered if this number is nonzero
     *                  and the BlockName isn't found in the input file.
     * @param charList  Vector of C strings containing the character 
     *                  strings to match
     * @param listLength Length of the charList vector.
     * @param constructLE Boolean indicating whether to construct the
     *                  Individual Line entry commands (this must be true)
     * @param varName   Variable name that is defined by this command.
     *                  This is only used for IO purposes.
     * @param parentBlock_input Pointer to the parent block. Set to
     *                 zero if this is no parent block
     */
    BE_MoleComp(const char *blockName, double **hndlAddr, 
		int numTimesRequired,
		char **charList, int listLength, int constructLE,
		const char *varName,
		BlockEntry *parentBlock_input = 0);

    //! Copy constructor
    /*!
     * @param right Object to be copied
     */
    BE_MoleComp(const BE_MoleComp& right);

    //! Copy assignment operator
    /*!
     * @param right Object to be copied
     */
    BE_MoleComp& operator=(const BE_MoleComp& right);

    //! Duplicator function
    /*!
     * This function duplicates the entry and returns a pointer
     * to a LineEntry
     */
    virtual BlockEntry* duplMyselfAsBlockEntry() const;  

    //! Destructor
    ~BE_MoleComp();

    //! This virtual function is used to wrap up the setup of the current
    //! block before returning to the parent block.
    /*!
     *
     *  The default behavior is listed below.
     *
     *   -  An Error exit will occur if 
     *      a blockArgTok string is actually supplied.
     *
     *   -  A checkRequirements check is made on all enclosed line
     *      entry objects. If they haven't been satisfied, an error
     *      is generated. 
     *
     *   -  A checkRequirements check is made on all enclosed blocks.
     *      If they haven't been satisfied, an error is generated. 
     *
     * Derived classes may override the default behavior. Usually
     * derived classes will call the base class Initialization function
     * and perhaps do some other processing.
     *
     *  @param ifp_input  File pointer to read additional keylines
     *                    in the recursive calls to read_block
     *                    Default = stdin
     *
     *  @param blockArgPtr pointer to the TOKEN structure representing
     *                     the argument list for the START BLOCK
     */
    void Wrapup(FILE *ifp_input, const TK_TOKEN *blockArgPtr);

  };
}
#endif
