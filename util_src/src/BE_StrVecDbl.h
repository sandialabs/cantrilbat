/**
 * @file BE_StrVecDbl.h
 *   Header for the BlockEntry of a set of doubles that fills up a standard vector
 *   (see \ref blockentryModule and class \link BEInput::BE_StrVecDbl BE_StrVecDbl\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef BE_STRVECDBL_H
#define BE_STRVECDBL_H

#include "BE_BlockEntry.h"

namespace BEInput
{

//==================================================================================================================================
//!  %BlockEntry of a vector of related std::vector<double> values into a common  vector of vector<double>
/*!
 *   This sets up the <b>BlockEntry</b> special case for the entry of
 *    multiple double inputs into a common vector.  It declares interfaces
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
 *  The <I>BlockName</I> variable  may be consist of multiple
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
 *  If the default value specified for the object
 *  is set to <tt>NO_DEFAULT_DBL</tt>, which is the initial value for
 *  #DefaultVal, then an error is thrown when <tt>default</tt>
 *  is encountered on the rhs.
 *
 *  During preprocessing each of the possible character KeyNames
 *  generates it's own #BEInput::LE_OneDbl LineEntry object. The address of
 *  that is written to is equal to the base address plus the
 *  index of that particular character KeyName.
 *
 *  During #Wrapup() of the block, all of the individual doubles
 *  in the underlying LineEntries are collected at the block level.
 *  This is useful for later querying the results on the block level.
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
 *  A maximum and minimum allowed double value may be set as a requirement.
 *  on each double entry.
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
 *         std::vector<double> molarVolume;
 *         globInput() : nSpecies(3), molarVolume(3, 0.0) {};
 *      } gI;
 *
 *      BlockEntry *besmd = new BlockEntry("Fluid Properties");
 *      int reqd = 1;
 *      int subReqd = 1;
 *      int listLength = 3;
 *      char *charList[3] = {"H2O(l)", "K+", "Cl-"};
 *      BE_StrVecDbl *d2 = new BE_StrVecDbl("Molar Volume",  &gI.molarVolume, reqd, subReqd,
 *                                          charList, listLength, true, "molarVolume", besmd);
 *      d2->set_default(0.01);
 *      d2->set_limits(10., 1.0E-9);
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
 *         Start Block Molar Volume
 *            K+     = 0.001
 *            Cl-    = 0.002
 *            H2O(l) = 0.1
 *         End  Block Molar Volume
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
 *   BlockEntry *le = besmd->searchBlockEntry("Molar Volume");
 *   const std::vector<double> *molarVolume_ptr =
 *         *(static_cast<const double *>(be->currentValueAsVoidP()));
 *
 *   BlockEntry *be = besmd->searchBlockEntry("Molar Volume");
 *   BE_StrVecDbl *be_dbl = dynamic_cast<BE_StrDbl *>(be);
 *   const std::vector<double> * molarVolume_ptr = be_dbl->currentTypedValue();
 *  @endcode
 *
 *  Using this alternative, we could have set the external address used in
 *  the constructor of <B>%BE_StrDbl</B> to 0.
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
class BE_StrVecDbl : public BlockEntry
{
public:

    //! Main constructor for a double keyline entry.
    /*!
     *   This sets up the BlockEntry special case.
     *
     *  When the keyline in the input file of the following form is
     *  found:
     *
     *  "charListName = [double]"
     *
     * The double value at address, addrVal, is assigned the value
     * read in from the input file.
     *
     * @param[in]            blockName           C character string setting up the name of the block to match
     * @param[in]            ptrVecDbl           Address of the std::vector of doubles, external to the
     *                                           object, which will get assigned the value of
     *                                           the expressions. Defaults to  0, in which case the values
     *                                           are storred solely within the Block structure and must be retrieved
     *                                           from there one read in))
     * @param[in]            numTimesRequired    Number of Required blocks in the input file.
     *                                           A fault is triggered if this number is nonzero
     *                                           and the BlockName isn't found in the input file.
     * @param[in]            numSubLERequired    Number of times each sub LineEntry
     *                                           generated within the block is required. 0 or 1 allowed.
     * @param[in]            charList            Vector of C strings containing the character strings to match
     * @param[in]            listLength          Length of the charList vector.
     * @param[in]            constructLE         Boolean indicating whether to construct the
     *                                           individual Line entry commands (this must be true)
     * @param[in]            varName             Variable name that is defined by this command. This is only used for IO purposes.
     * @param[in]            parentBlock_input   Pointer to the parent block. Set to zero if this is not a parent block
     */
    BE_StrVecDbl(const char* blockName, std::vector<double>* ptrVecDbl,
                 int numTimesRequired, int numSubLERequired,
                 char** charList, int listLength, int constructLE,
                 const char* varName, BlockEntry* parentBlock_input = 0);

    //! Copy constructor
    /*!
     * @param[in]            right               Object to be copied
     */
    BE_StrVecDbl(const BE_StrVecDbl& right);

    //! Copy assignment operator
    /*!
     * @param right Object to be copied
     *
     * @return    Returns a reference to the current object
     */
    BE_StrVecDbl& operator=(const BE_StrVecDbl& right);

    //! Duplicator function
    /*!
     * This function duplicates the entry and returns a pointer  to a LineEntry
     *
     *  @return     Returns a pointer to the duplicate object
     */
    virtual BlockEntry* duplMyselfAsBlockEntry() const;

    //! Destructor
    ~BE_StrVecDbl();

    //! Virtual function called at the start of internally processing the block when it is 
    //! encountered within the input deck
    /*!
     *  This function may be used to start the process of setting up
     *  internal data functions when the block is first called within the input deck
     *  This is also where the current block processes the arguments specified on the START BLOCK line.
     *
     *  The default behavior is listed below.
     *
     *   -  An Error exit will occur if
     *      a blockArgTok string is actually supplied.
     *
     *   -  We increment the NumProcessedBlocks counter, which
     *      keeps track of the number of times this particular object has been called for parent block.
     *
     *   -  A dependency check is made to make sure that the
     *      required dependencies for this block have been satisfied.
     *
     *  Derived classes may override the default behavior. Usually
     *  derived classes will call the base class Initialization function and perhaps do some other processing.
     *  Here, we set the default value of the vector both in the Handle to the external vector and in the vector within this object.
     *
     *  @param[in]           ifp_input           File pointer to read additional keylines in the recursive calls to read_block
     *                                           Default = stdin
     *
     *  @param[in]           blockArgPtr         pointer to the TOKEN structure representing the argument list for the START BLOCK
     */
    void Initialization(FILE* ifp_input, const TK_TOKEN* blockArgPtr);

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
    void Wrapup(FILE* ifp_input, const TK_TOKEN* blockArgPtr);

    //! Adjust the address of external target objects
    /*!
     *  External addresses are adjusted by a raw byte value.
     *  Note, the input is a incremental adjustment. The cumulative
     *  adjustment is kept within the BlockEntry class object.
     *
     * @param adjustAAA  Raw bytes to adjust all target addressses
     *                   for objects within the block.
     *
     * @note this puts a severe restriction on the type of
     *       object to be used for the target.
     */
    virtual void adjustAddress(LONG_PTR adjustAAA);

    //! Return the current value as a pointer to a std::vector of doubles
    /*!
     * This is a nonvirtual function since the return type is specific to this child.
     *
     * @return               Returns a const pointer to the current vector of doubles
     */
    const std::vector<double>* currentTypedValue() const;

    //! Return the current value as a const pointer to void
    /*!
     *  The calling function must know how to interpret the pointer to void. This is not very far fetched, since
     *  it should know that the value must be a boolean.
     *
     *  @return                                  Returns a pointer to current std::vector of doubles as a pointer to void.
     */
    virtual const void* currentValueAsVoidP() const;

    //! Set the default value of the entries in the vector
    /*!
     *  If this is not explicitly called, the default entry for the vector is set at 0.0.
     *
     * @param[in]            defValue            Boolean value that will be the default.
     */
    virtual void set_default(double defValue);

    //!  Set the maximum and minimum value of an entry
    /*!
     *  An error will be thrown inf the keyline value isn't
     *  between these two.
     *
     * @param maxValue maximum value.
     *           The default value is DBL_MAX
     * @param minValue minimum value
     *           The default value is -DBL_MAX
     */
    void set_limits(double maxValue, double minValue);

    //! Sets the string that will be used as the name of
    //! the variable.
    /*!
     *  @param ps  string name of the variable.
     */
    void set_PrintString(const char* ps);

    //!   This subroutine will set up a default set of LineEntries for this
    //!  block. The limits and default value are inherited from the
    //!  block element object.
    /*!
     *  The type of LineEntries are LE_OneDbl
     *  All LineEntries are made optional
     */
    void generateDefLE();

protected:

    /**
     *  following is the handle to the double vector that
     *  will receive the input. In other words, this
     *  stores the address of the vector of doubles. It
     *  may be malloced within this routine. However, it
     *  is never owned and thus never freed by this routine.
     */
    //double **HndlDblVec;

    std::vector<double>* HndlVecDbl_;


    //! Max value that this number can attain
    double MaxVal;

    //! Min value that the entry can attain
    double MinVal;

    //! Default value
    double DefaultVal;

    //! Number of times each LineEntry generated by this
    //! block is required. ( 0 or 1).
    int m_numTimesRequiredLE;

    //! Vector of C strings representing entries in the vector
    /*!
     * This has a length equal to m_ListLength
     */
    char** CharList;

    //! Length of the length of the vector of doubles and the
    //! corresponding character list.
    int ListLength;

    //! Current list entry that was last set
    int CurrListValue;

    //! Current value
    /*!
     *  Initially this gets set to the default. Then when
     * the keyline is called, this gets set to the value of
     * the argument of the keyline.
     */
    double CurrValue;

    //! PrintString -> name of the variable that this keyLine update,
    //! as it will appear on descriptive printouts.
    /*!
     * This defaults to the name of the keyline
     */
    char PrintString[MAX_INPUT_STR_LN+1];

    //!   This is true if the default line Entries for this class
    //!  have already been generated. False, if they have yet to be
    //!   generated
    int defaultLE_made;

    //! Vector of the current values
    std::vector<double> currentVecValues_;

};
}
#endif
