/**
 * @file BE_MultiBlockNested.h
 *  Header for the BlockEntry for a block which may occur
 *  multiple times within the input file.
 *  (see \ref blockentryModule and class 
 *  \link BEInput::BE_MultiBlockNested BE_MultiBlockNested\endlink).
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

#ifndef BEMULTIBLOCKNESTED_H
#define BEMULTIBLOCKNESTED_H

#include "BlockEntry.h"
#include "stdio.h"
#include "mdp_allo.h"

namespace BEInput {

  // EXPERIMENTAL CLASS

  //!  This class sets up the case where there may be multiple
  //!  blocks of the same type.
  /*!
   *  All input from within a block
   *  must be sent to elements of a single structure.
   *  This class will create a vector of "structures" for each
   *  block it encounters in the input file, mallocing additional
   *  elements of the structure as  needed.
   *  It's kind of tricky, so I will attempt to explain fully.
   *  The structure will be called storStruct in the discussion 
   *  below.
   *
   *  It can be used to handle the input of an indefinate number
   *  of identical blocks in the input file. You set up the block
   *  in the same way as any other block in the input file. However,
   *  all of the line entries below the Multiblock must write to the
   *  same structure. 
   *
   *  The vector of storStruct pointers will look like this in the
   *  calling program.
   *         struct storStruct **storVec;
   *
   *  On return there will be a vector of filled structures that can
   *  be accessed via the following method. The ith structure may be
   *  accessed via:
   *
   *        struct storStruct *storPtr = storVec[i];
   *        double dataP = storPtr->dataP;
   * 
   *  The class calculates an offset between one malloced storStruct
   *  structure and the next storStruct structure. It then adjusts the
   *  write addresses of all underlying line entries by the 
   *  value of the offset (in a very raw manor). 
   *  Therefore, subsequent writes will be directed
   *  to the next storStruct structure in the vector of storStruct
   *  structures if the initial line entry write addresses all 
   *  point to the initial element of the storStruct vector.
   *  This later requirement, of course, is in the hands of the user,
   *  so the user should take care.
   *
   *
   * Example:
   * 
   *
   */

  //!  BlockEntry of a vector of related double values into a common
  //!  vector of doubles
  /*!
   *   This sets up the <b>BlockEntry</b> special case for the entry of
   *    multiple
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
   *  the CurrValue is set to the default value, DefaultVal,
   *  specified for the object.
   *
   *  If the default value specified for the object 
   *  is set to <tt>NO_DEFAULT_DBL</tt>, which is the initial value for
   *  DefaultVal, then an error is thrown when <tt>default</tt> 
   *  is encountered on the rhs.
   *
   *  During preprocessing each of the possible character KeyNames
   *  generates it's own #BEInput::LE_OneDbl LineEntry object. The address of
   *  that is written to is equal to the base address plus the
   *  index of that particular character KeyName.
   *
   *  During #Wrapup() of the block, all of the individual doubles
   *  in the underlying LineEntries are collected at the block level.
   *  This is useful for later querying the results on the block
   *  level. 
   *
   *  These doubles are assigned to the m_CurrentVecValues[iList]
   *  entry and CurrValue of this object.
   *  And, it is optionally written out to an external address, as
   *  (*HndlDblVecAddr)[iList],
   *  that was supplied during the construction of the object.
   *  The object may then later be queried, using either the
   *  #currentTypedValue() or #currentValueAsVoidP()
   *  functions for the value of CurrValue and m_CurrentVecValues[i]
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
   *         double *molarVolume;
   *         globInput() : nSpecies(3), molarVolume(0) {};
   *      } gI;
   *
   *      BlockEntry *besmd = new BlockEntry("Fluid Properties");
   *      int reqd = 1;
   *      int subReqd = 1;
   *      int listLength = 3;
   *      char *charList[3] = {"H2O(l)", "K+", "Cl-"};
   *      BE_StrDbl *d2 = new BE_StrDbl("Molar Volume",
   *                                     &gI.molarVolume, reqd, subReqd,
   *                                     charList, listLength, true,
   *                                     "molarVolume", besmd);
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
   *   const double *molarVolume = 
   *         *(dynamic_cast<const double *>(be->currentValueAsVoidP()));
   * 
   *   BlockEntry *be = besmd->searchBlockEntry("Molar Volume");
   *   BE_StrDbl *be_dbl = dynamic_cast<BE_StrDbl *>(be);
   *   const double * molarVolume = be_dbl->currentTypedValue();
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
   *                 If the CurrValue of this entry is
   *                 in a specified range of values, then a
   *                 dependency check is made against the target
   *                 %BaseEntry.
   *
   *  This card may service dependency requests from other cards
   *  using the base service request, #ansDepCheck().
   *
   * @ingroup blockentryModule
   */
  class BE_MultiBlockNested : public BlockEntry {
  public:

    //! Constructor for the BE_Multiblock class.
    /*!  
     *
     * @param blockName Block name of the block. Note there can
     *                  be multiple blocks with the same name.
     *
     * @param hndlnumstructures This is a pointer to the location
     *          where the number of multiblocks processed will be
     *          storred. It is an optional entry, and not necessarily
     *          needed since the vector of pointers to storage structures
     *          is null terminated. However, the last storStruct 
     *          structure in the vector will be unfilled with data.
     *          So, it's convenient to have this entry.
     *
     * @param hndlStructVec This is the handle to the vector
     *         of pointers to storStruct structures. Typically, using the
     *         example above, you would use (&storVec) as the argument
     *         here.
     *
     * @param fptr This is a function pointer to a function that 
     *             returns (pointer to void) and takes a 
     *             (pointer to void) as its one argument.
     *             It creates a new malloced structure and returns
     *             a pointer to the structure. for example:
     *   
     *           void * newStorStruct(void *fdata) {
     *              void *ptr =  new storStruct(fdata);
     *              return (ptr);
     *           }
     * @param function_data_loc This is the data location of anything
     *               needed in the constructor call of the new operation
     *               in the function pointed to by fptr.
     *
     * @param numTR Number of times this block is required.
     *
     * @param parentBlock_input Pointer to the parent block. Set to
     *                 zero if this is no parent block
     *
     */
    BE_MultiBlockNested(const char *blockName,
		  int *hndlnumstructures,
		  int numTR,
		  BlockEntry *parentBlock_input = 0);

    //! Copy constructor
    /*!
     * @param right Object to be copied
     */
    BE_MultiBlockNested(const BE_MultiBlockNested& right);

    //! Duplicator function
    /*!
     *   This class operates on 
     *   multiple input blocks not by duplicating itself, but by
     *   "duplicating" the structure that the block writes to.
     */
    virtual BlockEntry* duplMyselfAsBlockEntry() const;

    //! Destructor for the MultiBlock class.
     /*! 
      *   Note, the class does not
      *    own the vector of storStruct structures, as is the paradigm
      *    throughout the block input utility. The calling program must
      *    free that structure.
      */
    ~BE_MultiBlockNested();

    //! Virtual function called at the start of internally processing
    //! the block
    /*!
     *  This function may be used to start the process of setting up
     *  internal data functions when the block is first called. 
     *  This is also where the current block processes the
     *  arguments specified on the START BLOCK line. 
     *
     *  The default behavior is listed below.
     *
     *   -  An Error exit will occur if 
     *      a blockARgTok string is actually supplied.
     *
     *   -  We increment the NumProcessedBlocks counter, which
     *      keeps track of the number of times this particular object has
     *      been called for parent block.
     *
     *   -  A dependency check is made to make sure that the 
     *      required dependencies for this block have been satisfied.
     *
     * Derived classes may override the default behavior. Usually
     * derived classes will call the base class Initialization function
     * and perhaps do some other processing.
     *
     *  @param ifp_input  File pointer to read additional keylines
     *                    in the recursive calls to read_block
     *                    Default = stdin
     *
     *  @param blockArgTok pointer to the TOKEN structure representing
     *                     the argument list for the START BLOCK
     */
    void Initialization(FILE *ifp_input, 
			const TK_TOKEN *blockArgTok);

   
    //! Virtual function called at the start of internally processing
    //! the block
    /*!
     *  Virtual program that supercedes the BlockEntry::Wrapup program.
     *  Actually, this is where a lot of the special work within
     *  BE_Multiblock is done. This routine checks to see that the current
     *  block input was correct. Then, it malloces both a new
     *  vector of storStruct pointers and then malloces a new
     *  storStruct pointer itself by calling m_fRetnNewStruct(),
     *  and then adjusting all of the underlying line element
     *  addresses to write into locations in that new structure.
     *  Therefore, when a new block is encounted in the input deck
     *  the BE_Multiblock object will be ready.
     *  This function may be used to start the process of setting up
     *  internal data functions when the block is first called. 
     *  This is also where the current block processes the
     *  arguments specified on the START BLOCK line. 
     *
     *  The default behavior is listed below.
     *
     *   -  An Error exit will occur if 
     *      a blockARgTok string is actually supplied.
     *
     *   -  We increment the NumProcessedBlocks counter, which
     *      keeps track of the number of times this particular object has
     *      been called for parent block.
     *
     *   -  A dependency check is made to make sure that the 
     *      required dependencies for this block have been satisfied.
     *
     * Derived classes may override the default behavior. Usually
     * derived classes will call the base class Initialization function
     * and perhaps do some other processing.
     *
     *  @param ifp_input  File pointer to read additional keylines
     *                    in the recursive calls to read_block
     *                    Default = stdin
     *
     *  @param blockArgTok pointer to the TOKEN structure representing
     *                     the argument list for the START BLOCK
     */
    virtual void Wrapup(FILE *ifp_input, const TK_TOKEN *blockArgTok);

    //! Expand the underlying external structure list by one
    /*!
     *   Then, it malloces both a new
     *  vector of storStruct pointers and then malloces a new
     *  storStruct pointer itself by calling m_fRetnNewStruct(),
     *  and then adjusting all of the underlying line element
     *  addresses to write into locations in that new structure.
     *  Therefore, when a new block is encounted in the input deck
     *  the BE_Multiblock object will be ready.
     *  This function may be used to start the process of setting up
     *  internal data functions when the block is first called.
     */
  private:
    LONG_PTR expandStructListByOne();

  public:
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

    /********************************************************************/
    /*           MEMBER DATA                                            */
    /********************************************************************/

  private:
    /**
     *  Number of multblocks processed.
     */
    int m_numStructures;


    /**
     *  This is a pointer to the location
     *  where the number of multiblocks processed will be
     *  storred. It is an optional entry, and not necessarily
     *  needed since the vector of pointers to storage structures
     *  is null terminated. However, the last storStruct 
     *  structure in the vector will be unfilled with data.
     *  So, it's convenient to have this entry.
     */
    int *hndlNumStructures;

    //! pointer to the original block
    BE_MultiBlockNested *originalBlockPtr_;

    
  };

}

#endif
