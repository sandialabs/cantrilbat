/**
 * @file BE_MultiBlockVec.h
 *   Declarations for the BlockEntry for a block which may occur multiple times within the input file
 *   (see \ref blockentryModule and class \link BEInput::BE_MultiBlockVec BE_MultiBlockVec\endlink).
 */

/*
 * Copywrite 2016 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef BEMULTIBLOCKVEC_H
#define BEMULTIBLOCKVEC_H

#include "BE_BlockEntry.h"

namespace BEInput
{

//!  BlockEntry of an arbitrary number of subblocks, putting the information into
//!  an external vector of pointers to storage structures
/*!
 *  All input from within a block can be sent to elements of a single structure,
 *  or queried from within the block structure after processing.
 *  This class will create a vector of "structures" for each
 *  block it encounters in the input file, mallocing additional
 *  elements of the structure as needed.
 *  It's kind of tricky, so I will attempt to explain fully.
 *  The structure, which is the templated type T, will be called storStruct in the discussion below.
 *
 *  It can be used to handle the input of an indefinate number
 *  of identical blocks in the input file. You set up the block
 *  in the same way as any other block in the input file. However,
 *  all of the line entries below the MultiblockVec must write to the
 *  same structure.
 *
 *  The vector of storStruct pointers will look like this in the calling program.
 *
 *         std::vector<storStruct*> storVec;
 *
 *  On return there will be a vector of filled structures that can
 *  be accessed via the following method. The ith structure may be
 *  accessed via:
 *
 *        storStruct* storPtr = storVec[i];
 *        double dataP = storPtr->dataP;
 *
 *  The class calculates an offset between one malloced storStruct
 *  structure and the next storStruct structure. It then adjusts the
 *  write addresses of all underlying line entries by the value of the offset (in a very raw manor).
 *  Therefore, subsequent writes will be directed to the next storStruct structure in the vector of storStruct
 *  pointers if the initial line entry write addresses all point to the initial element of the storStruct pointer vector.
 *  This later requirement, of course, is in the hands of the user, so the user should take care
 *  or all hell will break loose.
 *
 * <H2> Example of Usage </H2>
 *
 *  Example:
 *
 *   The example below sets up a required MultiBlock and keyline entry for entry of an arbitrary
 *   number of blocks
 *
 * @code
 *      struct globInput {
 *         int nSpecies;
 *         std::string pNames;
 *         int gasLiquid;
 *      }; 
 *      std::vector<globInput*> m_phases
 *      m_phases.push_back(new globInput());
 *      int numPhases;
 *
 *      BlockEntry *besmd = 
 *           new BE_MultiBlockVec<globInput>("Phase Definition", &numPhases, &m_phases, 1);
 *      rootBlock->addSubBlock(besmd);
 *      globInput* gI = m_phases[0];
 * 
 *      LE_OneInt* ins = new LE_OneInt("Number of Species", &(gI->nSpecies), 1, "nSpecies");
 *      ins->set_default(1);
 *      besmd->addLineEntry(ins);
 *   
 *      LE_OneStr* pName = new LE_OneStr("Phase Name",  &(gI->pNames), 15, 1, 0, "PhaseName");
 *      besmd->addLineEntry(pName);
 *
 *      LE_OneInt* ing = new LE_OneInt("Type of phase", &(gI->gasLiquid), 1, "GasLiquid");
 *      ing->set_default(1);
 *      besmd->addLineEntry(ing);
 *                                           
 * @endcode
 *
 *  An example of the input deck entry for this MultiBlockVec  follows.
 *  Note, the line may be enclosed in any number of nested blocks,
 *  as long as the last nested block is named "Fluid Properties".
 *  All white space differences and capitalization differences are ignored.
 *
 *  @code
 *    Start Block rootBlock
 *       Start Block Phase Definition
 *           Phase Name = WaterVapor
 *           Number of Species = 3
 *           Type of phase = 0
 *       End Block Phase Definition
 *
 *       Start Block Phase Definition
 *           Phase Name = Brine
 *           Number of Species = 33
 *           Type of phase = 1
 *       End Block   Phase Definition
 *
 *       Start Block Phase Definition
 *           Phase Name = CuO
 *           Number of Species = 1
 *           Type of phase = 2
 *       End Block   Phase Definition
 *    End Block rootBlock
 *  @endcode
 *
 *  It is an error for the Block "rootBlock"  to have zero "Phase Definition" blocks.
 *  However, it may have an arbitrary number of them as long as there is one.
 *  In this case, on exit, m_mphases will have a length of 3 with the globInput struct
 *  filled up with the appropriate inputs from the three blocks.
 *
 *  Alternatively, the user can process the rootBlock structure to extract all of the
 *  information input as well. There will be a total of 4 MultiBlockVec blocks
 *  defined in that structure, with the first three filled up with user information,
 *  and the last block empty, with a read count of zero.
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
 *  This card may service dependency requests from other cards using the base service request, #ansDepCheck().
 *
 * @ingroup blockentryModule
 */
template<typename T>
class BE_MultiBlockVec : public BlockEntry
{
public:

    //! Constructor for the BE_MultiBlockVec<T> class.
    /*!
     *
     * @param                blockName           Block name of the block. Note there can be multiple blocks with the same name.
     *
     * @param                hndlnumstructures   This is a pointer to the location where the number of multiblocks processed will be
     *                                           stored. It is an optional entry, and not necessarily needed.
     *                                           However, it's convenient. 
     *
     * @param                hndlStructVec       This is the handle to the vector of pointers to storStruct structures, T. 
     *                                           Typically, using the example above, you would use (&storVec) as the argument here.
     *
     *
     * @param                numTimesRequired    Number of times this block is required.
     *
     * @param                parentBlock_input   Pointer to the parent block. Set to zero if there is no parent block,
     *                                           or it will be hooked up later
     */
    BE_MultiBlockVec(const char* blockName, int* hndlnumstructures, std::vector<T*>* hndlStructVec, int numTimesRequired,
                     BlockEntry* parentBlock_input = 0);

    //! Copy constructor
    /*!
     * @param right Object to be copied
     */
    BE_MultiBlockVec(const BE_MultiBlockVec<T>& right);

    //! Duplicator function
    /*!
     *   This class operates on  multiple input blocks both by duplicating itself and by
     *   "duplicating" the structure that the block writes to.
     *
     *   @return                                 Returns a pointer to the duplicate object.
     */
    virtual BlockEntry* duplMyselfAsBlockEntry() const;

    //! Destructor for the MultiBlock class.
    /*!
     *   Note, the class does not own the vector of storStruct structures, as is the paradigm
     *   throughout the block input utility. The calling program must free that structure.
     *   Therefore, there is no memory to free.
     */
    ~BE_MultiBlockVec();

    //! Virtual function called at the start of internally processing the block
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
    void Initialization(FILE* ifp_input, const TK_TOKEN* blockArgTok);


    //! Virtual function called at the start of internally processing the block
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
     *                    in the recursive calls to read_block  Default = stdin
     *
     *  @param blockArgTok pointer to the TOKEN structure representing
     *                     the argument list for the START BLOCK
     */
    virtual void Wrapup(FILE* ifp_input, const TK_TOKEN* blockArgTok);

private:
    //! Expand the underlying external structure list by one
    /*!
     *   Then, it malloces both a new  vector of storStruct pointers and then malloces a new
     *  storStruct pointer itself by calling m_fRetnNewStruct(),
     *  and then adjusting all of the underlying line element
     *  addresses to write into locations in that new structure.
     *  Therefore, when a new block is encounted in the input deck
     *  the BE_Multiblock object will be ready.
     *  This function may be used to start the process of setting up
     *  internal data functions when the block is first called.
     *
     *    @return                                Returns the difference between the old pointer and the new pointer
     */
    LONG_PTR expandStructListByOne();

public:

    //! Return the current value as a pointer to T, the external storage structure
    /*!
     * This is a nonvirtual function since the return type
     * is specific to this child. The current block points to one external storage structure
     */
    const T* currentTypedValue() const;

    //! Return the current value as a const pointer to void. It's a pointer to T, 
    //! the external storage structure
    /*!
     *  The calling function must know how to interpret the pointer to void. This is not very far fetched, since
     *  it should know that the value must be a boolean.
     *
     *  @return                                  Returns a pointer to current storage structure
     */
    virtual const void* currentValueAsVoidP() const;

    /********************************************************************/
    /*           MEMBER DATA                                            */
    /********************************************************************/

private:
    //!  Number of multblocks processed.
    int m_numStructures;

    /**
     *  This is a pointer to the location where the number of multiblocks processed will be
     *  stored. It is an optional entry, and not necessarily needed since the vector of pointers to storage structures 
     *  has its own size dimension. However, it's a convenience.
     */
    int* hndlNumStructures;

    /**
     *  This is the handle to the vector of pointers to storStruct structures, T. This may be zero, in which case
     *  all input must be obtained from the BlockEntry structure at the end of processing.
     */
    std::vector<T*>* HndlStructVec;

};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
