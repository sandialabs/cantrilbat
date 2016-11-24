/**
 * @file BE_MultiBlockNested.cpp
 *  Definitions for the BlockEntry for a block which may occur multiple times within the input file
 *  (see \ref blockentryModule and class \link BEInput::BE_MultiBlockNested BE_MultiBlockNested\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "BE_MultiBlockNested.h"
#include "LE_LineEntry.h"
#include "mdp_allo.h"
#include "BE_UnitConversionLength.h"

using namespace mdpUtil;
//----------------------------------------------------------------------------------------------------------------------------------
namespace BEInput
{
//==================================================================================================================================
/*
 *   This sets up the line entry special case. We make sure to call the base class constructor here to do
 *   much of the initialization.
 */
BE_MultiBlockNested::BE_MultiBlockNested(const char* blockName, int* hndlnumstructures,
                                         int numTR, BlockEntry* parentBlock_input) :
    BlockEntry(blockName, numTR, parentBlock_input),
    m_numStructures(0),
    hndlNumStructures(hndlnumstructures),
    originalBlockPtr_(0)
{
    /*
     * Check to see HndlStructVec has been malloced ok
     */
    if (hndlNumStructures) {
        *hndlNumStructures = 0;
    }
    originalBlockPtr_ = this;
}
//==============================================================================================================
/*
 * BE_MultiBlockNested (const BE_MultiBlockNested&)
 */
BE_MultiBlockNested::BE_MultiBlockNested(const BE_MultiBlockNested& b) :
    BlockEntry(b),
    m_numStructures(b.m_numStructures),
    hndlNumStructures(b.hndlNumStructures),
    originalBlockPtr_(b.originalBlockPtr_)
{
}
//==================================================================================================================================
/*
 * BlockEntry* duplMyselfAsBlockEntry();     (virtual)
 *
 *  Duplicate the object in a list of BlockEntry objects
 */
BlockEntry* BE_MultiBlockNested::duplMyselfAsBlockEntry() const
{
    BE_MultiBlockNested* newBE = new  BE_MultiBlockNested(*this);
    return (BlockEntry*) newBE;
}
//==================================================================================================================================
/*
 * MultiBlockNested destructor: (virtual function)
 *
 * We malloced memory here, so we must explicitly call free.
 */
BE_MultiBlockNested::~BE_MultiBlockNested()
{
}
//==================================================================================================================================
/*
 *  Initialization() (virtual function)
 *
 *  This function is called when the block is seen in the input deck
 *  but before anything is done.
 */
void BE_MultiBlockNested::Initialization(FILE* ifp_input, const TK_TOKEN* blockArgTok)
{
    /*
     * Zero the linecounts of the blocks and lineEntries under the
     * Multiblock.
     */
    if (BlockLineInput) {
        for (int i = 0; i < numLineInput; i++) {
            LineEntry* le = BlockLineInput[i];
            le->zeroTimesCounter();
        }
    }
    if (SubBlocks) {
        for (int i = 0; i < numSubBlocks; i++) {
            BlockEntry* sb = SubBlocks[i];
            sb->ZeroLineCount();
        }
    }
    if (m_multiContribIndex != 0) {
        originalBlockPtr_->m_numTimesProcessed++;
    }

    BlockEntry::Initialization(ifp_input, blockArgTok);
}
//==================================================================================================================================
/*
 *  Wrapup() (virtual function)
 *
 *  This function is called when the end block line for the
 *  current block is read. Cleanup is done, and debugging printouts
 *  as well.
 */
void BE_MultiBlockNested::Wrapup(FILE* ifp_input, const TK_TOKEN* blockArgTok)
{
    // Expand the external structure list by one
    LONG_PTR addrDiff = expandStructListByOne();
    int mc = multiContribIndex();

    // Search for a block with the mc+1 multiContribIndex. If it exists, don't create it anew.
    int newmc = mc+1;
    BlockEntry* be = ParentBlock->match_block(&EntryName, newmc);
    if (!be) {
	// Create a copy of the current block
        BE_MultiBlockNested* newBlock = new BE_MultiBlockNested(*this);
	// Bump the multiContribIndex of the new block and add the new block into the existing structure
        newBlock->set_multiContribIndex(mc+1);
        newBlock->m_numTimesRequired = 0;
        newBlock->ZeroLineCount();
        ParentBlock->addSubBlock(newBlock);
        newBlock->adjustAddress(addrDiff);
    } else {
        be->ZeroLineCount();
    }
    // Call the underlying function
    BlockEntry::Wrapup(ifp_input, blockArgTok);
}
//==================================================================================================================================
/*
 * Expand the external structure by one (private function)
 */
LONG_PTR BE_MultiBlockNested::expandStructListByOne()
{

    /*
     * Create another entry in the structure list, null terminate the
     * list
     */
    if (m_numStructures <= m_multiContribIndex) {
        m_numStructures++;
	// Write the current value of multiContribIndex to the location specified by the calling program
        if (hndlNumStructures) {
            *hndlNumStructures = m_numStructures;
        }



        // Adjust the absolute addresses of the underlying entities
        //adjustAddress(addrDiff);
    }
    return (LONG_PTR) 0;
}
//==================================================================================================================================
int BE_MultiBlockNested::currentTypedValue() const
{
    return (m_multiContribIndex);
}
//==================================================================================================================================
/*
 * currentValueAsVoidP() (virtual)
 */
const void* BE_MultiBlockNested::currentValueAsVoidP() const
{
    return (const void*)(&m_multiContribIndex);
}
//==================================================================================================================================
/**************************************************************************/
}
