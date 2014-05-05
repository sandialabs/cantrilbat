/**
 * @file BE_MultiBlock.cpp
 *  Definitions for the BlockEntry for a block which may occur
 *  multiple times within the input file.
 *  (see \ref blockentryModule and class 
 *  \link BEInput::BE_MultiBlock BE_MultiBlock\endlink).
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

#include "BE_MultiBlock.h"
#include "mdp_allo.h"
#include "LE_OneDblUnits.h"

using namespace mdpUtil;

namespace BEInput {


/**********************************************************************
 *
 * BE_MultiBlock Constructor:
 *
 *   This sets up the line entry special case.
 *   We make sure to call the base class constructor here to do
 *   much of the initialization. 
 */
BE_MultiBlock::BE_MultiBlock(const char *blockName, 
			     int *hndlnumstructures,
			     void ***hndlStructVec,
			     void*  (*fptr)(void *),
			     void* function_data_loc,
			     int numTR,
			     BlockEntry *parentBlock_input) :
    BlockEntry(blockName, numTR, parentBlock_input),
    m_numStructures(0),
    hndlNumStructures(hndlnumstructures),
    HndlStructVec(hndlStructVec),
    m_fRetnNewStruct(fptr),
    m_function_data_loc(function_data_loc)
{
    /*
     * Check to see HndlStructVec has been malloced ok
     */
    if (HndlStructVec) {
	if (*HndlStructVec == 0) {
	    *HndlStructVec = mdp_alloc_ptr_1(2);
	}
	if (!fptr) {
	    throw BI_InputError("BE_MultiBlock::BE_MultiBlock",
				" Function pointer is null");
	}
	if ((*HndlStructVec)[0] == 0) {
	    (*HndlStructVec)[0] = m_fRetnNewStruct(m_function_data_loc);
	}
    }
    if (hndlNumStructures) {
	*hndlNumStructures = 0;
    }
}
//=================================================================================================
/*
 * BE_MultiBlock (const BE_MultiBlock&)
 */
BE_MultiBlock::BE_MultiBlock(const BE_MultiBlock &b) :
    BlockEntry(b),
    m_numStructures(b.m_numStructures),
    hndlNumStructures(b.hndlNumStructures),
    HndlStructVec(b.HndlStructVec),
    m_fRetnNewStruct(b.m_fRetnNewStruct),
    m_function_data_loc(b.m_function_data_loc)
{

}
//===================================================================================================
  /*
   * BlockEntry* duplMyselfAsBlockEntry();     (virtual)
   *
   *  Duplicate the object in a list of BlockEntry objects
   */  
  BlockEntry* BE_MultiBlock::duplMyselfAsBlockEntry() const {
    BE_MultiBlock* newBE = new  BE_MultiBlock(*this);
    return (BlockEntry *) newBE;
  }

  /*
   * MultiBlock destructor: (virtual function)
   *
   * We malloced memory here, so we must explicitly call free.
   */
  BE_MultiBlock::~BE_MultiBlock()
  {
  }

  /*
   *  Initialization() (virtual function)
   *
   *  This function is called when the block is seen in the input deck
   *  but before anything is done.
   */
  void BE_MultiBlock::Initialization(FILE *ifp_input, 
				     const TK_TOKEN *blockArgTok)
  {
    /*
     * Zero the linecounts of the blocks and lineEntries under the
     * Multiblock.
     */
    if (BlockLineInput) {
      for (int i = 0; i < numLineInput; i++) {
        LineEntry *le = BlockLineInput[i];
	le->zeroTimesCounter();
      }
    }
    if (SubBlocks) {
      for (int i = 0; i < numSubBlocks; i++) {
        BlockEntry *sb = SubBlocks[i];
        sb->ZeroLineCount();
      }
    }

    BlockEntry::Initialization(ifp_input, blockArgTok); 
  }

  /*
   *  Wrapup() (virtual function)
   *
   *  This function is called when the end block line for the
   *  current block is read. Cleanup is done, and debugging printouts
   *  as well.
   */
  void BE_MultiBlock::Wrapup(FILE *ifp_input, const TK_TOKEN *blockArgTok)
  {
    // Expand the external structure list by one
    LONG_PTR addrDiff = expandStructListByOne();
    int mc = multiContribIndex();

    int newmc = mc+1;
    BlockEntry *be = ParentBlock->match_block(&EntryName, newmc);
    if (!be) {
      BE_MultiBlock *newBlock = new BE_MultiBlock (*this);
      
      newBlock->set_multiContribIndex(mc+1);
      newBlock->ZeroLineCount();
      ParentBlock->addSubBlock(newBlock);
      newBlock->adjustAddress(addrDiff);
    } else {
      be->ZeroLineCount();
    }
    // Call the underlying function
    BlockEntry::Wrapup(ifp_input, blockArgTok);
  }

  /*
   * Expand the external structure by one (private function)
   */
  LONG_PTR BE_MultiBlock::expandStructListByOne()
  {
    // quick return if there are no external structures
    if (! HndlStructVec) {
      return (LONG_PTR) 0;
    }
    
    /*
     * Create another entry in the structure list, null terminate the
     * list
     */
    if (m_numStructures <= m_multiContribIndex) {
      m_numStructures++;
      if (hndlNumStructures) {
	*hndlNumStructures = m_numStructures;
      }

      // Expand the list
      mdp_realloc_ptr_1(HndlStructVec, m_numStructures + 2,
			m_numStructures);
    

      // get the pointer to the vector of structures
      void **StructList = *HndlStructVec;
      // Malloc a new structure onto the end of the list
      StructList[m_numStructures] = m_fRetnNewStruct(m_function_data_loc);
      // Make sure the list is null terminated
      StructList[m_numStructures + 1] = 0;
      LONG_PTR a1 = reinterpret_cast<LONG_PTR>(StructList[m_numStructures]);
      LONG_PTR a2 = reinterpret_cast<LONG_PTR>(StructList[m_numStructures -1]);
      LONG_PTR addrDiff = a1 - a2;
      // Adjust the absolute addresses of the underlying entities
      //adjustAddress(addrDiff);
      return addrDiff;
    }
    return (LONG_PTR) 0;
  }


  const void ** BE_MultiBlock::currentTypedValue() const {
    return const_cast<const void **>(*HndlStructVec);
  }

  /*
   * currentValueAsVoidP() (virtual)
   */
  const void * BE_MultiBlock::currentValueAsVoidP() const {
    return static_cast<const void *>(*HndlStructVec);
  }

  /**************************************************************************/
}
