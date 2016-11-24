/**
 * @file BE_MultiBlockVec.hpp
 *   Definitions for the BlockEntry for a block which may occur multiple times within the input file
 *   whose line entry output can be collected in a simple structure
 *   (see \ref blockentryModule and class \link BEInput::BE_MultiBlockVec BE_MultiBlockVec\endlink).
 */

/*
 * Copywrite 2016 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef BE_MULTIBLOCKVEC_HPP
#define BE_MULTIBLOCKVEC_HPP

#include "BE_MultiBlockVec.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace BEInput
{
//==================================================================================================================================
template<typename T>
BE_MultiBlockVec<T>::BE_MultiBlockVec(const char* blockName, int* hndlnumstructures, std::vector<T*> *hndlStructVec,
                                      int numTimesRequired, BlockEntry* parentBlock_input) :
    BlockEntry(blockName, numTimesRequired, parentBlock_input),
    m_numStructures(0),
    hndlNumStructures(hndlnumstructures),
    HndlStructVec(hndlStructVec)
{
    /*
     * Check to see HndlStructVec has been malloced ok
     */
    if (HndlStructVec) {
        if ((*HndlStructVec).size() == 0) {
            (*HndlStructVec).resize(1, nullptr);
        }
        if ((*HndlStructVec)[0] == 0) {
            (*HndlStructVec)[0] = new T();
        }
    }
    if (hndlNumStructures) {
        *hndlNumStructures = 0;
    }
}
//==================================================================================================================================
template<typename T>
BE_MultiBlockVec<T>::BE_MultiBlockVec(const BE_MultiBlockVec<T>& b) :
    BlockEntry(b),
    m_numStructures(b.m_numStructures),
    hndlNumStructures(b.hndlNumStructures),
    HndlStructVec(b.HndlStructVec)
{
}
//==================================================================================================================================
template<typename T>
BlockEntry* BE_MultiBlockVec<T>::duplMyselfAsBlockEntry() const
{
    BE_MultiBlockVec<T>* newBE = new  BE_MultiBlockVec<T>(*this);
    return static_cast<BlockEntry*>(newBE);
}
//==================================================================================================================================
template<typename T>
BE_MultiBlockVec<T>::~BE_MultiBlockVec()
{
}
//==================================================================================================================================
/*
 *  This function is called when the block is seen in the input deck but before anything is done.
 */
template<typename T>
void BE_MultiBlockVec<T>::Initialization(FILE* ifp_input, const TK_TOKEN* blockArgTok)
{
    /*
     * Zero the linecounts of the blocks and lineEntries under the MultiBlockVec
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
    BlockEntry::Initialization(ifp_input, blockArgTok);
}
//==================================================================================================================================
/*
 *  This function is called when the end block line for the
 *  current block is read. Cleanup is done, and debugging printouts as well.
 */
template<typename T>
void BE_MultiBlockVec<T>::Wrapup(FILE* ifp_input, const TK_TOKEN* blockArgTok)
{
    // Expand the external structure list by one
    LONG_PTR addrDiff = expandStructListByOne();
    int mc = multiContribIndex();
    BlockEntry* be = ParentBlock->match_block(&EntryName, mc + 1);
    if (!be) {
        BE_MultiBlockVec<T>* newBlock = new BE_MultiBlockVec<T>(*this);
        newBlock->set_multiContribIndex(mc + 1);
        newBlock->ZeroLineCount();
        ParentBlock->addSubBlock(newBlock);
        newBlock->adjustAddress(addrDiff);
    } else {
        be->ZeroLineCount();
    }
    BlockEntry::Wrapup(ifp_input, blockArgTok);
}
//==================================================================================================================================
template<typename T>
LONG_PTR BE_MultiBlockVec<T>::expandStructListByOne()
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

        // Expand the list by one
        std::vector<T*>& StructList = *HndlStructVec;
        StructList.resize(m_numStructures + 1, nullptr);

        // Malloc a new structure onto the end of the list
        StructList[m_numStructures] = new T();

        LONG_PTR a1 = reinterpret_cast<LONG_PTR>(StructList[m_numStructures]);
        LONG_PTR a2 = reinterpret_cast<LONG_PTR>(StructList[m_numStructures -1]);
        // addrDiff is the amount that the absolute addresses of the underlying entities in the structure 
        // have to be adjusted compared to the last entry
        LONG_PTR addrDiff = a1 - a2;
        return addrDiff;
    }
    return (LONG_PTR) 0;
}
//==================================================================================================================================
template<typename T>
const T* BE_MultiBlockVec<T>::currentTypedValue() const
{
    if (!HndlStructVec) {
        return (const T*) nullptr;
    }
    return const_cast<const T*>((*HndlStructVec)[multiContribIndex()]);
}
//==================================================================================================================================
template<typename T>
const void* BE_MultiBlockVec<T>::currentValueAsVoidP() const
{
    if (!HndlStructVec) {
        return nullptr;
    }
    return (const void*)((*HndlStructVec)[multiContribIndex()]);
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

