/**
 * @file BE_MolalityComp.cpp
 *  Definition for the BlockEntry of a set of doubles that fills up a vector of molalities
 *  (see \ref blockentryModule and class  \link BEInput::BE_MolalityComp BE_MolalityComp\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "BE_MolalityComp.h"
#include "LE_OneDbl.h"

#include <cmath>
#include <cstdlib>
//----------------------------------------------------------------------------------------------------------------------------------
namespace BEInput
{
//==================================================================================================================================
BE_MolalityComp::BE_MolalityComp(const char* blockName, double** hndlAAA, int numTimesRequired,
                                 char** charList, int listLength, int indexSolvent, double mwSolvent,
                                 const char* varName, BlockEntry* parentBlock_input) :
    BE_StrDbl(blockName, hndlAAA, numTimesRequired, 0, charList,listLength, 1, varName, parentBlock_input),
    m_indexSolvent(indexSolvent),
    m_mwSolvent(mwSolvent)
{
    BE_StrDbl::set_default(0.0);
    set_limits(1.0E10, 0.0);
    /*
     * Fill in the default value for the solvent molality.
     * Note by the return from BE_StrDbl(), the space for the
     * storage of the molalities has been malloced.
     */
    double* addr = *HndlDblVec;
    double defV =  1000. / m_mwSolvent;
    addr[m_indexSolvent] = defV;

    LineEntry* le =  BlockLineInput[indexSolvent];
    LE_OneDbl* led = dynamic_cast<LE_OneDbl*>(le);
    if (led) {
        led->set_default(defV);
    } else {
        printf("Didn't set default value\n");
        exit(-1);
    }
    DefaultVal = NO_DEFAULT_DOUBLE;
}
//==================================================================================================================================
BE_MolalityComp::BE_MolalityComp(const char* blockName, double* const fixedAddr, int numTimesRequired,
                                 char** charList, int listLength, int indexSolvent, double mwSolvent,
                                 const char* varName, BlockEntry* parentBlock_input) :
    BE_StrDbl(blockName, fixedAddr, numTimesRequired, 0, charList,listLength, 1, varName, parentBlock_input),
    m_indexSolvent(indexSolvent),
    m_mwSolvent(mwSolvent)
{
    BE_StrDbl::set_default(0.0);
    set_limits(1.0E10, 0.0);
    /*
     * Fill in the default value for the solvent molality.
     * Note by the return from BE_StrDbl(), the space for the
     * storage of the molalities has been malloced.
     */
    double defV =  1000. / m_mwSolvent;
    if (HndlDblVec) {
      double* addr = *HndlDblVec;
      addr[m_indexSolvent] = defV;
    }

    LineEntry* le = BlockLineInput[indexSolvent];
    LE_OneDbl* led = dynamic_cast<LE_OneDbl*>(le);
    if (led) {
        led->set_default(defV);
    } else {
        printf("Didn't set default value\n");
        exit(-1);
    }
    DefaultVal = NO_DEFAULT_DOUBLE;
}
//==================================================================================================================================
BE_MolalityComp::BE_MolalityComp(const BE_MolalityComp& b) :
    BE_StrDbl(b),
    m_indexSolvent(b.m_indexSolvent),
    m_mwSolvent(b.m_mwSolvent)
{
}
//===============================================================================================================================
BE_MolalityComp& BE_MolalityComp::operator=(const BE_MolalityComp& right)
{
    if (&right != this) {
        BE_StrDbl::operator=(right);
        m_indexSolvent = right.m_indexSolvent;
        m_mwSolvent    = right.m_mwSolvent;
    }
    return *this;
}
//====================================================================================================================================
BlockEntry* BE_MolalityComp::duplMyselfAsBlockEntry() const
{
    BE_MolalityComp* newBE = new BE_MolalityComp(*this);
    return (BlockEntry*) newBE;
}
//==================================================================================================================================
BE_MolalityComp::~BE_MolalityComp()
{
#ifdef DEBUG_DESTRUCTOR
    printf("~BE_MolalityComp called for %s\n", BlockName.orig_str);
#endif
}
//==================================================================================================================================
void BE_MolalityComp::set_default(double def)
{
    if (defaultLE_made) {
        LE_OneDbl* led;
        LineEntry* le;
        for (int i = 0; i < ListLength; i++) {
            if (i != m_indexSolvent) {
                le =  BlockLineInput[i];
                led = dynamic_cast<LE_OneDbl*>(le);
                if (led) {
                    led->set_default(def);
                }
            }
        }
    }
    double defV = 1000. / m_mwSolvent;
    if (HndlDblVec) {
        double* addr = *HndlDblVec;
        addr[m_indexSolvent] = defV;
    }
    LineEntry* le = BlockLineInput[m_indexSolvent];
    LE_OneDbl* led = dynamic_cast<LE_OneDbl*>(le);
    if (led) {
        led->set_default(defV);
    } else {
        printf("Didn't set default value\n");
        exit(-1);
    }
    DefaultVal = NO_DEFAULT_DOUBLE;
}
//==================================================================================================================================
/*
 *  This function is called when the end block line for the
 *  current block is read. Cleanup is done, and debugging printouts as well.
 */
void BE_MolalityComp::Wrapup(FILE* ifp_input, const TK_TOKEN* blockArgTok)
{
    BE_StrDbl::Wrapup(ifp_input, blockArgTok);
    double correct = 1000. / m_mwSolvent;
    double* AddrVal =  m_currentVecValues;
    double mSol = AddrVal[m_indexSolvent];
    if (fabs(correct - mSol) > 1.0E-7) {
        throw BI_InputError("BE_MolalityComp::Wrapup", "Solvent molality has been incorrectly overwritten");
    }
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
