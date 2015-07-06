/**
 * @file BE_StrDbl.cpp
 *  Definitions for the BlockEntry of a set of doubles that fills up a vector
 *  (see \ref blockentryModule and class
 *  \link BEInput::BE_StrDbl BE_StrDbl\endlink).
 */
/*
 * $Author: hkmoffa $
 * $Revision: 508 $
 * $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "BE_StrDbl.h"
#include "LE_OneDbl.h"
#include "mdp_allo.h"

#include <cmath>
#include <cfloat>

using namespace mdpUtil;

namespace BEInput
{

/*
 * BE_StrDbl Constructor:
 *
 *   This sets up the line entry special case.
 *   We make sure to call the base class constructor here to do
 *   much of the initialization.
 */
BE_StrDbl::BE_StrDbl(const char* blockName, double** hndlAAA,
                     int numTimesRequired, int numSubLERequired,
                     char** charList, int ll, int constructLE,
                     const char* varName,  BlockEntry* parentBlock_input) :
    BlockEntry(blockName, numTimesRequired,
               parentBlock_input),
    HndlDblVec(hndlAAA),
    MaxVal(DBL_MAX),
    MinVal(-DBL_MAX),
    //DefaultVal(NO_DEFAULT_DOUBLE),
    DefaultVal(0.0),
    m_numTimesRequiredLE(numSubLERequired),
    CharList(0),
    ListLength(ll),
    CurrListValue(0),
    CurrValue(0.0),
    defaultLE_made(0),
    m_currentVecValues(0)
{
    char* item(0);
    /*
     * Set the limits and the default values
     */
    PrintString[0] = '\0';
    /*
     * Set up the block's character list
     */
    if (charList) {
        CharList =
            mdp_alloc_VecFixedStrings(ListLength, MAX_INPUT_STR_LN+1);
        for (int i = 0; i < ListLength; i++) {
            item = charList[i];
            if (!item) throw BI_InputError("BE_StrDbl::BE_StrDbl",
                                               "char list item is null");
            strncpy(CharList[i], item, MAX_INPUT_STR_LN);
            CharList[i][MAX_INPUT_STR_LN] = '\0';
        }
    }
    if (varName) {
        strncpy(PrintString, varName, MAX_INPUT_STR_LN);
    } else {
        strncpy(PrintString, blockName, MAX_INPUT_STR_LN);
    }
    PrintString[MAX_INPUT_STR_LN] = '\0';
    /*
     * Ok, now we are ready to instantiate the list of lineEntries for
     * this class
     */
    if (constructLE) {
        generateDefLE();
    }
    /*
     * Generate the space for the output if it hasn't already.
     * Note, in contrast to other mallocs, this is not owned
     * by class.
     */
    if (HndlDblVec) {
        if (*HndlDblVec == 0 && ListLength > 0) {
            *HndlDblVec = mdp_alloc_dbl_1(ListLength, DefaultVal);
        }
    }
    m_currentVecValues = mdp_alloc_dbl_1(ListLength, DefaultVal);
}

/*
 * BE_StrDbl(const BE_StrDbl&):
 *
 * copy constructor
 */
BE_StrDbl::BE_StrDbl(const BE_StrDbl& b) :
    BlockEntry(b),
    HndlDblVec(b.HndlDblVec),
    MaxVal(b.MaxVal),
    MinVal(b.MinVal),
    DefaultVal(b.DefaultVal),
    m_numTimesRequiredLE(b.m_numTimesRequiredLE),
    CharList(0),
    ListLength(b.ListLength),
    CurrListValue(b.CurrListValue),
    CurrValue(b.CurrValue),
    m_currentVecValues(0)
{
    if (ListLength > 0) {
        char* item;
        CharList =
            mdp_alloc_VecFixedStrings(ListLength, MAX_INPUT_STR_LN+1);
        for (int i = 0; i < ListLength; i++) {
            item = b.CharList[i];
            if (!item) throw BI_InputError("BE_StrDbl::BE_StrDbl",
                                               "char list item is null");
            strncpy(CharList[i], item, MAX_INPUT_STR_LN);
            CharList[i][MAX_INPUT_STR_LN] = '\0';
        }
        m_currentVecValues = mdp_alloc_dbl_1(ListLength, DefaultVal);
        mdp_copy_dbl_1(m_currentVecValues, b.m_currentVecValues, ListLength);
    }
    strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
}

/*
 * BE_StrDbl& operator=(const BE_StrDbl &b) :
 *
 *  assignment operator
 */
BE_StrDbl& BE_StrDbl::operator=(const BE_StrDbl& b)
{
    if (&b != this) {
        BlockEntry::operator=(b);
        HndlDblVec  = b.HndlDblVec;
        MaxVal     = b.MaxVal;
        MinVal     = b.MinVal;
        DefaultVal = b.DefaultVal;
        m_numTimesRequiredLE = b.m_numTimesRequiredLE;
        mdp_safe_free((void**) &CharList);
        ListLength = b.ListLength;
        CurrListValue = b.CurrListValue;
        CurrValue  = b.CurrValue;
        strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
        mdp_safe_free((void**) &m_currentVecValues);
        if (ListLength > 0) {
            char* item;
            CharList =
                mdp_alloc_VecFixedStrings(ListLength, MAX_INPUT_STR_LN+1);
            for (int i = 0; i < ListLength; i++) {
                item = b.CharList[i];
                if (!item) throw BI_InputError("LE_StrListDbl::LE_StrListDbl",
                                                   "char list item is null");
                strncpy(CharList[i], item, MAX_INPUT_STR_LN);
                CharList[i][MAX_INPUT_STR_LN] = '\0';
            }

            m_currentVecValues = mdp_alloc_dbl_1(ListLength, DefaultVal);
            mdp_copy_dbl_1(m_currentVecValues, b.m_currentVecValues, ListLength);
        }
    }
    return *this;
}
/*
 * BlockEntry* duplMyselfAsBlockEntry();  (virtual)
 *
 *  Duplicate yourself in a list of the base class
 */
BlockEntry* BE_StrDbl::duplMyselfAsBlockEntry() const
{
    BE_StrDbl* newBE = new BE_StrDbl(*this);
    return (BlockEntry*) newBE;
}

/*
 * BE_StrDbl destructor: (virtual function)
 *
 * We malloced memory here, so we must explicitly call free.
 */
BE_StrDbl::~BE_StrDbl()
{
#ifdef DEBUG_DESTRUCTOR
    printf("~BE_StrDbl called for %s\n", EntryName.orig_str);
#endif
    mdp_safe_free((void**) &CharList);
    mdp_safe_free((void**) &m_currentVecValues);
}

/*
 *  generateDefLE():
 *
 *   This subroutine will set up a default set of line entries for this
 *   block. The limits and default value are inherited from the
 *   block element object.
 */
void BE_StrDbl::generateDefLE()
{
    if (*HndlDblVec == 0) {
        *HndlDblVec = mdp_alloc_dbl_1(ListLength, 0.0);
    }
    double* AddrVal = *HndlDblVec;
    for (int i = 0; i < ListLength; i++) {
        LE_OneDbl* led =
            new LE_OneDbl(CharList[i], AddrVal + i,
                          m_numTimesRequiredLE, CharList[i]);
        if (DefaultVal != NO_DEFAULT_DOUBLE) {
            led->set_default(DefaultVal);
        }
        led->set_limits(MaxVal, MinVal);
        addLineEntry(led);
    }
    defaultLE_made = 1;
}

/*
 *  Initialization() (virtual function)
 *
 *  This function is called when the block is seen in the input deck
 *  but before anything is done.
 */
void BE_StrDbl::Initialization(FILE* ifp_input, const TK_TOKEN* blockArgTok)
{
    BlockEntry::Initialization(ifp_input, blockArgTok);
    /*
     *  Set the default values into the double array, if there is a
     *  default
     */
    if (DefaultVal != NO_DEFAULT_DOUBLE) {
        if (HndlDblVec) {
            double* AddrVal = *HndlDblVec;
            for (int i = 0; i < ListLength; i++) {
                AddrVal[i] = DefaultVal;
            }
        }
        for (int i = 0; i < ListLength; i++) {
            m_currentVecValues[i] = DefaultVal;
        }
    }
}

/*
 *  Wrapup() (virtual function)
 *
 *  This function is called when the end block line for the
 *  current block is read. Cleanup is done, and debugging printouts
 *  as well.
 */
void BE_StrDbl::Wrapup(FILE* ifp_input, const TK_TOKEN* blockArgTok)
{
    for (int i = 0; i < ListLength; i++) {
        LineEntry* LE = searchLineEntry(CharList[i]);
        LE_OneDbl* led = dynamic_cast<LE_OneDbl*>(LE);
        m_currentVecValues[i] = led->currentTypedValue();
    }

    /*
     *  Set the default values into the double array, if there is a
     *  default
     */
    if (BaseEntry::s_PrintProcessedLine && PrintString[0]) {
        printf("\t====> Summary of the block, %s:\n", EntryName.orig_str);
        for (int i = 0; i < ListLength; i++) {
            printf("\t====>   %s[ %s : %d ] = %g\n", PrintString, CharList[i],
                   i, m_currentVecValues[i]);
        }
        printf("\t====> End of Summmary of the block, %s:\n",
               EntryName.orig_str);
    }
    BlockEntry::Wrapup(ifp_input, blockArgTok);
}

/*
 * adjustAddress:
 *
 *  Note, the input is a incremental adjustment. The cumulative
 *  adjustment is kept in the class object.
 *
 */
void BE_StrDbl::adjustAddress(LONG_PTR adjustAAA)
{
    LONG_PTR pold = reinterpret_cast<LONG_PTR>(*HndlDblVec);
    LONG_PTR ll = reinterpret_cast<LONG_PTR>(HndlDblVec);
    ll += adjustAAA;
    HndlDblVec = reinterpret_cast<double**>(ll);
    LONG_PTR pnew = reinterpret_cast<LONG_PTR>(*HndlDblVec);
    LONG_PTR pdiff = pnew - pold;
    BlockEntry::adjustAddress(pdiff);
}

const double* BE_StrDbl::currentTypedValue() const
{
    return m_currentVecValues;
}

/*
 * currentValueAsVoidP() (virtual)
 */
const void* BE_StrDbl::currentValueAsVoidP() const
{
    return static_cast<const void*>(m_currentVecValues);
}

/*
 * set_default():
 *
 *   Set the Default value for entries in the list. Note, this value
 *  is imposed on all array elements when the block is initialized,
 *  if a default is specified.
 */
void BE_StrDbl::set_default(double def)
{
    DefaultVal = def;
    if (defaultLE_made) {
        LE_OneDbl* led;
        LineEntry* le;
        for (int i = 0; i < ListLength; i++) {
            le =  BlockLineInput[i];
            led = dynamic_cast<LE_OneDbl*>(le);
            if (led) {
                led->set_default(def);
            }
        }
    }
}

/*
 * set_limits()
 *
 */
void BE_StrDbl::set_limits(double maxV, double minV)
{
    MaxVal = maxV;
    MinVal = minV;
    if (defaultLE_made) {
        LE_OneDbl* led;
        LineEntry* le;
        for (int i = 0; i < ListLength; i++) {
            le =  BlockLineInput[i];
            led = dynamic_cast<LE_OneDbl*>(le);
            if (led) {
                led->set_limits(maxV, minV);
            }
        }
    }
}

/*
 * set_PrintString()
 *
 */
void BE_StrDbl::set_PrintString(const char* ps)
{
    if (ps) {
        strncpy(PrintString, ps, MAX_INPUT_STR_LN);
        PrintString[MAX_INPUT_STR_LN] = '\0';
    }
}
/**********************************************************************/
}
