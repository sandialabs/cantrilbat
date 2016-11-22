/**
 * @file BE_StrVecDbl.cpp
 *   Definitions for the BlockEntry of a set of doubles that fills up a standard vector
 *   (see \ref blockentryModule and class \link BEInput::BE_StrVecDbl BE_StrVecDbl\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "BE_StrVecDbl.h"
#include "LE_OneDbl.h"
#include "mdp_allo.h"

#include <cstring>
#include <cfloat>
#include <vector>

using namespace mdpUtil;
//----------------------------------------------------------------------------------------------------------------------------------
namespace BEInput
{
//===================================================================================================================================
/*
 * BE_StrDbl Constructor:
 *
 *   This sets up the line entry special case.
 *   We make sure to call the base class constructor here to do
 *   much of the initialization.
 */
BE_StrVecDbl::BE_StrVecDbl(const char* blockName, std::vector<double>* refVecDbl, int numTimesRequired, int numSubLERequired,
                           char** charList, int ll, int constructLE, const char* varName,  BlockEntry* parentBlock_input) :
    BlockEntry(blockName, numTimesRequired, parentBlock_input),
    HndlVecDbl_(refVecDbl),
    MaxVal(DBL_MAX),
    MinVal(-DBL_MAX),
    DefaultVal(0.0),
    m_numTimesRequiredLE(numSubLERequired),
    CharList(0),
    ListLength(ll),
    CurrListValue(0),
    CurrValue(0.0),
    defaultLE_made(0)
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
        CharList = mdp_alloc_VecFixedStrings(ListLength, MAX_INPUT_STR_LN+1);
        for (int i = 0; i < ListLength; i++) {
            item = charList[i];
            if (!item) {
                throw BI_InputError("BE_StrVecDbl::BE_StrVecDbl", "char list item is null");
            }
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
     * Ok, now we are ready to instantiate the list of lineEntries for this class
     */
    if (constructLE) {
        generateDefLE();
    }
    /*
     * Generate the space for the output if it hasn't already been allocated appropriately.
     * If not external posting has been specified, we do nothing.
     */
    if (HndlVecDbl_) {
        if (ListLength > 0) {
            if ((*HndlVecDbl_).size() < (size_t) ListLength) {
                (*HndlVecDbl_).resize(ListLength, DefaultVal);
            }
        }
    }
    currentVecValues_.resize(ListLength, DefaultVal);
}
//==================================================================================================================================
BE_StrVecDbl::BE_StrVecDbl(const BE_StrVecDbl& b) :
    BlockEntry(b),
    HndlVecDbl_(b.HndlVecDbl_),
    MaxVal(b.MaxVal),
    MinVal(b.MinVal),
    DefaultVal(b.DefaultVal),
    m_numTimesRequiredLE(b.m_numTimesRequiredLE),
    CharList(0),
    ListLength(b.ListLength),
    CurrListValue(b.CurrListValue),
    CurrValue(b.CurrValue)
{
    if (ListLength > 0) {
        char* item;
        CharList = mdp_alloc_VecFixedStrings(ListLength, MAX_INPUT_STR_LN+1);
        for (int i = 0; i < ListLength; i++) {
            item = b.CharList[i];
            if (!item) throw BI_InputError("BE_StrVecDbl::BE_StrVecDbl", "char list item is null");
            strncpy(CharList[i], item, MAX_INPUT_STR_LN);
            CharList[i][MAX_INPUT_STR_LN] = '\0';
        }
    }
    currentVecValues_ = b.currentVecValues_;
    strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
}
//==================================================================================================================================
BE_StrVecDbl& BE_StrVecDbl::operator=(const BE_StrVecDbl& b)
{
    if (&b != this) {
        BlockEntry::operator=(b);
        HndlVecDbl_  = b.HndlVecDbl_;
        MaxVal     = b.MaxVal;
        MinVal     = b.MinVal;
        DefaultVal = b.DefaultVal;
        m_numTimesRequiredLE = b.m_numTimesRequiredLE;
        mdp_safe_free((void**) &CharList);
        ListLength = b.ListLength;
        CurrListValue = b.CurrListValue;
        CurrValue  = b.CurrValue;
        strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
        if (ListLength > 0) {
            char* item;
            CharList = mdp_alloc_VecFixedStrings(ListLength, MAX_INPUT_STR_LN+1);
            for (int i = 0; i < ListLength; i++) {
                item = b.CharList[i];
                if (!item) throw BI_InputError("LE_StrListDbl::LE_StrListDbl", "char list item is null");
                strncpy(CharList[i], item, MAX_INPUT_STR_LN);
                CharList[i][MAX_INPUT_STR_LN] = '\0';
            }
        }
        currentVecValues_ = b.currentVecValues_;
    }
    return *this;
}
//==================================================================================================================================
BlockEntry* BE_StrVecDbl::duplMyselfAsBlockEntry() const
{
    BE_StrVecDbl* newBE = new BE_StrVecDbl(*this);
    return (BlockEntry*) newBE;
}
//==================================================================================================================================
BE_StrVecDbl::~BE_StrVecDbl()
{
#ifdef DEBUG_DESTRUCTOR
    printf("~BE_StrVecDbl called for %s\n", EntryName.orig_str);
#endif
    mdp_safe_free((void**) &CharList);
}
//==================================================================================================================================
/*
 *   This subroutine will set up a default set of line entries for this
 *   block. The limits and default value are inherited from the block element object.
 */
void BE_StrVecDbl::generateDefLE()
{
    double* AddrVal = &(currentVecValues_[0]);
    if (HndlVecDbl_) {
        if (ListLength > 0) {
            if ((*HndlVecDbl_).size() < (size_t) ListLength) {
                (*HndlVecDbl_).resize(ListLength, DefaultVal);
            }
        }
        AddrVal = &((*HndlVecDbl_)[0]);
    }
    for (int i = 0; i < ListLength; i++) {
        LE_OneDbl* led = new LE_OneDbl(CharList[i], AddrVal + i, m_numTimesRequiredLE, CharList[i]);
        if (DefaultVal != NO_DEFAULT_DOUBLE) {
            led->set_default(DefaultVal);
        }
        led->set_limits(MaxVal, MinVal);
        addLineEntry(led);
    }
    defaultLE_made = 1;
}
//==================================================================================================================================
/*
 *  This function is called when the block is seen in the input deck but before anything is done.
 */
void BE_StrVecDbl::Initialization(FILE* ifp_input, const TK_TOKEN* blockArgTok)
{
    BlockEntry::Initialization(ifp_input, blockArgTok);
    /*
     *  Set the default values into the double array, if there is a default
     */
    if (DefaultVal != NO_DEFAULT_DOUBLE) {
        if (HndlVecDbl_) {
            double* AddrVal = &((*HndlVecDbl_)[0]);
            for (int i = 0; i < ListLength; i++) {
                AddrVal[i] = DefaultVal;
            }
        }
        for (int i = 0; i < ListLength; i++) {
            currentVecValues_[i] = DefaultVal;
        }
    }
}
//==================================================================================================================================
/*
 *  This function is called when the end block line for the
 *  current block is read. Cleanup is done, and debugging printouts as well.
 */
void BE_StrVecDbl::Wrapup(FILE* ifp_input, const TK_TOKEN* blockArgTok)
{
    // Collect the line entries and fill up the currentVecValues_ vector
    for (int i = 0; i < ListLength; i++) {
        LineEntry* LE = searchLineEntry(CharList[i]);
        LE_OneDbl* led = dynamic_cast<LE_OneDbl*>(LE);
        currentVecValues_[i] = led->currentTypedValue();
    }
    /*
     *  Set the default values into the double array, if there is a  default
     */
    if (BaseEntry::s_PrintProcessedLine && PrintString[0]) {
        printf("\t====> Summary of the block, %s:\n", EntryName.orig_str);
        for (int i = 0; i < ListLength; i++) {
            printf("\t====>   %s[ %s : %d ] = %g\n", PrintString, CharList[i], i, currentVecValues_[i]);
        }
        printf("\t====> End of Summmary of the block, %s:\n", EntryName.orig_str);
    }
    BlockEntry::Wrapup(ifp_input, blockArgTok);
}
//==================================================================================================================================
/*
 *  Note, the input is a incremental adjustment. The cumulative
 *  adjustment is kept in the class object.
 */
void BE_StrVecDbl::adjustAddress(LONG_PTR adjustAAA)
{
    if (HndlVecDbl_) {
        LONG_PTR ll = reinterpret_cast<LONG_PTR>(HndlVecDbl_);
        ll += adjustAAA;
        HndlVecDbl_ = reinterpret_cast<std::vector<double> *>(ll);
    }
}
//==================================================================================================================================
const std::vector<double>* BE_StrVecDbl::currentTypedValue() const
{
    return & currentVecValues_;
}
//==================================================================================================================================
const void* BE_StrVecDbl::currentValueAsVoidP() const
{
    return static_cast<const void*>(& currentVecValues_);
}
//==================================================================================================================================
/*
 *   Set the Default value for entries in the list. Note, this value
 *   is imposed on all array elements when the block is initialized, if a default is specified.
 */
void BE_StrVecDbl::set_default(double def)
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
//==================================================================================================================================
void BE_StrVecDbl::set_limits(double maxV, double minV)
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
//==================================================================================================================================
void BE_StrVecDbl::set_PrintString(const char* ps)
{
    if (ps) {
        strncpy(PrintString, ps, MAX_INPUT_STR_LN);
        PrintString[MAX_INPUT_STR_LN] = '\0';
    }
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
