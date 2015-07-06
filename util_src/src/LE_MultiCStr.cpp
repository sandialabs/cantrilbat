/**
 * @file LE_MultiCStr.cpp
 *  Definitions for the LineEntry of multiple C strings.
 *  (see \ref blockentryModule and class
 *  \link BEInput::LE_MultiCStr LE_MultiCStr\endlink).
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

#include "LE_MultiCStr.h"
#include "mdp_allo.h"

#include <climits>

using namespace mdpUtil;

namespace BEInput
{

//====================================================================================================================
/*
 * LE_MultiCStr Constructor:
 *
 *   This sets up the line entry special case.
 *   We make sure to call the base class constructor here to do
 *   much of the initialization.
 *
 *   The maximum number of tokens allowed in the string is equal
 *   to the limit specified in the tok_input_util.h file.
 *   The minimum number of tokens is equal to one. In other words we
 *   require something in the string as a default.
 *
 *   aaa is the handle for the resulting array of character strings
 *   to be returned by this function.
 *   *aaa points to a malloced array of pointers to c character strings.
 *   This array is always null terminated.
 */
LE_MultiCStr::LE_MultiCStr(const char* lineName, char** *aaa, int maxval,
                           int minval, int numRL, const char* varName) :
    LineEntry(lineName, numRL),
    AddrVal(aaa),
    MaxVal(maxval),
    MinVal(minval),
    DefaultVal(0),
    CurrValue(0),
    LastCurrValue(0)
{
    /*
     * Set the limits and the default values
     */
    CurrValue = (char**) mdp_alloc_ptr_1(1);
    PrintString[0] = '\0';
    if (AddrVal) {
        if (*AddrVal == 0) {
            *AddrVal = (char**)mdp_alloc_ptr_1(1);
        } else {
            printf("Error *aaa for LE_MultiCStr is non-zero at start\n");
        }
    }
    if (varName) {
        strncpy(PrintString, varName, MAX_INPUT_STR_LN);
    } else {
        strncpy(PrintString, lineName, MAX_INPUT_STR_LN);
    }
    PrintString[MAX_INPUT_STR_LN] = '\0';
}
//====================================================================================================================
/*
 * LE_MultiCStr(const LE_MultiCStr&):
 *
 * copy Constructor:
 */
LE_MultiCStr::LE_MultiCStr(const LE_MultiCStr& b) :
    LineEntry(b),
    AddrVal(b.AddrVal),
    MaxVal(b.MaxVal),
    MinVal(b.MinVal)
{
    DefaultVal = mdp_copy_string(b.DefaultVal);
    CurrValue = (char**) mdp_alloc_ptr_1(m_numTimesProcessed+1);
    for (int i = 0; i < m_numTimesProcessed; i++) {
        CurrValue[i]  = mdp_copy_string(b.CurrValue[i]);
    }
    if (m_numTimesProcessed == 0) {
        LastCurrValue = 0;
    } else {
        LastCurrValue = CurrValue[m_numTimesProcessed-1];
    }
    strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
}
//====================================================================================================================
/*
 * ~LE_MultiCStr():
 *
 *  Destructor eliminates the space possibly malloced by this structure
 *  and owned by this structure. Note, no space under AddrVal is
 *  owned by this structure by definition.
 */
LE_MultiCStr::~LE_MultiCStr()
{
    mdp_safe_free((void**)&DefaultVal);
    if (CurrValue) {
        int i = 0;
        while (CurrValue[i] != 0) {
            mdp_safe_free((void**) &CurrValue[i]);
            i++;
        }
    }
    mdp_safe_free((void**)&CurrValue);
}
//====================================================================================================================
/*
 * LE_MultiCStr& operator=(const LE_MultiCStr &b) :
 *
 *  assignment operator
 */
LE_MultiCStr& LE_MultiCStr::operator=(const LE_MultiCStr& b)
{
    if (&b != this) {
        LineEntry::operator=(b);
        AddrVal    = b.AddrVal;
        MaxVal     = b.MaxVal;
        MinVal     = b.MinVal;

        mdp_safe_copy_string(&DefaultVal, b.DefaultVal);
        int i = 0;
        while (CurrValue[i] != 0) {
            mdp_safe_free((void**) &CurrValue[i]);
            i++;
        }
        mdp_safe_free((void**)&CurrValue);
        CurrValue = (char**) mdp_alloc_ptr_1(m_numTimesProcessed+1);
        for (int i = 0; i < m_numTimesProcessed; i++) {
            CurrValue[i]  = mdp_copy_string(b.CurrValue[i]);
        }
        LastCurrValue =  CurrValue[m_numTimesProcessed-1];
        strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
    }
    return *this;
}
//====================================================================================================================
/*
 * LineEntry* duplMyselfAsLineEntry() (virtual)
 *
 * Duplication as a base class
 */
LineEntry* LE_MultiCStr::duplMyselfAsLineEntry() const
{
    LE_MultiCStr* newLE = new LE_MultiCStr(*this);
    return (LineEntry*) newLE;
}
//====================================================================================================================
/*
 * process_lineEntry():
 *
 *   Here we interpret the token as a single string, and then
 *   assign it to the address we set up.
 */
void LE_MultiCStr::process_LineEntry(const TK_TOKEN* lineArgTok)
{
    for (int i = 0; i < NumEntryDependencies; i++) {
        BI_Dependency* idep = EntryDepList[i];
        if (idep->ResultType() ==  BIDRT_PT_ERROR) {
            if (!idep->checkDependencySatisfied()) {
                throw BI_InputError("LE_MultiCStr::process_LineEntry",
                                    "Dependency not satisfied");
            }
        }
        if (idep->ResultType() ==  BIDRT_ANTITHETICAL_ERROR) {
            if (idep->checkDependencySatisfied()) {
                throw BI_InputError("LE_MultiCStr::process_LineEntry",
                                    "Mutual Exclusive Dependency is satisfied");
            }
        }
    }

    BOOLEAN error = 0;
    /*
     * Malloc a copy of the string.
     */
    char* str = 0;
    str = tok_to_string(lineArgTok, MaxVal, MinVal, DefaultVal, &error);
    if (error) {
        throw BI_InputError("LE_MultiCStr::process_LineEntry",
                            "tok_to_string interpretation");
    }

    /*
     * Check to see whether there is an existing string in the
     * address. If there is, then free it, as we have stated that
     * it is always a malloced quantity
     */
    int len = m_numTimesProcessed+1;

    /*
     * It's possible that the entire structure is being reprocessed.
     * In that case len is reinitialized. There will be a
     * memory leak if the list is released. Since the list always
     * has a null at the end, we can make sure that all of the
     * entries that aren't needed for the current list are released here.
     */
    char** bPtrVec = 0;
    if (AddrVal) {
        bPtrVec = *AddrVal + (len-1);
        if (bPtrVec) {
            for (int i = 0; bPtrVec[i] != 0; i++) {
                mdp_safe_free((void**) &(bPtrVec[i]));
            }
        }
    }
    bPtrVec = CurrValue + (len-1);
    for (int i = 0; bPtrVec[i] != 0; i++) {
        mdp_safe_free((void**) &(bPtrVec[i]));
    }

    /*
     * Add one to the size of the list, putting a null on the end.
     */
    if (AddrVal) {
        mdp_realloc_ptr_1((void***)AddrVal, len+1, len);
    }
    mdp_realloc_ptr_1((void***)&CurrValue, len+1, len);

    /*
     * Save a copy of str in the aPtrVec vector. Note, this
     * vector will have to be freed one component at a time
     * if there is to be no memory leak. However, this
     * procedure is not part of the input routines.
     */
    char** aPtrVec = 0;
    aPtrVec = CurrValue;
    aPtrVec[len - 1] = mdp_copy_string(str);
    aPtrVec[len] = 0;

    if (AddrVal) {
        aPtrVec = *AddrVal;
        aPtrVec[len - 1] = str;
        aPtrVec[len] = 0;
    } else {
        mdp_safe_free((void**)&str);
    }

    LastCurrValue = aPtrVec[len - 1];
    m_numTimesProcessed++;
}
//====================================================================================================================
/*
 * LE_MultiCStr::print_ProcessedLine() (virtual function):
 *
 *   This routine will print out a processed line
 *
 *  The default behavior is to print the original line with a "=>"
 *  prefix to indicate that action has been taken on it.
 */
void
LE_MultiCStr::print_ProcessedLine(const TK_TOKEN* lineArgTok) const
{
    LineEntry::print_ProcessedLine(lineArgTok);
    if (strlen(PrintString) > 0) {
        printf(" ====> %s = ", PrintString);
        if (CurrValue) {
            printf(" %s\n", LastCurrValue);
        } else {
            printf(" NULL\n");
        }
    }
}
//====================================================================================================================
/*
 * print_usage() (virtual function)
 *
 */
void LE_MultiCStr::print_usage(int ilvl) const
{
    print_indent(ilvl);
    printf("%s = (Cstring)", EntryName.orig_str);
    if (MaxVal != INT_MAX || MinVal != -INT_MAX) {
        printf(" with Token limits(%d, %d)", MaxVal, MinVal);
    }
    if (DefaultVal) {
        printf(" with default %s", DefaultVal);
    }
    if (m_numTimesRequired >= 1) {
        printf(" (REQUIRED LINES = %d)", m_numTimesRequired);
    }
    printf("\n");
    if (m_numTimesProcessed) {
        for (int i = 0; i < m_numTimesProcessed; i++) {
            print_indent(ilvl + 1);
            printf(" ====> %s = ", PrintString);
            char* cc = CurrValue[i];
            if (cc) {
                printf(" %s\n", cc);
            } else {
                printf(" NULL\n");
            }
        }
    }
}
//====================================================================================================================
/*
 * Adjust base address to store the value (virtual function)
 *
 */
void LE_MultiCStr::adjustAddress(LONG_PTR addrAdjustment)
{
    if (AddrVal) {
        if (addrAdjustment != 0) {
            LONG_PTR ll = reinterpret_cast<LONG_PTR>(AddrVal);
            ll += addrAdjustment;
            AddrVal = reinterpret_cast<char***>(ll);
        }
    }
}
//====================================================================================================================
const char** LE_MultiCStr::currentTypedValue() const
{
    return ((const char**) CurrValue);
}
//====================================================================================================================
const void* LE_MultiCStr::currentValueAsVoidP() const
{
    return static_cast<const void*>(&CurrValue);
}
//====================================================================================================================
/*
 * set_default()
 *   Set the default value for this entry
 */
void LE_MultiCStr::set_default(const char* def)
{
    DefaultVal = mdp_copy_string(def);
}
//====================================================================================================================
/*
 * set_limits()
 *   This sets a limit on the maximum and minimum number of tokens
 *   allowed in a valid string.
 */
void LE_MultiCStr::set_limits(int maxV, int minV)
{
    MaxVal = maxV;
    MinVal = minV;
}
//====================================================================================================================
/*
 * set_PrintString()
 *
 */
void LE_MultiCStr::set_PrintString(const char* ps)
{
    if (ps) {
        strncpy(PrintString, ps, MAX_INPUT_STR_LN);
        PrintString[MAX_INPUT_STR_LN] = '\0';
    }
}
//====================================================================================================================

}
