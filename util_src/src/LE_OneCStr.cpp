/**
 * @file LE_OneCStr.cpp
 *  Header for the LineEntry of a single C string
 *  (see \ref blockentryModule and
 *    class \link BEInput::LE_OneCStr LE_OneCStr\endlink).
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

#include "LE_OneCStr.h"
#include "mdp_allo.h"

#include <climits>
#include <cstdlib>

using namespace mdpUtil;

namespace BEInput
{

/*
 *
 * LE_OneCStr Constructor:
 *
 *   This sets up the line entry special case.
 *   We make sure to call the base class constructor here to do
 *   much of the initialization.
 *
 *   The maximum number of tokens allowed in the string is equal
 *  to the limit specified in the tok_input_util.h file.
 *   The minimum number of tokens is equal to one. In other words we
 *  require something in the string as a default.
 */
LE_OneCStr::LE_OneCStr(const char* lineName, char** aaa, int maxval,
                       int minval, int numRL, const char* varName) :
    LineEntry(lineName, numRL),
    AddrVal(aaa),
    MaxVal(maxval),
    MinVal(minval),
    DefaultVal(0),
    CurrValue(0)
{
    /*
     * Set the limits and the default values
     */
    CurrValue = DefaultVal;
    PrintString[0] = '\0';
    if (varName) {
        strncpy(PrintString, varName, MAX_INPUT_STR_LN);
    } else {
        strncpy(PrintString, lineName, MAX_INPUT_STR_LN);
    }
    PrintString[MAX_INPUT_STR_LN] = '\0';
    if (numRL > 1) {
        printf("LE_OneCStr ERROR: num required times is 0 or 1\n");
        exit(-1);
    }
}

/**********************************************************************
 *
 * LE_OneCStr(const LE_OneCStr&):
 *
 * copy Constructor:
 */
LE_OneCStr::LE_OneCStr(const LE_OneCStr& b) :
    LineEntry(b),
    AddrVal(b.AddrVal),
    MaxVal(b.MaxVal),
    MinVal(b.MinVal)
{
    DefaultVal = mdp_copy_string(b.DefaultVal);
    CurrValue  = mdp_copy_string(b.CurrValue);
    strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
}

/**********************************************************************
 *
 * ~LE_OneCStr():
 *
 *  Destructor eliminates the space possibly malloced by this structure
 */
LE_OneCStr::~LE_OneCStr()
{
    mdp_safe_free((void**)&DefaultVal);
    mdp_safe_free((void**)&CurrValue);
}

/**********************************************************************
 *
 * LE_OneCStr& operator=(const LE_OneCStr &b) :
 *
 *  assignment operator
 */
LE_OneCStr& LE_OneCStr::operator=(const LE_OneCStr& b)
{
    if (&b != this) {
        LineEntry::operator=(b);
        AddrVal    = b.AddrVal;
        MaxVal     = b.MaxVal;
        MinVal     = b.MinVal;
        mdp_safe_copy_string(&DefaultVal, b.DefaultVal);
        mdp_safe_copy_string(&CurrValue, b.CurrValue);
        strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
    }
    return *this;
}

/*********************************************************************
 *
 * LineEntry* duplMyselfAsLineEntry() (virtual)
 *
 * Duplication as a base class
 */
LineEntry* LE_OneCStr::duplMyselfAsLineEntry() const
{
    LE_OneCStr* newLE = new LE_OneCStr(*this);
    return (LineEntry*) newLE;
}

/**********************************************************************
 *
 * process_lineEntry():
 *
 *   Here we interpret the token as a single string, and then
 *   assign it to the address we set up.
 */
void LE_OneCStr::process_LineEntry(const TK_TOKEN* lineArgTok)
{
    for (int i = 0; i < NumEntryDependencies; i++) {
        BI_Dependency* idep = EntryDepList[i];
        if (idep->ResultType() == BIDRT_PT_ERROR) {
            if (!idep->checkDependencySatisfied()) {
                throw BI_InputError("LE_OneCStr::process_LineEntry",
                                    "Dependency not satisfied");
            }
        }
        if (idep->ResultType() ==  BIDRT_ANTITHETICAL_ERROR) {
            if (idep->checkDependencySatisfied()) {
                throw BI_InputError("LE_OneCStr::process_LineEntry",
                                    "Mutual Exclusive Dependency is satisfied");
            }
        }
    }

    BOOLEAN error = 0;
    /*
     * Malloc a copy of the string.
     */
    char* str = tok_to_string(lineArgTok, MaxVal, MinVal,
                              DefaultVal, &error);
    if (error) {
        throw BI_InputError("LE_OneCStr::process_LineEntry",
                            "tok_to_string interpretation");
    }

    /*
     * Check to see whether there is an existing string in the
     * address. If there is, then free it, as we have stated that
     * it is always a malloced quantity
     */
    if (*AddrVal) {
        mdp_safe_free((void**) AddrVal);
    }
    *AddrVal = str;
    CurrValue = mdp_copy_string(str);
    m_numTimesProcessed++;
}

/*
 * LE_OneCStr::print_ProcessedLine() (virtual function):
 *
 *   This routine will print out a processed line
 *
 *  The default behavior is to print the original line with a "=>"
 *  prefix to indicate that action has been taken on it.
 */
void
LE_OneCStr::print_ProcessedLine(const TK_TOKEN* lineArgTok) const
{
    LineEntry::print_ProcessedLine(lineArgTok);
    if (strlen(PrintString) > 0) {
        printf(" ====> %s = ", PrintString);
        if (CurrValue) {
            printf(" %s\n", CurrValue);
        } else {
            printf(" NULL\n");
        }
    }
}

/*
 * print_usage() (virtual function)
 *
 */
void LE_OneCStr::print_usage(int ilvl) const
{
    print_indent(ilvl);
    printf("%s = (Cstring)", EntryName.orig_str);
    if (MaxVal != INT_MAX || MinVal != -INT_MAX) {
        printf(" with Token limits(%d, %d)", MaxVal, MinVal);
    }
    if (DefaultVal) {
        printf(" with default %s", DefaultVal);
    }
    if (m_numTimesRequired == 1) {
        printf(" (REQUIRED LINE)");
    }
    if (m_numTimesRequired > 1) {
        printf(" (REQUIRED LINES = %d)", m_numTimesRequired);
    }
    printf("\n");
}

/**********************************************************************
 *
 * Adjust base address to store the value (virtual function)
 *
 */
void LE_OneCStr::adjustAddress(LONG_PTR addrAdjustment)
{
    if (addrAdjustment != 0) {
        LONG_PTR ll = reinterpret_cast<LONG_PTR>(AddrVal);
        ll += addrAdjustment;
        AddrVal = reinterpret_cast<char**>(ll);
    }
}

const char* LE_OneCStr::currentTypedValue() const
{
    return (CurrValue);
}

const void* LE_OneCStr::currentValueAsVoidP() const
{
    return static_cast<const void*>(&CurrValue);
}
//====================================================================================================================
/*
 *
 * set_default()
 *   Set the default value for this entry
 */
void LE_OneCStr::set_default(const char* def)
{
    mdp_safe_copy_string(&DefaultVal, def);
    if (m_numTimesProcessed == 0) {
        mdp_safe_copy_string(&CurrValue, def);
        if (AddrVal) {
            mdp_safe_copy_string(AddrVal, def);
        }
    }
}
//====================================================================================================================
/*
 *
 * set_limits()
 *   This sets a limit on the maximum and minimum number of tokens
 *   allowed in a valid string.
 */
void LE_OneCStr::set_limits(int maxV, int minV)
{
    MaxVal = maxV;
    MinVal = minV;
}

/**********************************************************************
 *
 * set_PrintString()
 *
 */
void LE_OneCStr::set_PrintString(const char* ps)
{
    if (ps) {
        strncpy(PrintString, ps, MAX_INPUT_STR_LN);
        PrintString[MAX_INPUT_STR_LN] = '\0';
    }
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

}
