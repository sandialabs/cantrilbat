/**
 * @file LE_OneBoolInt.cpp
 *  Definition for the LineEntry object LE_OneBoolInt
 */
/* $Author: hkmoffa $
* $Revision: 508 $
* $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
*/
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "LE_OneBoolInt.h"
#include "mdp_allo.h"

namespace BEInput
{

/**********************************************************************
 *
 * LE_OneBoolInt():
 *
 * Constructor:
 *
 *   This sets up the line entry special case.
 *   We make sure to call the base class constructor here to do
 *   much of the initialization.
 */
LE_OneBoolInt::LE_OneBoolInt(const char* lineName, int* aaa, int numRL,
                             const char* varName) :
    LineEntry(lineName, numRL),
    AddrVal(aaa)
{
    /*
     * Set the limits and the default values
     */
    DefaultVal = NO_DEFAULT_INT;
    CurrValue = DefaultVal;
    PrintString[0] = '\0';
    if (varName) {
        strncpy(PrintString, varName, MAX_INPUT_STR_LN);
    } else {
        strncpy(PrintString, lineName, MAX_INPUT_STR_LN);
    }
    PrintString[MAX_INPUT_STR_LN] = '\0';
}

/**********************************************************************
 *
 * LE_OneBoolInt(const LE_OneBoolInt&):
 *
 * copy Constructor:
 */
LE_OneBoolInt::LE_OneBoolInt(const LE_OneBoolInt& b) :
    LineEntry(b),
    AddrVal(b.AddrVal),
    DefaultVal(b.DefaultVal),
    CurrValue(b.CurrValue)
{
    strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
}

/**********************************************************************
 *
 * LE_OneBoolInt& operator=(const LE_OneBoolInt &b) :
 *
 *  assignment operator
 */
LE_OneBoolInt& LE_OneBoolInt::operator=(const LE_OneBoolInt& b)
{
    if (&b != this) {
        LineEntry::operator=(b);
        AddrVal    = b.AddrVal;
        DefaultVal = b.DefaultVal;
        CurrValue  = b.CurrValue;
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
LineEntry* LE_OneBoolInt::duplMyselfAsLineEntry() const
{
    LE_OneBoolInt* newLE = new LE_OneBoolInt(*this);
    return (LineEntry*) newLE;
}

/**********************************************************************
 *
 * process_lineEntry():
 *
 *   Here we interpret the token as a single boolean, and then
 *   assign it to the address we set up.
 *   Note, the address is an integer address.
 */
void LE_OneBoolInt::process_LineEntry(const TK_TOKEN* lineArgTok)
{
    if (NumEntryDependencies > 0) {
        for (int i = 0; i < NumEntryDependencies; i++) {
            BI_Dependency* idep = EntryDepList[i];
            if (idep->ResultType() == BIDRT_PT_ERROR) {
                if (!idep->checkDependencySatisfied()) {
                    throw BI_InputError("LE_OneStr::process_LineEntry",
                                        "Dependency not satisfied");
                }
            }
            if (idep->ResultType() ==  BIDRT_ANTITHETICAL_ERROR) {
                if (idep->checkDependencySatisfied()) {
                    throw BI_InputError("LE_OneStr::process_LineEntry",
                                        "Mutual Exclusive Dependency is satisfied");
                }
            }
        }
    }

    BOOLEAN error = 0;
    BOOLEAN value = 1;
    if (lineArgTok->ntokes > 0) {
        value = tok_to_boolean(lineArgTok, DefaultVal, &error);
        if (error) {
            throw BI_InputError("LE_OneBoolInt::process_LineEntry",
                                "tok_to_boolean interpretation");
        }
    }
    if (AddrVal) {
        *(AddrVal) = value;
    }
    CurrValue = value;
    m_numTimesProcessed++;

    /*
     * Now that we have the value, check to see if any dependencies based
     * on that value are broken
     */
    for (int i = 0; i < NumEntryDependencies; i++) {
        BI_Dependency* idep = EntryDepList[i];
        if (idep->ResultType() == BIDRT_RTINTMM_ERROR) {
            if (idep->checkCardIntMaxMin(CurrValue)) {
                if (!idep->checkDependencySatisfied()) {
                    throw BI_InputError("LE_OneBoolInt::process_LineEntry",
                                        "Triggered dependency check from the local input"
                                        " resulted in a denpendency failure");
                }
            }
        }
    }

}

/*
 * LE_OneBoolInt::print_ProcessedLine() (virtual function):
 *
 *   This routine will print out a processed line
 *
 *  The default behavior is to print the original line with a "=>"
 *  prefix to indicate that action has been taken on it.
 */
void
LE_OneBoolInt::print_ProcessedLine(const TK_TOKEN* lineArgTok) const
{
    LineEntry::print_ProcessedLine(lineArgTok);
    if (strlen(PrintString) > 0) {
        printf(" ====> %s = %d\n", PrintString, CurrValue);
    }
}

/*
 * print_usage() (virtual function)
 *
 */
void LE_OneBoolInt::print_usage(int ilvl) const
{
    print_indent(ilvl);
    printf("%s = (true | false)", EntryName.orig_str);
    if (DefaultVal != NO_DEFAULT_INT) {
        printf(" with default %d", DefaultVal);
    }
    if (m_numTimesRequired == 1) {
        printf(" (REQUIRED LINE)");
    }
    if (m_numTimesRequired > 1) {
        printf(" (REQUIRED LINES = %d)", m_numTimesRequired);
    }
    printf("\n");
}
//====================================================================================================================
/*
 *
 * Adjust base address to store the value (virtual function)
 *
 */
void LE_OneBoolInt::adjustAddress(LONG_PTR addrAdjustment)
{
    if (addrAdjustment != 0) {
        LONG_PTR ll = reinterpret_cast<LONG_PTR>(AddrVal);
        ll += addrAdjustment;
        AddrVal = reinterpret_cast<int*>(ll);
    }
}
//====================================================================================================================
int LE_OneBoolInt::currentTypedValue() const
{
    return CurrValue;
}
//====================================================================================================================
/*
 * Virtual Function
 */
const void* LE_OneBoolInt::currentValueAsVoidP() const
{
    return static_cast<const void*>(&CurrValue);
}
//====================================================================================================================
/*
 *
 * set_default()
 *
 */
void LE_OneBoolInt::set_default(int def)
{
    DefaultVal = def;
    /*
     * Check the storred value
     */
    bool aval = DefaultVal;
    if (AddrVal) {
        aval = *AddrVal;
    }
    if (m_numTimesProcessed == 0) {
        CurrValue = def;
        if (aval != def) {
            printf(" WARNING: LE_OneBool::set_default() for %s: Current Value, %d, doesn't agree "
                   "with default value %d. Changing value at target address.\n",
                   PrintString, aval, DefaultVal);
            if (AddrVal) {
                *AddrVal = def;
            }
        }
    }

}

/**********************************************************************
 *
 * set_PrintString()
 *
 */
void LE_OneBoolInt::set_PrintString(const char* ps)
{
    if (ps) {
        strncpy(PrintString, ps, MAX_INPUT_STR_LN);
        PrintString[MAX_INPUT_STR_LN] = '\0';
    }
}
/**********************************************************************/
}
