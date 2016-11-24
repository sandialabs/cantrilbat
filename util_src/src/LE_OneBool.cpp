/**
 * @file LE_OneBool.cpp
 * Definitions for the LineEntry object, LE_OneBool.
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "LE_OneBool.h"
#include "mdp_allo.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace BEInput
{
//==================================================================================================================================
LE_OneBool::LE_OneBool(const char* lineName, bool* aaa, int numRL, const char* varName) :
    LineEntry(lineName, numRL),
    AddrVal(aaa)
{
    DefaultVal = NO_DEFAULT_INT;
    CurrValue = false;
    PrintString[0] = '\0';
    if (varName) {
        strncpy(PrintString, varName, MAX_INPUT_STR_LN);
        PrintString[MAX_INPUT_STR_LN] = '\0';
    } else {
        strncpy(PrintString, lineName, MAX_INPUT_STR_LN);
        PrintString[MAX_INPUT_STR_LN] = '\0';
    }
}
//==================================================================================================================================
LE_OneBool::LE_OneBool(const LE_OneBool& b) :
    LineEntry(b),
    AddrVal(b.AddrVal),
    DefaultVal(b.DefaultVal),
    CurrValue(b.CurrValue)
{
    strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
}
//==================================================================================================================================
LE_OneBool& LE_OneBool::operator=(const LE_OneBool& b)
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
//==================================================================================================================================
LineEntry* LE_OneBool::duplMyselfAsLineEntry() const
{
    LE_OneBool* newLE = new LE_OneBool(*this);
    return (LineEntry*) newLE;
}
//==================================================================================================================================
void LE_OneBool::process_LineEntry(const TK_TOKEN* lineArgTok)
{
    if (NumEntryDependencies > 0) {
        for (int i = 0; i < NumEntryDependencies; i++) {
            BI_Dependency* idep = EntryDepList[i];
            if (idep->ResultType() == BIDRT_PT_ERROR) {
                if (!idep->checkDependencySatisfied()) {
                    throw BI_InputError("LE_OneBool::process_LineEntry", "Dependency not satisfied");
                }
            }
            if (idep->ResultType() ==  BIDRT_ANTITHETICAL_ERROR) {
                if (idep->checkDependencySatisfied()) {
                    throw BI_InputError("LE_OneBool::process_LineEntry", "Mutual Exclusive Dependency is satisfied");
                }
            }
        }
    }

    BOOLEAN error = 0;
    BOOLEAN value = 1;
    if (lineArgTok->ntokes > 0) {
        value = tok_to_boolean(lineArgTok, DefaultVal, &error);
        if (error) {
            throw BI_InputError("LE_OneBool::process_LineEntry", "tok_to_boolean interpretation");
        }
    }
    if (value == 0) {
        if (AddrVal) {
            *(AddrVal) = false;
        }
        CurrValue = false;
    } else {
        if (AddrVal) {
            *(AddrVal) = true;
        }
        CurrValue = true;
    }
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
                    throw
                    BI_InputError("LE_OneBool::process_LineEntry",
                                  "Triggered dependency check from the local input resulted in a denpendency failure");
                }
            }
        }
    }
}
//==================================================================================================================================
/*
 *  The default behavior is to print the original line with a "=>"
 *  prefix to indicate that action has been taken on it.
 */
void
LE_OneBool::print_ProcessedLine(const TK_TOKEN* lineArgTok) const
{
    LineEntry::print_ProcessedLine(lineArgTok);
    if (strlen(PrintString) > 0) {
        printf(" ====> %s = %d\n", PrintString, CurrValue);
    }
}
//==================================================================================================================================
void LE_OneBool::print_usage(int ilvl) const
{
    print_indent(ilvl);
    printf("%s = (true | false)", EntryName.orig_str);
    if (DefaultVal != NO_DEFAULT_INT) {
        if (DefaultVal) {
            printf(" with default true");
        } else {
            printf(" with default false");
        }
    }
    if (m_numTimesRequired == 1) {
        printf(" (REQUIRED LINE)");
    }
    if (m_numTimesRequired > 1) {
        printf(" (REQUIRED LINES = %d)", m_numTimesRequired);
    }
    if (m_numTimesRequired == 0) {
        printf(" (OPTIONAL LINE)");
    }
    printf("\n");
}
//==================================================================================================================================
void LE_OneBool::adjustAddress(LONG_PTR addrAdjustment)
{
    if (AddrVal) {
        if (addrAdjustment != 0) {
            LONG_PTR ll = reinterpret_cast<LONG_PTR>(AddrVal);
            ll += addrAdjustment;
            AddrVal = reinterpret_cast<bool*>(ll);
        }
    }
}
//==================================================================================================================================
bool LE_OneBool::currentTypedValue() const
{
    return CurrValue;
}
//==================================================================================================================================
const void* LE_OneBool::currentValueAsVoidP() const
{
    return static_cast<const void*>(&CurrValue);
}
//==================================================================================================================================
void LE_OneBool::set_default(bool def)
{
    DefaultVal = def;
    bool aval = DefaultVal;
    if (AddrVal) {
        aval = *AddrVal;
    }
    if (m_numTimesProcessed == 0) {
        CurrValue = def;
        if (aval != def) {
            printf(" WARNING: LE_OneBool::set_default() for %s: Current Value, %d, doesn't agree "
                   "with default value %d. Changing value at target address.\n", PrintString, aval, DefaultVal);
            if (AddrVal) {
                *AddrVal = def;
            }
        }
    }
}
//==================================================================================================================================
void LE_OneBool::set_PrintString(const char* ps)
{
    if (ps) {
        strncpy(PrintString, ps, MAX_INPUT_STR_LN);
        PrintString[MAX_INPUT_STR_LN] = '\0';
    }
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------

