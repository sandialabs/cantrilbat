/**
 * @file LE_PickList.cpp
 *  Definitions for the LineEntry of an item chosen from a list
 *  (see \ref blockentryModule and class
 *  \link BEInput::LE_PickList LE_PickList\endlink).
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

#include "LE_PickList.h"
#include "mdp_allo.h"

#include "string.h"

using namespace TKInput;
using namespace mdpUtil;

namespace BEInput
{
//====================================================================================================================
/*
 * LE_PickList Constructor:
 *
 *   This sets up the line entry special case.
 *   We make sure to call the base class constructor here to do
 *   much of the initialization.
 */
LE_PickList::LE_PickList(const char* lineName, int* aaa,
                         const char** charList, int listLength, int numRL,
                         const char* varName) :
    LineEntry(lineName, numRL),
    AddrVal(aaa),
    CharList(0),
    ListLength(listLength)
{
    const char* item;
    /*
     * Set the limits and the default values
     */
    MaxVal = ListLength - 1;
    MinVal = 0;
    DefaultVal = NO_DEFAULT_INT;
    CurrValue = DefaultVal;
    if (charList) {
        CharList =
            mdp_alloc_VecFixedStrings(ListLength, MAX_INPUT_STR_LN+1);
        for (int i = 0; i < ListLength; i++) {
            item = charList[i];
            if (!item) throw BI_InputError("LE_PickList::LE_PickList",
                                               "char list item is null");
            strncpy(CharList[i], item, MAX_INPUT_STR_LN);
            CharList[i][MAX_INPUT_STR_LN] = '\0';
        }
    }
    PrintString[0] = '\0';
    if (varName) {
        strncpy(PrintString, varName, MAX_INPUT_STR_LN);
        PrintString[MAX_INPUT_STR_LN] = '\0';
    } else {
        strncpy(PrintString, lineName, MAX_INPUT_STR_LN);
        PrintString[MAX_INPUT_STR_LN] = '\0';
    }
}
//====================================================================================================================
/*
 * LE_PickList(const LE_PickList&):
 *
 * copy Constructor:
 */
LE_PickList::LE_PickList(const LE_PickList& b) :
    LineEntry(b),
    AddrVal(b.AddrVal),
    MaxVal(b.MaxVal),
    MinVal(b.MinVal),
    DefaultVal(b.DefaultVal),
    CharList(0),
    ListLength(b.ListLength),
    CurrValue(b.CurrValue)
{
    if (ListLength > 0) {
        char* item;
        CharList = mdp_alloc_VecFixedStrings(ListLength, MAX_INPUT_STR_LN+1);
        for (int i = 0; i < ListLength; i++) {
            item = b.CharList[i];
            if (!item) throw BI_InputError("LE_PickList::LE_PickList",
                                               "char list item is null");
            strncpy(CharList[i], item, MAX_INPUT_STR_LN);
            CharList[i][MAX_INPUT_STR_LN] = '\0';
        }
    }
    strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
}
//====================================================================================================================
/*
 *
 * LE_PickList& operator=(const LE_PickList &b) :
 *
 *  assignment operator
 */
LE_PickList& LE_PickList::operator=(const LE_PickList& b)
{
    if (&b != this) {
        LineEntry::operator=(b);
        AddrVal    = b.AddrVal;
        MaxVal     = b.MaxVal;
        MinVal     = b.MinVal;
        ListLength = b.ListLength;
        DefaultVal = b.DefaultVal;
        CurrValue  = b.CurrValue;
        strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
        mdp_safe_free((void**) &CharList);
        if (ListLength > 0) {
            char* item;
            CharList =
                mdp_alloc_VecFixedStrings(ListLength, MAX_INPUT_STR_LN+1);
            for (int i = 0; i < ListLength; i++) {
                item = b.CharList[i];
                if (!item) throw BI_InputError("LE_PickList::LE_PickList",
                                                   "char list item is null");
                strncpy(CharList[i], item, MAX_INPUT_STR_LN);
                CharList[i][MAX_INPUT_STR_LN] = '\0';
            }
        }
    }
    return *this;
}
//====================================================================================================================
/*
 *
 * LineEntry* duplMyselfAsLineEntry() (virtual)
 *
 * Duplication as a base class
 */
LineEntry* LE_PickList::duplMyselfAsLineEntry() const
{
    LE_PickList* newLE = new LE_PickList(*this);
    return (LineEntry*) newLE;
}
//====================================================================================================================
/*
 *
 * LE_PickList destructor:
 *
 * We malloced memory here, so we must explicitly call free.
 */
LE_PickList::~LE_PickList()
{
    mdp_safe_free((void**) &CharList);
}
//====================================================================================================================
/*
 *
 * process_lineEntry():
 *
 *   Here we interpret the token as a single string, and then
 *   compare it to the list of strings in the list that was set up.
 */
void LE_PickList::process_LineEntry(const TK_TOKEN* lineArgTok)
{
    std::string ename = EntryName.orig_str;
    for (int i = 0; i < NumEntryDependencies; i++) {
        BI_Dependency* idep = EntryDepList[i];
        std::string dname = idep->TargetBaseEntry()->keyname();
        if (idep->ResultType() ==  BIDRT_USETOPROCESS) {
            int value;
            if (idep->checkDepOneInt(value)) {
                DefaultVal = value;
            } else {
                throw BI_InputError("LE_PickList::process_LineEntry  on \"" + ename + "\"",
                                    "Dependency not satisfied  wrt keyname \"" + dname + "\"");
            }
        }
        if (idep->ResultType() ==  BIDRT_PT_ERROR) {
            if (!idep->checkDependencySatisfied()) {
                throw BI_InputError("LE_PickList::process_LineEntry on \"" + ename + "\"",
                                    "Dependency not satisfied wrt keyname \"" + dname + "\"");
            }
        }

        if (idep->ResultType() ==  BIDRT_ANTITHETICAL_ERROR) {
            if (idep->checkDependencySatisfied()) {
                throw BI_InputError("LE_PickList::process_LineEntry on \"" + ename + "\"",
                                    "Mutual Exclusive Dependency is satisfied wrt keyname \"" + dname + "\"");
            }
        }
    }

    BOOLEAN error = 0;
    char* cvalue = 0;
    if (DefaultVal != NO_DEFAULT_INT) {
        cvalue = tok_to_string(lineArgTok, MAXTOKENS, 0, CharList[DefaultVal], &error);
    } else {
        cvalue = tok_to_string(lineArgTok, MAXTOKENS, 0, NO_DEFAULT_STR,  &error);
    }
    int ifind = -1;
    for (int i = 0; i < ListLength; i++) {
        if (strstrmatch(cvalue, CharList[i])) {
            ifind = i;
            break;
        }
    }
    if (ifind >= 0) {
        CurrValue = ifind;
        *AddrVal = CurrValue;
    } else {
        std::string tmp = " with list entries:\n";
        for (int i = 0; i < ListLength; i++) {
            tmp += std::string("\t\t") + CharList[i] + std::string("\n");
        }
        tmp += std::string("\t\t");
        throw BI_UnknownListEntry("LE_PickList::process_LineEntry" + tmp, cvalue);
    }
    m_numTimesProcessed++;

    /*
     * Now that we have the value, check to see if any dependencies based
     * on that value are broken
     */
    for (int i = 0; i < NumEntryDependencies; i++) {
        BI_Dependency* idep = EntryDepList[i];
        std::string dname = idep->TargetBaseEntry()->keyname();
        if (idep->ResultType() == BIDRT_RTINTMM_ERROR) {
            if (idep->checkCardIntMaxMin(CurrValue)) {
                if (!idep->checkDependencySatisfied()) {
                    throw BI_InputError("LE_PickList::process_LineEntry on \"" + ename + "\"",
                                        "Triggered dependency check from the local input"
                                        " resulted in a dependency failure wrt keyname \"" + dname + "\"");
                }
            }
        }
    }

    mdp_safe_free((void**) &cvalue);
}
//====================================================================================================================
/*
 * LE_PickList::print_ProcessedLine() (virtual function):
 *
 *   This routine will print out a processed line. It is only called
 *   if the static member of the class LineEntry, PrintProcessedLine,
 *   is true.
 *
 *  The default behavior is to print the original line with a "=>"
 *  prefix to indicate that action has been taken on it.
 */
void
LE_PickList::print_ProcessedLine(const TK_TOKEN* lineArgTok) const
{
    LineEntry::print_ProcessedLine(lineArgTok);
    if (strlen(PrintString) > 0) {
        printf(" ====> %s = %s, list item %d\n", PrintString,
               CharList[CurrValue], CurrValue);
    }
}
//====================================================================================================================
/*
 *
 * print_usage() (virtual function)
 *
 */
void LE_PickList::print_usage(int ilvl) const
{
    print_indent(ilvl);
    printf("%s = (LIST)", EntryName.orig_str);
    if (DefaultVal != NO_DEFAULT_INT) {
        printf(" with default %s : %d",  CharList[DefaultVal], DefaultVal);
    }
    if (m_numTimesRequired == 1) {
        printf(" (REQUIRED LINE)");
    }
    if (m_numTimesRequired > 1) {
        printf(" (REQUIRED LINES = %d)", m_numTimesRequired);
    }
    printf(" where:\n");
    print_indent(ilvl+1);
    printf("LIST = \n");
    for (int i = 0; i < ListLength; i++) {
        print_indent(ilvl+2);
        printf("%s : %d\n", CharList[i], i);
    }
}
//====================================================================================================================
/*
 * Adjust base address to store the value (virtual function)
 *
 */
void LE_PickList::adjustAddress(LONG_PTR addrAdjustment)
{
    if (addrAdjustment != 0) {
        LONG_PTR ll = reinterpret_cast<LONG_PTR>(AddrVal);
        ll += addrAdjustment;
        AddrVal = reinterpret_cast<int*>(ll);
    }
}
//====================================================================================================================
int LE_PickList::currentTypedValue() const
{
    return (CurrValue);
}
//====================================================================================================================
const void* LE_PickList::currentValueAsVoidP() const
{
    return static_cast<const void*>(&CurrValue);
}
//====================================================================================================================
/*
 * set_default()
 *
 */
void LE_PickList::set_default(int def)
{
    DefaultVal = def;
    /*
     * Check the storred value
     */
    int aval = DefaultVal;
    if (AddrVal) {
        aval = *AddrVal;
    }
    if (m_numTimesProcessed == 0) {
        CurrValue = def;
        if (aval != def) {
            printf(" WARNING: LE_PickList::set_default() for %s: Current Value, %d, doesn't agree "
                   "with default value %d. Changing value at target address.\n",
                   PrintString, aval, DefaultVal);
            if (AddrVal) {
                *AddrVal = def;
            }
        }
    }
}
//====================================================================================================================
/*
 *
 * set_default()
 *
 */
void LE_PickList::set_default(const char* defString)
{
    int ifind = -1;
    for (int i = 0; i < ListLength; i++) {
        if (strstrmatch(defString, CharList[i])) {
            ifind = i;
            break;
        }
    }
    if (ifind >= 0) {
        DefaultVal = ifind;
        if (m_numTimesProcessed == 0) {
            CurrValue = ifind;
            int aval = DefaultVal;
            if (AddrVal) {
                aval = *AddrVal;
            }
            if (aval != DefaultVal) {
                printf(" WARNING: LE_PickList::set_default() for %s: Current target addr Value, %d, doesn't agree "
                       "with default value %d. Changing value at target address.\n",
                       PrintString, aval, DefaultVal);
                *AddrVal = DefaultVal;
            }
        }
    } else {
        throw BI_UnknownListEntry("LE_PickList::set_default unknown", defString);
    }
}
//====================================================================================================================
/*
 *
 * set_limits()
 *
 */
void LE_PickList::set_limits(int maxV, int minV)
{
    MaxVal = maxV;
    MinVal = minV;
}
//====================================================================================================================
/*
 *
 *set_PrintString()
 *
 */
void LE_PickList::set_PrintString(const char* ps)
{
    if (ps) {
        strncpy(PrintString, ps, MAX_INPUT_STR_LN);
        PrintString[MAX_INPUT_STR_LN] = '\0';
    }
}
//====================================================================================================================
/*
 *
 * DepCanService():
 *
 *  This call checks to whether requests can be fulfilled by this
 *  Line Entry object. Basically, this object can return whether it
 *  has been called enough times, and it can return the integer
 *  value that was last processed.
 */
bool LE_PickList::DepCanService(BIDSR_TYPE BIDSR_value) const
{
    bool retn = BaseEntry::DepCanService(BIDSR_value);
    if (!retn) {
        if (BIDSR_value == BIDSR_ONEINT) {
            return true;
        }
    }
    return retn;
}
//====================================================================================================================
/*
 *
 * ansDepCheckOneInt(): (virtual)
 *
 * This call returns a bool indicating whether it has been called before.
 * And if it has been, it returns the last value picked/processed.
 */
bool LE_PickList::ansDepCheckOneInt(int& returnInt) const
{
    bool retn = ansDepCheck();
    if (!retn) {
        returnInt = DefaultVal;
        return retn;
    }
    returnInt = CurrValue;
    return retn;
}
//====================================================================================================================
}

