/**
 * @file LE_OneDblUnits.cpp
 *  Definitions for the LineEntry of a single double with unit conversion
 *  (see \ref blockentryModule and
 *   class \link BEInput::LE_OneDblUnits LE_OneDblUnits\endlink).
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

#include "LE_OneDblUnits.h"
#include "mdp_allo.h"
#include "BE_UnitConversion.h"
#include "mdp_stringUtils.h"

#include <cfloat>

using std::string;
using namespace TKInput;

namespace BEInput
{

/*
 *
 * LE_OneDblUnits Constructor:
 *
 *   This sets up the line entry special case.
 *   We make sure to call the base class constructor here to do
 *   much of the initialization.
 */
LE_OneDblUnits::LE_OneDblUnits(const char* lineName, double* aaa,
                               int numRL,
                               const char* varName,
                               BE_UnitConversion* uc) :
    LineEntry(lineName, numRL),
    AddrVal(aaa),
    m_uc(uc)
{
    /*
     * Set the limits and the default values
     */
    MaxVal = DBL_MAX;
    MinVal = -DBL_MAX;
    DefaultVal = NO_DEFAULT_DOUBLE;
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
 * LE_OneDblUnits(const LE_OneDblUnits&):
 *
 * copy Constructor:
 * Notes -> shallow copy for units converter.
 */
LE_OneDblUnits::LE_OneDblUnits(const LE_OneDblUnits& b) :
    LineEntry(b),
    AddrVal(b.AddrVal),
    MaxVal(b.MaxVal),
    MinVal(b.MinVal),
    DefaultVal(b.DefaultVal),
    CurrValue(b.CurrValue),
    m_uc(0)
{
    m_uc = (b.m_uc)->duplMyselfAsUnitConversion();
    strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
}

/**********************************************************************
 *
 * LE_OneDblUnits& operator=(const LE_OneDblUnits &b) :
 *
 *  assignment operator
 */
LE_OneDblUnits& LE_OneDblUnits::operator=(const LE_OneDblUnits& b)
{
    if (&b != this) {
        LineEntry::operator=(b);
        AddrVal    = b.AddrVal;
        MaxVal     = b.MaxVal;
        MinVal     = b.MinVal;
        DefaultVal = b.DefaultVal;
        CurrValue  = b.CurrValue;
        m_uc       = (b.m_uc)->duplMyselfAsUnitConversion();

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
LineEntry* LE_OneDblUnits::duplMyselfAsLineEntry() const
{
    LE_OneDblUnits* newLE = new LE_OneDblUnits(*this);
    return (LineEntry*) newLE;
}

/**********************************************************************
 *
 * LE_OneDblUnits Destructor:
 *
 *  delete units converter.
 */
LE_OneDblUnits::~LE_OneDblUnits()
{
    if (m_uc) {
        delete m_uc;
        m_uc = 0;
    }
}

/**********************************************************************
 *
 * process_lineEntry():
 *
 *   Here we interpret the token as a single double, and then
 *   assign it to the address we set up.
 */
void LE_OneDblUnits::
process_LineEntry(const TK_TOKEN* lineArgTok)
{
    // do the default processing first
    LineEntry::process_LineEntry(lineArgTok);

    double value = DBL_MAX;
    BOOLEAN error = 0;
    if (lineArgTok) {
        if (lineArgTok->ntokes == 1 || !m_uc) {
            value = tok_to_double(lineArgTok, MaxVal, MinVal, DefaultVal,
                                  &error);
        } else {
            char* vString = lineArgTok->tok_ptrV[0];
            // If there is a max and min, we can't compare until we do the units conversion.
            //	value = str_to_double(vString, MaxVal, MinVal, DefaultVal, &error);
            value = str_to_double(vString, DBL_MAX, -DBL_MAX, DefaultVal, &error);
            char* uString = lineArgTok->tok_ptrV[1];
            if (lineArgTok->ntokes > 2) {
                throw BI_InputError("LE_OneDblUnits::process_LineEntry",
                                    "spaces aren't allowed - interpretation error");
            }
            string uuString(uString);
            double fchange = m_uc->toSI(uuString);
            value *= fchange;
            if (value > MaxVal || value < MinVal) {
                throw BI_InputError("LE_OneDblUnits OutofBounds:",
                                    mdpUtil::fp2str(value) + " is not within max = " + mdpUtil::fp2str(MaxVal) + ", min = " + mdpUtil::fp2str(MinVal));
            }
        }
    } else {
        error = true;
    }
    if (error) {
        throw BI_InputError("LE_OneDblUnits::process_LineEntry",
                            "tok_to_double interpretation");
    }
    *AddrVal = value;
    CurrValue = value;
}

/*
 *   This routine will print out a processed line
 *
 *  The default behavior is to print the original line with a "=>"
 *  prefix to indicate that action has been taken on it.
 */
void
LE_OneDblUnits::print_ProcessedLine(const TK_TOKEN* lineArgTok) const
{
    LineEntry::print_ProcessedLine(lineArgTok);
    if (strlen(PrintString) > 0) {
        printf(" ====> %s = %g\n", PrintString, CurrValue);
    }
}

/**********************************************************************
 *
 * print_usage() (virtual function)
 *
 */
void LE_OneDblUnits::print_usage(int ilvl) const
{
    print_indent(ilvl);
    printf("%s = (double)", EntryName.orig_str);
    if (m_uc) {
        printf(" [Units]");
    }
    if (MaxVal != DBL_MAX || MinVal != -DBL_MAX) {
        printf(" with limits (%g, %g)", MinVal, MaxVal);
    }
    if (DefaultVal != NO_DEFAULT_DOUBLE) {
        printf(" with default %g", DefaultVal);
    }
    if (m_numTimesRequired == 1) {
        printf(" (REQUIRED LINE)");
    }
    if (m_numTimesRequired > 1) {
        printf(" (REQUIRED LINES = %d)", m_numTimesRequired);
    }
    printf("\n");
    if (m_uc) {
        printf("\tUnits = ");
        string pu = m_uc->returnUsage();
        printf("%s", pu.c_str());
        printf("\n");
    }
}

/**********************************************************************
 *
 * Adjust base address to store the value (virtual function)
 *
 */
void LE_OneDblUnits::adjustAddress(LONG_PTR addrAdjustment)
{
    if (addrAdjustment != 0) {
        LONG_PTR ll = reinterpret_cast<LONG_PTR>(AddrVal);
        ll += addrAdjustment;
        AddrVal = reinterpret_cast<double*>(ll);
    }
}

double LE_OneDblUnits::currentTypedValue() const
{
    return CurrValue;
}

/*
 * currentValueAsVoidP() (virtual)
 */
const void* LE_OneDblUnits::currentValueAsVoidP() const
{
    return static_cast<const void*>(&CurrValue);
}
//====================================================================================================================
/*
 * set_default()
 *
 */
void LE_OneDblUnits::set_default(double def)
{
    DefaultVal = def;
    /*
     * Check the storred value
     */
    double val = DefaultVal;
    if (AddrVal) {
        val = *AddrVal;
    }
    if (m_numTimesProcessed == 0) {
        /*
         * Install the default value into the current value if it
         * makes sense to do so.
         */
        CurrValue = DefaultVal;
        if (val != DefaultVal) {
            printf(" WARNING: LE_OneDblUnits::set_default() for %s: Current Value, %g, doesn't agree "
                   "with default value %g. Changing value at target address.\n",
                   PrintString, val, DefaultVal);
            *AddrVal = DefaultVal;
        }
    }
}
//====================================================================================================================
/*
 *
 * set_limits()
 *
 */
void LE_OneDblUnits::set_limits(double maxV, double minV)
{
    MaxVal = maxV;
    MinVal = minV;
}

/*
 *
 *set_PrintString()
 *
 */
void LE_OneDblUnits::set_PrintString(const char* ps)
{
    if (ps) {
        strncpy(PrintString, ps, MAX_INPUT_STR_LN);
        PrintString[MAX_INPUT_STR_LN] = '\0';
    }
}

}
