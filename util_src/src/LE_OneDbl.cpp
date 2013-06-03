/**
 * @file LE_OneDbl.cpp
 *  Definitions for the LineEntry of a single double
 *  (see \ref blockentryModule and 
 *   class \link BEInput::LE_OneDbl LE_OneDbl\endlink).
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

#include "LE_OneDbl.h"
#include "mdp_allo.h"

#include "stdio.h"
#include "string.h"
#include "limits.h"
#include "float.h"

namespace BEInput {

  /*
   *   This sets up the line entry special case for the entry of a single
   *   double input.
   *  
   *   @param lineName c-string containing the keyline id.
   *   @param aaa      Address of the location where the input double 
   *                   will be stored.
   *   @param numRL    Number of required times this line should be in the
   *                   input deck. 
   *   @param varName  c-string used for id purposes for IO.
   *
   *  Example:
   *
   *   The example below sets up a required keyline entry, putting the value in
   *   the global address globInput.viscosity.
   *
   *  Setup:
   *      struct globInput {
   *         double viscosity;
   *      };
   *
   *      BlockEntry *besmd;
   *      int reqd = 1;
   *      LE_OneDbl *d2 = new LE_OneDbl("Viscosity", 
   *                                    &globInput.viscosity, 
   *                                    reqd, "viscosity");
   *      d2->set_default(34.3.);
   *      d2->set_limits(1.0E20, 0.0);
   *      besmd->addLineEntry(d2);
   *
   *  Entry:
   *       Viscosity = 34.3
   *
   *  Running:
   *     besmd->read_block(fp *inputfile);
   *
   *  Dependencies that may be set on this card that makes sense.
   *          BIDRT_PT_ERROR
   *          BIDRT_ANTITHETICAL_ERROR
   *          BIDRT_ZERONUMTIMESREQUIRED
   *          BIDRT_ONENUMTR
   */
  LE_OneDbl::LE_OneDbl(const char *lineName, double *aaa, int numRL,
		       const char *varName) :
    LineEntry(lineName, numRL),
    AddrVal(aaa)
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
      PrintString[MAX_INPUT_STR_LN] = '\0';
    } else {
      strncpy(PrintString, lineName, MAX_INPUT_STR_LN);
      PrintString[MAX_INPUT_STR_LN] = '\0';
    }
  }

  /*
   * LE_OneDbl(const LE_OneBoolInt&):
   *
   * copy Constructor:
   */
  LE_OneDbl::LE_OneDbl(const LE_OneDbl& b) : 
    LineEntry(b),
    AddrVal(b.AddrVal),
    MaxVal(b.MaxVal),
    MinVal(b.MinVal),
    DefaultVal(b.DefaultVal),
    CurrValue(b.CurrValue)
  {
    strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
  }

  /*
   *
   * process_lineEntry():
   *
   *   Here we interpret the token as a single double, and then
   *   assign it to the address we set up.
   */
  void LE_OneDbl::process_LineEntry(const TK_TOKEN *lineArgTok)
  {
 
    /*
     * Handle common dependencies and incr NumTimesProcessed.
     */
    LineEntry::process_LineEntry(lineArgTok);
 
    BOOLEAN error = 0;
    double value = tok_to_double(lineArgTok, MaxVal, MinVal, DefaultVal, 
				 &error);
    if (error) {
      throw BI_InputError("LE_OneDbl::process_LineEntry",
			  "tok_to_double interpretation");
    }
    if (AddrVal) {
      *AddrVal = value;
    }
    CurrValue = value;
  }

  /*
   *
   * LE_OneDbl& operator=(const LE_OneDbl &b) :
   *
   *  assignment operator
   */
  LE_OneDbl& LE_OneDbl::operator=(const LE_OneDbl &b) {
    if (&b != this) {
      LineEntry::operator=(b);
      AddrVal    = b.AddrVal;
      MaxVal     = b.MaxVal;
      MinVal     = b.MinVal;
      DefaultVal = b.DefaultVal;
      CurrValue  = b.CurrValue;
      strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
    }
    return *this;
  }

  /*
   *
   * LineEntry* duplMyselfAsLineEntry() (virtual)
   *
   * Duplication as a base class
   */
  LineEntry* LE_OneDbl::duplMyselfAsLineEntry() const {
    LE_OneDbl* newLE = new LE_OneDbl(*this);
    return (LineEntry*) newLE;
  }

  // Empty destructor
  LE_OneDbl::~LE_OneDbl() {
  }

  /*
   *
   *   This routine will print out a processed line
   *
   *  The default behavior is to print the original line with a "=>"
   *  prefix to indicate that action has been taken on it.
   */
  void
  LE_OneDbl::print_ProcessedLine(const TK_TOKEN *lineArgTok) const {
    LineEntry::print_ProcessedLine(lineArgTok);
    if (strlen(PrintString) > 0) {
      printf(" ====> %s = %g\n", PrintString, CurrValue);
    }
  }

  /*
   * print_usage() (virtual function)
   *
   */
  void LE_OneDbl::print_usage(int ilvl) const {
    print_indent(ilvl);
    printf("%s = (double)", EntryName.orig_str);
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
  }

  /*
   * Adjust the destination base address to store the input value 
   * (virtual function)
   */
  void LE_OneDbl::adjustAddress(LONG_PTR addrAdjustment) {
    if (AddrVal != 0) {
      if (addrAdjustment != 0) {
	LONG_PTR ll = reinterpret_cast<LONG_PTR>(AddrVal);
	ll += addrAdjustment;
	AddrVal = reinterpret_cast<double *>(ll);
      }
    }
  }

  double LE_OneDbl::currentTypedValue() const {
    return CurrValue;
  }

  /*
   * currentValueAsVoidP() (virtual)
   */
  const void * LE_OneDbl::currentValueAsVoidP() const {
    return static_cast<const void *>(&CurrValue);
  }
  //====================================================================================================================
  /*
   * set_default()
   *
   * Note a warning will be printed if the storred value in the 
   * destination address isn't equal to the default value input here.
   */
  void LE_OneDbl::set_default(double def)
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
        printf(" WARNING: LE_OneDbl set_default() for %s: Current Value, %g, doesn't agree "
   	       "with default value %g. Changing value at target address.\n",
	       PrintString, val, DefaultVal);
        *AddrVal = DefaultVal;
      }
    }
  }
  //====================================================================================================================
  /*
   * set_limits()
   *
   */
  void LE_OneDbl::set_limits(double maxV, double minV)
  {
    MaxVal = maxV;
    MinVal = minV;
  }

  /*
   *set_PrintString()
   *
   */
  void LE_OneDbl::set_PrintString(const char *ps)
  {
    if (ps) {
      strncpy(PrintString, ps, MAX_INPUT_STR_LN);
      PrintString[MAX_INPUT_STR_LN] = '\0';
    }
  }
}
