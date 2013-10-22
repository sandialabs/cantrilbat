/**
 * @file LE_OneInt.cpp
 *  Definitions for the LineEntry of a single integer
 *  (see \ref blockentryModule and class 
 *  \link BEInput::LE_OneInt LE_OneInt\endlink).
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
#include "LE_OneInt.h"
#include "mdp_allo.h"

#include "stdio.h"
#include "string.h"
#include "limits.h"
#include "float.h"

namespace BEInput {

  /*
   *   Entry of a single integer value.
   *
   *   This sets up the line entry special case for the entry of a single
   *   int input.
   *
   *   @param lineName c-string containing the keyline id.
   *   @param aaa      Address of the location where the input int 
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
   *         int numberIterations;
   *      };
   *
   *      BlockEntry *besmd;
   *      int reqd = 1;
   *      LE_OneInt *i2 = new LE_OneInt("Number of Iterations",
   *                                     &globInput.numberIterations, 
   *                                     reqd, "numberIterations");
   *      i2->set_default(10.);
   *      i2->set_limits(1000, 0);
   *      besmd->addLineEntry(di);
   *
   *  Entry:
   *       Number of Iterations = 5
   *
   *  Running:
   *     besmd->read_block(fp *inputfile);
   *
   *  Dependencies that may be set on this card that makes sense.
   *          BIDRT_PT_ERROR
   *          BIDRT_ANTITHETICAL_ERROR
   *          BIDRT_ZERONUMTIMESREQUIRED
   *          BIDRT_ONENUMTR
   *          BIDRT_USETOPROCESS   If set on this BaseEntry,
   *                 then the required BaseEntry must have been already
   *                 processed, and the default value for this card
   *                 is set by the value of the required BaseEntry.
   *          BIDRT_RTINTMM_ERROR 
   *                 if the CurrentValue of this entry is 
   *                 in a specified range of values, then a 
   *                 dependency check is made against the required 
   *                 BaseEntry.
   */
  LE_OneInt::LE_OneInt(const char *lineName, int *aaa, int numRL,
		       const char *varName) :
    LineEntry(lineName, numRL),
    AddrVal(aaa)
  {
    /*
     * Set the limits and the default values
     */
    MaxVal = INT_MAX;
    MinVal = -INT_MAX;
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

  /*
   *
   * LE_OneInt(const LE_OneInt&) :
   *
   *   copy constructor
   */
  LE_OneInt::LE_OneInt(const LE_OneInt& b) :
    LineEntry(b),
    AddrVal(b.AddrVal),
    MaxVal(b.MaxVal),
    MinVal(b.MinVal),
    DefaultVal(b.MinVal),
    CurrValue(b.MinVal)
  {
    strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
  }

  /*
   *
   * LE_OneInt& operator=(const LE_OneInt&) :
   *
   *  assignment operator
   */
  LE_OneInt& LE_OneInt::operator=(const LE_OneInt& b)
  {  
    if (&b != this) {
      LineEntry::operator=(b);
      AddrVal = b.AddrVal;
      MaxVal  = b.MaxVal;
      MinVal  = b.MinVal;
      DefaultVal = b.DefaultVal;
      CurrValue = b.CurrValue;
      strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
    }
    return *this;
  }

  /*
   *
   * LineEntry* duplMyselfAsLineEntry();   (virtual)
   *
   *  Duplicate myself as the base class
   */
  LineEntry* LE_OneInt::duplMyselfAsLineEntry() const {
    LE_OneInt* newLE = new LE_OneInt(*this);
    return (LineEntry*) newLE;
  }

  /*
   *
   * process_lineEntry():
   *
   *   Here we interpret the token as a single integer, and then
   *   assign it to the address we set up.
   */
  void LE_OneInt::
  process_LineEntry(const TK_TOKEN *lineArgTok)
  {

    /*
     * Handle common dependencies and incr NumTimesProcessed.
     */
    LineEntry::process_LineEntry(lineArgTok);

    for (int i = 0; i < NumEntryDependencies; i++) {
      BI_Dependency *idep = EntryDepList[i];
      if (idep->ResultType() == BIDRT_USETOPROCESS) {
	int value;
	if (idep->checkDepOneInt(value)) {
	  DefaultVal = value;
	} else {
	  throw BI_InputError("LE_OneInt::process_LineEntry",
			      "Dependency not satisfied");
	}
      }
    }

    BOOLEAN error = 0;
    int value = tok_to_int(lineArgTok, MaxVal, MinVal, DefaultVal, 
			   &error);
    if (error) {
      throw BI_InputError("LE_OneInt::process_LineEntry",
			  "tok_to_int interpretation");
    }
    if (AddrVal) {
      *AddrVal = value;
    }
    CurrValue = value;

    /*
     * Now that we have the value, check to see if any dependencies based
     * on that value are broken
     */
    for (int i = 0; i < NumEntryDependencies; i++) {
      BI_Dependency *idep = EntryDepList[i];
      if (idep->ResultType() == BIDRT_RTINTMM_ERROR) {
	if (idep->checkCardIntMaxMin(CurrValue)) {
	  if (!idep->checkDependencySatisfied()) {
	    throw BI_InputError("LE_OneInt::process_LineEntry",
				"Triggered dependency check from the local input"
				" resulted in a dependency failure");
	  }
	} 
      }
    }
  }

  /*
   *
   * LE_OneInt::print_ProcessedLine() (virtual function):
   *
   *   This routine will print out a processed line
   *
   *  The default behavior is to print the original line with a "=>"
   *  prefix to indicate that action has been taken on it.
   */
  void
  LE_OneInt::print_ProcessedLine(const TK_TOKEN *lineArgTok) const {
    printf(" => %s", EntryName.orig_str);
    if (strlen(PrintString) > 0) {
      printf(" ====> %s = %d\n", PrintString, CurrValue);
    }
  }

  /*
   * print_usage() (virtual function)
   *
   */
  void LE_OneInt::print_usage(int ilvl) const {
    print_indent(ilvl);
    printf("%s = (integer)", EntryName.orig_str);
    if (MaxVal != INT_MAX || MinVal != -INT_MAX) {
      printf(" with limits (%d, %d)", MinVal, MaxVal);
    }
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
    if (m_numTimesProcessed) {
      print_indent(ilvl+1);
      print_ProcessedLine(0);
    }
  }


  /*
   *
   * LE_OneInt::adjustAddress() : (virtual)
   *
   * Adjust base address to store the value (virtual function)
   * Note, the address adjustment is done in byte units, and not
   * via pointer arithmetic.
   */
  void LE_OneInt::adjustAddress(LONG_PTR addrAdjustment) {
    if (AddrVal) {
      if (addrAdjustment != 0) {
	LONG_PTR ll = reinterpret_cast<LONG_PTR>(AddrVal);
	ll += addrAdjustment;
	AddrVal = reinterpret_cast<int *>(ll);
      }
    }
  }

  int LE_OneInt::currentTypedValue() const {
    return(CurrValue);
  }

  const void * LE_OneInt::currentValueAsVoidP() const {
    return static_cast<const void *>(&CurrValue);
  }
  //====================================================================================================================
  /*
   *
   * set_default()
   *
   */
  void LE_OneInt::set_default(int def)
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
      CurrValue = DefaultVal;
      if (aval != DefaultVal) {
        printf(" WARNING: LE_Int::set_default() for %s: Current Value, %d, doesn't agree "
  	       "with default value %d. Changing value at target address.\n",
	       PrintString, aval, DefaultVal);
	if (AddrVal) {
	  *AddrVal = DefaultVal;
	}
      }
    }
  }
  //====================================================================================================================
  /*
   *
   * set_limits()
   *
   */
  void LE_OneInt::set_limits(int maxV, int minV)
  {
    MaxVal = maxV;
    MinVal = minV;
  }

  /*
   * set_PrintString()
   *
   */
  void LE_OneInt::set_PrintString(const char *ps)
  {
    if (ps) {
      strncpy(PrintString, ps, MAX_INPUT_STR_LN);
      PrintString[MAX_INPUT_STR_LN] = '\0';
    }
  }

  /*
   *
   * DepCanService():
   *
   *  This call checks to whether requests can be fulfilled by this
   *  Line Entry object. Basically, this object can return whether it
   *  has been called enough times, and it can return the integer
   *  value that was last processed.
   */
  bool LE_OneInt::DepCanService(BIDSR_TYPE BIDSR_value) const {
    bool retn = BaseEntry::DepCanService(BIDSR_value);
    if (!retn) {
      if (BIDSR_value == BIDSR_ONEINT) return true;
    }
    return retn;
  }

  /*
   *
   * ansDepCheckOntInt(): (virtual)
   *
   * This call returns a bool indicating whether it has been called before.
   * And if it has been, it returns the last value processed.
   */
  bool LE_OneInt::ansDepCheckOneInt(int &returnInt) const {
    bool retn = ansDepCheck();
    if (!retn) {
      returnInt = DefaultVal;
      return retn;
    }
    returnInt = CurrValue;
    return retn;
  }

}
