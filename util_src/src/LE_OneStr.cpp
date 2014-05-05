/**
 * @file LE_OneStr.cpp
 *   Header for the LineEntry of a single string
 *  (see \ref blockentryModule and class \link BEInput::LE_OneStr LE_OneStr\endlink).
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


#include "LE_OneStr.h"
#include "mdp_allo.h"

#include <string>
#include "climits"

using namespace mdpUtil;

namespace BEInput {

  /*
   * LE_OneStr Constructor:
   *
   *   This sets up the line entry special case.
   *   We make sure to call the base class constructor here to do
   *   much of the initialization. 
   *
   *   The maximum number of tokens allowed in the string is equal
   *   to the limit specified in the tok_input_util.h file.
   *   The minimum number of tokens is equal to one. In other words we
   *   require something in the string as a default.
   */
  LE_OneStr::LE_OneStr(const char *lineName, std::string *aaa, int maxval,
		       int minval, int numRL, const char *varName) :
    LineEntry(lineName, numRL),
    AddrVal(aaa),
    MaxVal(maxval),
    MinVal(minval)
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
  }

  /*
   *
   * LE_OneStr(const LE_OneStr&):
   *
   * copy Constructor:
   */
  LE_OneStr::LE_OneStr(const LE_OneStr& b) : 
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
   * LE_OneStr& operator=(const LE_OneStr &b) :
   *
   *  assignment operator
   */
  LE_OneStr& LE_OneStr::operator=(const LE_OneStr &b) {
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
   * LineEntry* duplMyselfAsLineEntry() (virtual)
   *
   * Duplication as a base class
   */
  LineEntry* LE_OneStr::duplMyselfAsLineEntry() const {
    LE_OneStr* newLE = new LE_OneStr(*this);
    return (LineEntry*) newLE;
  }

  /*
   * process_lineEntry():
   *
   *   Here we interpret the token as a single string, and then
   *   assign it to the address we set up.
   */
  void LE_OneStr::process_LineEntry(const TK_TOKEN *lineArgTok)
  {
    for (int i = 0; i < NumEntryDependencies; i++) {
      BI_Dependency *idep = EntryDepList[i];
      if (idep->ResultType() ==  BIDRT_USETOPROCESS) {
	int value;
	if (idep->checkDepOneInt(value)) {
	  MaxVal = value;
	  MinVal = value;
	} else {
	  throw BI_InputError("LE_OneStr::process_LineEntry",
			      "Dependency not satisfied");
	}
      }
      if (idep->ResultType() ==  BIDRT_PT_ERROR) {
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

    BOOLEAN error = 0;
    char *str = tok_to_string(lineArgTok, MaxVal, MinVal, DefaultVal.c_str(),
			      &error);
    if (error) {
      throw BI_InputError("LE_OneStr::process_LineEntry",
			  "tok_to_string interpretation");
    }
    if (AddrVal) {
      *AddrVal = str;
    }
    CurrValue = str;
    m_numTimesProcessed++;

    mdp_safe_free((void **) &str);
  }

  /*
   * LE_OneStr::print_ProcessedLine() (virtual function):
   *
   *   This routine will print out a processed line
   *
   *  The default behavior is to print the original line with a "=>"
   *  prefix to indicate that action has been taken on it.
   *  Then, we print out a further line stating we have
   *  set a variable
   *   ====> PrintString = currValue
   */
  void
  LE_OneStr::print_ProcessedLine(const TK_TOKEN *lineArgTok) const {
    LineEntry::print_ProcessedLine(lineArgTok);
    if (strlen(PrintString) > 0) {
      printf(" ====> %s = %s\n", PrintString, CurrValue.c_str());
    }
  }

  /*
   *
   * print_usage() (virtual function)
   *
   */
  void LE_OneStr::print_usage(int ilvl) const {
    print_indent(ilvl);
    printf("%s = (string)", EntryName.orig_str);
    if (MaxVal != INT_MAX || MinVal != -INT_MAX) {
      printf(" with Token limits(%d, %d)", MaxVal, MinVal);
    }
    if (DefaultVal != "") {
      printf(" with default %s", DefaultVal.c_str());
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
   *
   * Adjust base address to store the value (virtual function)
   *
   */
  void LE_OneStr::adjustAddress(LONG_PTR addrAdjustment) {
    if (AddrVal) {
      if (addrAdjustment != 0) {
	LONG_PTR ll = reinterpret_cast<LONG_PTR>(AddrVal);
	ll += addrAdjustment;
	AddrVal = reinterpret_cast<std::string *>(ll);
      }
    }
  }

  std::string LE_OneStr::currentTypedValue() const {
    return(std::string(CurrValue));
  }

  const void * LE_OneStr::currentValueAsVoidP() const {
    return static_cast<const void *>(&CurrValue);
  }
  //====================================================================================================================
  /*
   *
   * set_default()
   *   Set the default value for this entry
   */
  void LE_OneStr::set_default(std::string def)
  {
    DefaultVal = def;
    if  (m_numTimesProcessed == 0) {
      CurrValue = def;
      std::string aval = DefaultVal;
      if (AddrVal) {
        aval = *AddrVal;
      }
      if (aval != DefaultVal) {
        printf(" WARNING: LE_OneStr::set_default() for %s: Current target addr Value, %s, doesn't agree "
		 "with default value %s. Changing value at target address.\n",
		 PrintString, aval.c_str(), DefaultVal.c_str());
        *AddrVal = DefaultVal;
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
  void LE_OneStr::set_limits(int maxV, int minV)
  {
    MaxVal = maxV;
    MinVal = minV;
  }

  /*
   *
   * set_PrintString()
   *
   */
  void LE_OneStr::set_PrintString(const char *ps)
  {
    if (ps) {
      strncpy(PrintString, ps, MAX_INPUT_STR_LN);
      PrintString[MAX_INPUT_STR_LN] = '\0';
    }
  }

}
