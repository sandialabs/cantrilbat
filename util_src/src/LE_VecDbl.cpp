/**
 * @file LE_VecDbl.cpp
 *  Definitions for the LineEntry of a vector of  doubles
 *  (see \ref blockentryModule and class 
 *  \link BEInput::LE_VecDbl LE_VecDbl\endlink).
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

#include "LE_VecDbl.h"
#include "mdp_allo.h"

#include <cfloat>

using namespace TKInput;
using namespace mdpUtil;

namespace BEInput {

  /*
   * LE_VecDbl Constructor:
   *
   *   This sets up the line entry special case.
   *   We make sure to call the base class constructor here to do
   *   much of the initialization. 
   */
  LE_VecDbl::LE_VecDbl(const char *lineName, double **hndlAAA, 
		       int vecLength,
		       int numRL, const char *varName) :
    LineEntry(lineName, numRL),
    HndlDblVec(hndlAAA),
    m_CurrentVecValues(0)
  {
    /*
     * Set the limits and the default values
     */
    MaxVal = DBL_MAX;
    MinVal = -DBL_MAX;
    DefaultVal = NO_DEFAULT_DOUBLE;
    CurrValue = DefaultVal;
    VecLength = vecLength;
    CurrIndex = 0;
    PrintString[0] = '\0';
    if (varName) {
      strncpy(PrintString, varName, MAX_INPUT_STR_LN);
    } else {
      strncpy(PrintString, lineName, MAX_INPUT_STR_LN);
    }
    PrintString[MAX_INPUT_STR_LN] = '\0';
    if (numRL != 0 && numRL != 1) {
      throw BI_InputError("LE_VecDbl::constructor",
			  "numRL must be 0 or 1");     
    }
    if (vecLength <= 0) {
      throw BI_InputError("LE_VecDbl::constructor",
			  "vecLength must be greater than 0");  
    }
    if (*HndlDblVec == 0 && vecLength > 0) {
      *HndlDblVec = mdp_alloc_dbl_1(vecLength, DefaultVal);
    }
    m_CurrentVecValues =  mdp_alloc_dbl_1(vecLength, DefaultVal);
  }


  /*
   * LE_VecDbl(const LE_VecDbl&):
   *
   * copy Constructor:
   *  -> we set the CurrIndex to zero here. And, note, we still
   *     address and write to the same memory locations.
   */
  LE_VecDbl::LE_VecDbl(const LE_VecDbl& b) : 
    LineEntry(b),
    HndlDblVec(b.HndlDblVec),
    MaxVal(b.MaxVal),
    MinVal(b.MinVal),
    DefaultVal(b.DefaultVal),
    VecLength(b.VecLength),
    CurrIndex(0),
    CurrValue(b.CurrValue),
    m_CurrentVecValues(0)
  {   
    strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
    mdp_realloc_dbl_1(&m_CurrentVecValues, VecLength, 0, DefaultVal);
    mdp_copy_dbl_1(b.m_CurrentVecValues, m_CurrentVecValues,
		   VecLength);
  }


  /*
   * LE_VecDbl& operator=(const LE_VecDbl&) :
   *
   *  assignment operator
   */
  LE_VecDbl& LE_VecDbl::operator=(const LE_VecDbl& b)
  {  
    if (&b != this) {
      LineEntry::operator=(b);
      HndlDblVec = b.HndlDblVec;
      MaxVal  = b.MaxVal;
      MinVal  = b.MinVal;
      DefaultVal = b.DefaultVal;
      CurrIndex = b.CurrIndex;
      CurrValue = b.CurrValue;
      strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
      mdp_realloc_dbl_1(&m_CurrentVecValues, VecLength, 0, DefaultVal);
      mdp_copy_dbl_1(b.m_CurrentVecValues, m_CurrentVecValues,
		     VecLength);
    }
    return *this;
  }
  LE_VecDbl::~LE_VecDbl() {
    mdp_safe_free((void **) &m_CurrentVecValues);
  }
  
  /*
   * LineEntry* duplMyselfAsLineEntry();   (virtual)
   *
   *  Duplicate myself as the base class
   */
  LineEntry* LE_VecDbl::duplMyselfAsLineEntry() const {
    LE_VecDbl* newLE = new LE_VecDbl(*this);
    return (LineEntry*) newLE;
  }

  /*
   * process_lineEntry():
   *
   *   Here we interpret the token as a single double, and then
   *   assign it to the address we set up.
   */
  void LE_VecDbl::process_LineEntry(const TK_TOKEN *lineArgTok)
  {
    int i;
    std::string ename = EntryName.orig_str;
    for (i = 0; i < NumEntryDependencies; i++) {
      BI_Dependency *idep = EntryDepList[i];
      std::string dname = idep->TargetBaseEntry()->keyname();
      if (idep->ResultType() == BIDRT_USETOPROCESS) {
	int value;
	if (idep->checkDepOneInt(value)) {
	  if (value > VecLength) {
	    mdp_realloc_dbl_1(HndlDblVec, value, VecLength, DefaultVal);
	    mdp_realloc_dbl_1(&m_CurrentVecValues, value,
			      VecLength, DefaultVal);
	  }
	  VecLength = value;
	} else {
	  throw BI_InputError("LE_VecDbl::process_LineEntry on \"" + ename + "\"",
			      "Dependency not satisfied  wrt keyname \"" + dname + "\"");
	}
      }
      if (idep->ResultType() == BIDRT_PT_ERROR) {
	if (!idep->checkDependencySatisfied()) {
	  throw BI_InputError("LE_VecDbl::process_LineEntry on \"" + ename + "\"",
			      "Dependency not satisfied wrt keyname \"" + dname + "\"");
	}
      }
      if (idep->ResultType() ==  BIDRT_ANTITHETICAL_ERROR) {
	if (idep->checkDependencySatisfied()) {
	  throw BI_InputError("LE_VecDbl::process_LineEntry on \"" + ename + "\"",
			      "Mutual Exclusive Dependency is satisfied  wrt keyname \"" + dname + "\"");
	}
      }
    }

    BOOLEAN error = 0;
    double *AddrVal = *HndlDblVec;
    for (i = 0; i < lineArgTok->ntokes; i++) {
      char *strVal = lineArgTok->tok_ptrV[i];
      double value = str_to_double(strVal, MaxVal, MinVal, DefaultVal, 
				   &error);
      if (error) {
	throw BI_InputError("LE_VecDbl::process_LineEntry on \"" + ename + "\"",
			    "str_to_double interpretation");
      }
      if (CurrIndex >= VecLength) {
	throw BI_InputError("LE_VecDbl::process_LineEntry on \"" + ename + "\"",
			    "Too many entries");	
      }
      *(AddrVal + CurrIndex) = value;
      CurrIndex++;
      CurrValue = value;
      m_CurrentVecValues[CurrIndex] = CurrValue;
    }
    if (CurrIndex == VecLength) {
      m_numTimesProcessed++;
      CurrIndex = 0;
    }

  }

  /*
   * LE_VecDbl::print_ProcessedLine() (virtual function):
   *
   *   This routine will print out a processed line
   *
   *  The default behavior is to print the original line with a "=>"
   *  prefix to indicate that action has been taken on it.
   */
  void
  LE_VecDbl::print_ProcessedLine(const TK_TOKEN *lineArgTok) const {
    printf(" => %s", EntryName.orig_str);
    if (strlen(PrintString) > 0) {
      printf(" ====> %s = %g\n", PrintString, CurrValue);
    }
  }

  /*
   * set_PrintString()
   *
   */
  void LE_VecDbl::set_PrintString(const char *ps)
  {
    if (ps) {
      strncpy(PrintString, ps, MAX_INPUT_STR_LN);
      PrintString[MAX_INPUT_STR_LN] = '\0';
    }
  }

  /*
   * LE_VecDbl::adjustAddress() : (virtual)
   *
   * Adjust base address to store the value (virtual function)
   * Note, the address adjustment is done in byte units, and not
   * via pointer arithmetic.
   */
  void LE_VecDbl::adjustAddress(LONG_PTR addrAdjustment) {
    if (addrAdjustment != 0) {
      LONG_PTR ll = reinterpret_cast<LONG_PTR>(HndlDblVec);
      ll += addrAdjustment;
      HndlDblVec = reinterpret_cast<double **>(ll);
    }
  }

  const double * LE_VecDbl::currentTypedValue() const {
    return m_CurrentVecValues;
  }

  /*
   * currentValueAsVoidP() (virtual)
   */
  const void * LE_VecDbl::currentValueAsVoidP() const {
    return static_cast<const void *>(m_CurrentVecValues);
  }

  /*
   * checkRequirements(): (virtual function)
   *
   * Check for requirements being met at the end of input
   *
   *  Defaults are specified in the .h file. They are:
   *   throwSpecificError = false
   *
   */
  bool LE_VecDbl::checkRequirements(bool throwSpecificError) 
  {
    if (m_numTimesRequired) {
      if (m_numTimesRequired != m_numTimesProcessed) {
	return false;
      }
    }
    if (CurrIndex != 0) {
      if (throwSpecificError) {
	throw BI_MissingRequiredVec("LE_VecDbl::checkRequirements",
				    EntryName.orig_str,
				    VecLength,
				    CurrIndex);
      } else {
	return false;
      }
    }
    return true;
  }

  /*
   *
   * print_usage() (virtual function)
   *
   */
  void LE_VecDbl::print_usage(int ilvl) const {
    print_indent(ilvl);
    printf("%s = double1 double2\n", EntryName.orig_str);
    print_indent(ilvl);
    printf("%s = (double3 ... double_VecLength", EntryName.orig_str);
    printf(" where VecLength = %d ", VecLength);
    if (MaxVal != DBL_MAX || MinVal != -DBL_MAX) {
      printf(" with limits (%g, %g)", MinVal, MaxVal);
    }
    if (DefaultVal != NO_DEFAULT_DOUBLE) {
      printf(" with default %g", DefaultVal);
    }
    if (m_numTimesRequired > 1) {
      printf(" (REQUIRED LINE)");
    }
    if (m_numTimesRequired > 1) {
      printf(" (REQUIRED LINES = %d)", m_numTimesRequired);
    }
    printf("\n");
  }

  /**********************************************************************
   *
   * set_default()
   *
   */
  void LE_VecDbl::set_default(double def)
  {
    DefaultVal = def; 
    if (m_numTimesProcessed == 0) {
      CurrValue = def;
    }
  }

  /**********************************************************************
   *
   * set_limits()
   *
   */
  void LE_VecDbl::set_limits(double maxV, double minV)
  {
    MaxVal = maxV;
    MinVal = minV;
  }

}
