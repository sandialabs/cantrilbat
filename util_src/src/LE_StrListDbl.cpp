/**
 * @file LE_StrListDbl.cpp
 *  Definitions for the LineEntry of a list of doubles
 *  (see \ref blockentryModule and class 
 *  \link BEInput::LE_StrListDbl LE_StrListDbl\endlink).
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

#include "LE_StrListDbl.h"
#include "mdp_allo.h"

#include "stdio.h"
#include "string.h"

using namespace TKInput;
using namespace mdpUtil;

namespace BEInput {

  /*
   * LE_StrListDbl Constructor:
   *
   *   This sets up the line entry special case.
   *   We make sure to call the base class constructor here to do
   *   much of the initialization. 
   */
  LE_StrListDbl::LE_StrListDbl(const char *lineName, double **hndlAAA, 
			       char **charList, int listLength, int numRL,
			       const char *varName) :
    LineEntry(lineName, numRL),
    HndlDblVec(hndlAAA),
    m_CharList(0),
    m_ListLength(listLength),
    m_CurrentVecValues(0)
  {
    char *item(0);
    /*
     * Set the limits and the default values
     */
    MaxVal = m_ListLength - 1;
    MinVal = 0;
    DefaultVal = NO_DEFAULT_DOUBLE;
    CurrValue = DefaultVal;
    CurrListValue = NO_DEFAULT_INT;
    if (charList) {
      m_CharList = 
	mdp_alloc_VecFixedStrings(m_ListLength, MAX_INPUT_STR_LN+1);
      for (int i = 0; i < m_ListLength; i++) {
	item = charList[i];
	if (!item) throw BI_InputError("LE_StrListDbl::LE_StrListDbl",
				       "char list item is null");
	strncpy(m_CharList[i], item, MAX_INPUT_STR_LN);
	m_CharList[i][MAX_INPUT_STR_LN] = '\0';
      }
    }
    PrintString[0] = '\0';
    if (varName) {
      strncpy(PrintString, varName, MAX_INPUT_STR_LN);
    } else {
      strncpy(PrintString, lineName, MAX_INPUT_STR_LN);
    }
    PrintString[MAX_INPUT_STR_LN] = '\0';
    // If we have been given an address to store the vector
    // external to the object, make sure it is malloced.
    if (HndlDblVec) {
      if (*HndlDblVec == 0 && m_ListLength > 0) {
	*HndlDblVec = mdp_alloc_dbl_1(m_ListLength, DefaultVal);
      }
    }
    // Always store a complete list of values.
    m_CurrentVecValues = mdp_alloc_dbl_1(m_ListLength, DefaultVal);
  }

  /*
   * LE_StrListDbl(const LE_StrListDbl&) :
   *
   *   copy constructor
   */
  LE_StrListDbl::LE_StrListDbl(const LE_StrListDbl& b) :
    LineEntry(b),
    HndlDblVec(b.HndlDblVec),
    MaxVal(b.MaxVal),
    MinVal(b.MinVal),
    DefaultVal(b.MinVal),
    m_CharList(0),
    m_ListLength(b.m_ListLength),
    CurrListValue(b.CurrListValue),
    CurrValue(b.MinVal),
    m_CurrentVecValues(0)
  {
    if (m_ListLength > 0) {
      char *item;
      m_CharList = 
	mdp_alloc_VecFixedStrings(m_ListLength, MAX_INPUT_STR_LN+1);
      for (int i = 0; i < m_ListLength; i++) {
	item = b.m_CharList[i];
	if (!item) throw BI_InputError("LE_StrListDbl::LE_StrListDbl",
				       "char list item is null");
	strncpy(m_CharList[i], item, MAX_INPUT_STR_LN);
	m_CharList[i][MAX_INPUT_STR_LN] = '\0';
      }
    }
    strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
    m_CurrentVecValues = mdp_alloc_dbl_1(m_ListLength, DefaultVal);
    mdp_copy_dbl_1(m_CurrentVecValues, b.m_CurrentVecValues, m_ListLength);
  }

  /*
   * LE_StrListDbl& operator=(const LE_StrListDbl &b) :
   *
   *  assignment operator
   */
  LE_StrListDbl& LE_StrListDbl::operator=(const LE_StrListDbl &b) {
    if (&b != this) {
      LineEntry::operator=(b);
      HndlDblVec = b.HndlDblVec;
      MaxVal     = b.MaxVal;
      MinVal     = b.MinVal;
      DefaultVal = b.DefaultVal;
      m_ListLength = b.m_ListLength;
      CurrListValue = b.CurrListValue;
      CurrValue  = b.CurrValue;
      strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
      mdp_safe_free((void **) &m_CharList);
      mdp_safe_free((void **) &m_CurrentVecValues);
      if (m_ListLength > 0) {
	char *item;
	m_CharList = 
	  mdp_alloc_VecFixedStrings(m_ListLength, MAX_INPUT_STR_LN+1);
	for (int i = 0; i < m_ListLength; i++) {
	  item = b.m_CharList[i];
	  if (!item) throw BI_InputError("LE_StrListDbl::LE_StrListDbl",
					 "char list item is null");
	  strncpy(m_CharList[i], item, MAX_INPUT_STR_LN);
	  m_CharList[i][MAX_INPUT_STR_LN] = '\0';
	}

	m_CurrentVecValues = mdp_alloc_dbl_1(m_ListLength, DefaultVal);
	mdp_copy_dbl_1(m_CurrentVecValues, b.m_CurrentVecValues, m_ListLength);
      }
    }
    return *this;
  }

  /*
   * LineEntry* duplMyselfAsLineEntry() (virtual)
   *
   * Duplication as a base class
   */
  LineEntry* LE_StrListDbl::duplMyselfAsLineEntry() const {
    LE_StrListDbl* newLE = new LE_StrListDbl(*this);
    return (LineEntry*) newLE;
  }

  /*
   * LE_StrListDbl destructor:
   *
   * We malloced memory here, so we must explicitly call free.
   */
  LE_StrListDbl::~LE_StrListDbl()
  {
    mdp_safe_free((void **) &m_CharList);
    mdp_safe_free((void **) &m_CurrentVecValues);
  }

  /*
   * process_lineEntry():
   *
   *  We expect two tokens. The first token is used to look up the
   *  index in the list (case sensitive search).  The second token 
   *  is interpreted as a double and then the index address is assigned
   *  that value. 
   */
  void LE_StrListDbl::
  process_LineEntry(const TK_TOKEN *lineArgTok)
  {
    for (int i = 0; i < NumEntryDependencies; i++) {
      BI_Dependency *idep = EntryDepList[i];
      if (idep->ResultType() == BIDRT_PT_ERROR) {
	if (!idep->checkDependencySatisfied()) {
	  throw BI_InputError("LE_StrListDbl::process_LineEntry",
			      "Dependency not satisfied");
	}
      }
      if (idep->ResultType() ==  BIDRT_ANTITHETICAL_ERROR) {
	if (idep->checkDependencySatisfied()) {
	  throw BI_InputError("LE_StrListDbl::process_LineEntry",
			      "Mutual Exclusive Dependency is satisfied");
	}
      }
    }

    BOOLEAN error = 0;
    char *cvalue = 0;
    if (lineArgTok->ntokes < 2) {
      throw BI_InputError("LE_StrListDbl::process_LineEntry",
			  "Need at least two tokens");
    }
    cvalue = lineArgTok->tok_ptrV[0];
    int ifind = -1;
    for (int i = 0; i < m_ListLength; i++) {
      if (!strncmp(cvalue, m_CharList[i], MAX_INPUT_STR_LN)) {
	ifind = i;
	break;
      }
    }
    if (ifind >= 0) {
      CurrListValue = ifind;
    } else {
      throw BI_UnknownListEntry("LE_StrListDbl::process_LineEntry", cvalue);
    }
    double value = str_to_double(lineArgTok->tok_ptrV[1], MaxVal, MinVal, 
				 DefaultVal, &error);
    if (error) {
      throw BI_InputError("LE_StrListDbl::process_LineEntry",
			  "tok_to_double error");
    }
    if (HndlDblVec) {
      if (*HndlDblVec == 0) {
	*HndlDblVec = mdp_alloc_dbl_1(m_ListLength, DefaultVal);
      }
      double *AddrVal = *HndlDblVec;
      AddrVal[CurrListValue] = value;
    }
    CurrValue = value;
    m_CurrentVecValues[CurrListValue] = value;
    m_numTimesProcessed++;
  }

  /*
   * LE_StrListDbl::print_ProcessedLine() (virtual function):
   *
   *   This routine will print out a processed line
   *
   *  The default behavior is to print the original line with a "=>"
   *  prefix to indicate that action has been taken on it.
   */
  void
  LE_StrListDbl::print_ProcessedLine(const TK_TOKEN *lineArgTok) const {
    LineEntry::print_ProcessedLine(lineArgTok);
    if (strlen(PrintString) > 0) {
      printf(" ====> %s[ %s : %d] = %g\n", PrintString, 
	     m_CharList[CurrListValue], CurrListValue, CurrValue);
    }
  }

  /*
   * print_usage() (virtual function)
   *
   */
  void LE_StrListDbl::print_usage(int ilvl) const {
    print_indent(ilvl);
    printf("%s = (LIST) (double)", EntryName.orig_str);
    if (m_numTimesRequired == 1) {
      printf(" (REQUIRED LINE)");
    }
    if (m_numTimesRequired > 1) {
      printf(" (REQUIRED LINES = %d)", m_numTimesRequired);
    }
    printf(" where:\n");
    print_indent(ilvl+1); printf("LIST = \n");
    for (int i = 0; i < m_ListLength; i++) {
      print_indent(ilvl+2);
      printf("%s : %d\n", m_CharList[i], i);
    }
  }

  /*
   * LE_StrListDbl::adjustAddress(): (virtual)
   *
   * Adjust base address of the vector of doubles.
   * Note, the address adjustment is done in byte units, and not
   * via pointer arithmetic.
   */
  void LE_StrListDbl::adjustAddress(LONG_PTR addrAdjustment) {
    if (addrAdjustment != 0) {
      LONG_PTR ll = reinterpret_cast<LONG_PTR>(HndlDblVec);
      ll += addrAdjustment;
      HndlDblVec = reinterpret_cast<double **>(ll);
    }
  }

  const double * LE_StrListDbl::currentTypedValue() const {
    return m_CurrentVecValues;
  }

  /*
   * currentValueAsVoidP() (virtual)
   */
  const void * LE_StrListDbl::currentValueAsVoidP() const {
    return static_cast<const void *>(m_CurrentVecValues);
  }

  /*
   * set_default()
   *
   */
  void LE_StrListDbl::set_default(double def) {
    DefaultVal = def;
  }

  /*
   * set_limits()
   *
   */
  void LE_StrListDbl::set_limits(double maxV, double minV) {
    MaxVal = maxV;
    MinVal = minV;
  }

  /*
   * set_PrintString()
   *
   */
  void LE_StrListDbl::set_PrintString(const char *ps) {
    if (ps) {
      strncpy(PrintString, ps, MAX_INPUT_STR_LN);
      PrintString[MAX_INPUT_STR_LN] = '\0';
    }
  }

}
