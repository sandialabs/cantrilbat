/**
 * @file LE_StdVecDblVarLength.cpp
 *  Definitions for the LineEntry of a vector of doubles
 *  (see \ref blockentryModule and class 
 *  \link BEInput::LE_StdVecDblVarLength LE_StdVecDblVarLength\endlink).
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

#include "LE_StdVecDblVarLength.h"

#include <cstring>
#include <cfloat>

#ifndef DATA_PTR
#define DATA_PTR(x)  &(x[0])
#endif


using std::vector;

namespace BEInput {

  /*
   * LE_StdVecDblVarLength Constructor:
   *
   *   This sets up the line entry special case.
   *   We make sure to call the base class constructor here to do
   *   much of the initialization. 
   */
  LE_StdVecDblVarLength::
  LE_StdVecDblVarLength(const char *lineName, std::vector<double> * const hndlAAA, 
			int vecLength, int numRL, const char *varName) :
    LineEntry(lineName, numRL),
    m_vectorDest_ptr(hndlAAA),
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
    if (vecLength <= 0) {
      throw BI_InputError("LE_StdVecDblVarLength::constructor",
			  "vecLength must be greater than 0");  
    }
    if (m_vectorDest_ptr) {
      if ((int) (*m_vectorDest_ptr).size() < vecLength) {
        (*m_vectorDest_ptr).resize(vecLength, DefaultVal);
      }
    }
    m_CurrentVecValues.resize(VecLength, DefaultVal);
  }

  /*
   * LE_StdVecDblVarLength(const LE_StdVecDblVarLength&):
   *
   * copy Constructor:
   *  -> we set the CurrIndex to zero here. And, note, we still
   *     address and write to the same memory locations.
   */
  LE_StdVecDblVarLength::LE_StdVecDblVarLength(const LE_StdVecDblVarLength& b) : 
    LineEntry(b),
    m_vectorDest_ptr(b.m_vectorDest_ptr),
    MaxVal(b.MaxVal),
    MinVal(b.MinVal),
    DefaultVal(b.DefaultVal),
    VecLength(b.VecLength),
    CurrIndex(0),
    CurrValue(b.CurrValue),
    m_CurrentVecValues(0)
  {   
    strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
    m_CurrentVecValues = b.m_CurrentVecValues;
  }

  /*
   * LE_StdVecDblVarLength& operator=(const LE_StdVecDblVarLength&) :
   *
   *  assignment operator
   */
  LE_StdVecDblVarLength&
  LE_StdVecDblVarLength::operator=(const LE_StdVecDblVarLength& b)
  {  
    if (&b != this) {
      LineEntry::operator=(b);
      m_vectorDest_ptr = b.m_vectorDest_ptr;
      MaxVal  = b.MaxVal;
      MinVal  = b.MinVal;
      DefaultVal = b.DefaultVal;
      VecLength = b.VecLength;
      CurrIndex = b.CurrIndex;
      CurrValue = b.CurrValue;
      strncpy(PrintString, b.PrintString, MAX_INPUT_STR_LN+1);
      m_CurrentVecValues = b.m_CurrentVecValues;
    }
    return *this;
  }

  // Destructor
  LE_StdVecDblVarLength::~LE_StdVecDblVarLength() {
  }

  /*
   * LineEntry* duplMyselfAsLineEntry();   (virtual)
   *
   *  Duplicate myself as the base class
   */
  LineEntry* LE_StdVecDblVarLength::duplMyselfAsLineEntry() const {
    LE_StdVecDblVarLength* newLE = new LE_StdVecDblVarLength(*this);
    return (LineEntry*) newLE;
  }

  //===============================================================================================================
  void LE_StdVecDblVarLength::zeroTimesCounter()
  {
    CurrIndex = 0;
    LineEntry::zeroTimesCounter();
  }
  //===============================================================================================================
  /*
   * process_lineEntry():
   *
   *   Here we interpret the token as a string of a variable number
   * of  doubles, and then  assign them sequentially  to a vector of 
   * doubles.
   */
  void LE_StdVecDblVarLength::process_LineEntry(const TK_TOKEN *lineArgTok)
  {
    int i;
    for (i = 0; i < NumEntryDependencies; i++) {
      BI_Dependency *idep = EntryDepList[i];
      if (idep->ResultType() ==  BIDRT_USETOPROCESS) {
	int value;
	if (idep->checkDepOneInt(value)) {
	  if (value > VecLength) {
	    if (m_vectorDest_ptr) {
	      (*m_vectorDest_ptr).resize(value, DefaultVal);
	    }
            m_CurrentVecValues.resize(value, DefaultVal);
	  }
	  VecLength = value;
	} else {
	  throw BI_InputError("LE_StdVecDblVarLength::process_LineEntry",
			      "Dependency not satisfied");
	}
      }
      if (idep->ResultType() ==  BIDRT_PT_ERROR) {
	if (!idep->checkDependencySatisfied()) {
	  throw BI_InputError("LE_StdVecDblVarLength::process_LineEntry",
			      "Dependency not satisfied");
	}
      }

      if (idep->ResultType() ==  BIDRT_ANTITHETICAL_ERROR) {
	if (idep->checkDependencySatisfied()) {
	  throw BI_InputError("LE_StdVecDblVarLength::process_LineEntry",
			      "Mutual Exclusive Dependency is satisfied");
	}
      }
    }

    BOOLEAN error = 0;
    double *AddrVal = 0;
    if (m_vectorDest_ptr) {
      AddrVal = DATA_PTR((*m_vectorDest_ptr));
    }

    for (i = 0; i < lineArgTok->ntokes; i++) {
      char *strVal = lineArgTok->tok_ptrV[i];
      double value = TKInput::str_to_double(strVal, MaxVal, 
					    MinVal, DefaultVal, 
					    &error);
      if (error) {
	throw BI_InputError("LE_StdVecDblVarLength::process_LineEntry",
			    "str_to_double interpretation");
      }
      if (CurrIndex >= VecLength) {
	if (m_vectorDest_ptr) {
	  (*m_vectorDest_ptr).resize(VecLength + lineArgTok->ntokes, DefaultVal);
	  AddrVal =  DATA_PTR((*m_vectorDest_ptr));;
	}
	m_CurrentVecValues.resize(VecLength + lineArgTok->ntokes, DefaultVal);
	VecLength += lineArgTok->ntokes;
      }
      if (AddrVal) {
	*(AddrVal + CurrIndex) = value;
      }
      m_CurrentVecValues[CurrIndex] = value;
      CurrIndex++;
      CurrValue = value;
    }
 
    m_numTimesProcessed = CurrIndex;
  }

  /*
   * LE_StdVecDblVarLength::print_ProcessedLine() (virtual function):
   *
   *   This routine will print out a processed line
   *
   *  The default behavior is to print the original line with a "=>"
   *  prefix to indicate that action has been taken on it.
   */
  void
  LE_StdVecDblVarLength::print_ProcessedLine(const TK_TOKEN *lineArgTok) const {
    printf(" => %s", EntryName.orig_str);
    if (strlen(PrintString) > 0) {
      printf(" ====> %s =", PrintString);
      for (int i = 0; i < (int) m_CurrentVecValues.size(); i++) {
        printf(" %g", m_CurrentVecValues[i]);
        if (i < (int) m_CurrentVecValues.size() - 1) {
          printf(",");
        }
      }
      printf("\n");
    }
  }

  /*
   * set_PrintString()
   *
   */
  void LE_StdVecDblVarLength::set_PrintString(const char *ps)
  {
    if (ps) {
      strncpy(PrintString, ps, MAX_INPUT_STR_LN);
      PrintString[MAX_INPUT_STR_LN] = '\0';
    }
  }

  /*
   * LE_StdVecDblVarLength::adjustAddress() : (virtual)
   *
   * Adjust base address to store the value (virtual function)
   * Note, the address adjustment is done in byte units, and not
   * via pointer arithmetic.
   */
  void LE_StdVecDblVarLength::adjustAddress(LONG_PTR addrAdjustment) {
    if (addrAdjustment != 0) {
      LONG_PTR ll = reinterpret_cast<LONG_PTR>(m_vectorDest_ptr);
      ll += addrAdjustment;
      m_vectorDest_ptr = reinterpret_cast<std::vector<double> *>(ll);
    }
  }

  const std::vector<double> & LE_StdVecDblVarLength::currentTypedValue() const {
    return m_CurrentVecValues;
  }

  /*
   * currentValueAsVoidP() (virtual)
   */
  const void * LE_StdVecDblVarLength::currentValueAsVoidP() const {
    return static_cast<const void *>(DATA_PTR(m_CurrentVecValues));
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
  bool LE_StdVecDblVarLength::checkRequirements(bool throwSpecificError) 
  {
    if (m_numTimesRequired) {
      if (m_numTimesRequired != m_numTimesProcessed) {
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
  void LE_StdVecDblVarLength::print_usage(int ilvl) const {
    print_indent(ilvl);
    printf("%s =  double_1 double_2                                      (OPTIONAL LINE)\n", EntryName.orig_str);
    print_indent(ilvl);
    printf("%s = (double_3 ... double_VecLength", EntryName.orig_str);
    printf(" where VecLength = arbitrary) ");
    if (MaxVal != DBL_MAX || MinVal != -DBL_MAX) {
      printf(" with limits (%g, %g)", MinVal, MaxVal);
    }
    if (DefaultVal != NO_DEFAULT_DOUBLE) {
      printf(" with default %g", DefaultVal);
    }
    if (m_numTimesRequired >= 1) {
      printf(" (REQUIRED DOUBLES = %d)", m_numTimesRequired);
    } else {
      printf(" (OPTIONAL LINE)");
    }
    printf("\n");
  }

  /**********************************************************************
   *
   * set_default()
   *
   */
  void LE_StdVecDblVarLength::set_default(double def)
  {
    DefaultVal = def;  
    if (m_numTimesProcessed == 0) {
      CurrValue = def;
      for (size_t i = 0; i < (size_t) VecLength; i++) {
        (*m_vectorDest_ptr)[i] = def;
         m_CurrentVecValues[i] = def;
      }
    }
  }

  /**********************************************************************
   *
   * set_limits()
   *
   */
  void LE_StdVecDblVarLength::set_limits(double maxV, double minV)
  {
    MaxVal = maxV;
    MinVal = minV;
  }

 
}
