/**
 * @file LineEntry.cpp
 *
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "LE_LineEntry.h"
#include "mdp_allo.h"

#include <string>
using namespace std;

namespace BEInput {
  /*
   *  LineEntry()
   *
   *  Constructor
   *
   *  This is the constructor for the LineEntry object. We enter the
   *  number of required lines of this object and the name of the line
   *  entry here.
   */
  LineEntry::LineEntry(const char *lineName, int numRL) :
    BaseEntry(lineName, numRL)
  {
  }

  /*
   *
   * LineEntry(const LineEntry &b) :
   *
   * Copy constructor():
   */
  LineEntry::LineEntry(const LineEntry& b) :
    BaseEntry(b)
  {
  }

  /*
   *
   * LineEntry& operator=(const LineEntry& b)
   *
   * assignment operator for the LineEntry class
   */
  LineEntry& LineEntry::operator=(const LineEntry& b)
  {
    if (&b != this) {
      BaseEntry::operator=(b);
    }
    return *this;
  }

  /*
   *
   * BaseEntry* duplMyselfAsBaseEntry() : (virtual)
   *
   * Duplicate as a base class function. Note, doing it this way means
   * that derived classes only have to carry along the LineEntry*
   * function (I Think - check)
   */
  BaseEntry* LineEntry::duplMyselfAsBaseEntry() const
  {
    LineEntry* newLE = duplMyselfAsLineEntry();
    return ((BaseEntry *)newLE);
  }

  /*
   *
   * print_indent (LineIndent static function)
   */
  void LineEntry::print_indent(int ilvl)
  {
    for (int i = 0; i < ilvl; i++) {
      printf("    ");
    }
  }

  /*
   * checkRequirements():
   *
   * Check for all requirements being met at the end of input.
   *
   * @param throwSpecificError If true then you should throw a 
   *        specific error condition, if you have one. If not,
   *        an generic error condition will be thrown on return.
   */
  bool LineEntry::checkRequirements(bool throwSpecificError) {
    /*
     * The first thing we do is to check for satisfaction
     * of any Dependency Result Types that might zero
     * out the Number of times required field.
     */
    for (int i = 0; i < NumEntryDependencies; i++) {
      BI_Dependency* dep = EntryDepList[i];
      if (dep->ResultType() == BIDRT_ZERONUMTIMESREQUIRED) {
	if (dep->checkDependencySatisfied()) {
	  m_numTimesRequired = 0;
	}
      }
    }
    /*
     * The next thing we do is to check for satisfaction
     * of any Dependency Result Types that might turn on
     * the Number of times required field.
     */
    for (int i = 0; i < NumEntryDependencies; i++) {
      BI_Dependency* dep = EntryDepList[i];
      if (dep->ResultType() == BIDRT_ONENUMTR) {
	if (dep->checkDependencySatisfied()) {
	  m_numTimesRequired = 1;
	}
      }
    }

    /*
     * The next thing we do is to check for satisfaction
     * of any Dependency Result Types that might
     * be antithetical to the current card.
     */
    for (int i = 0; i < NumEntryDependencies; i++) {
      BI_Dependency* dep = EntryDepList[i];
      if (dep->ResultType() == BIDRT_ANTITHETICAL_ERROR) {
        if (m_numTimesProcessed > 0) {
	  if (dep->checkDependencySatisfied()) {
	    string ename = EntryName.orig_str;
	    string dname = dep->TargetBaseEntry()->keyname();
	    throw BI_InputError("LineEntry::checkRequirements() on  \"" + ename + "\"",
				"Mutually exclusive BaseEntry objects "
				"were invoked: \"" + dname + "\"");
	  }
	}
      }
    }

    if (m_numTimesRequired) {
      if (m_numTimesRequired != m_numTimesProcessed) {
	return false;
      }
    }
    return true;
  }

  /*
   * LineEntry Destructor (virtual function):
   *
   *  note this is a virtual function
   */
  LineEntry::~LineEntry()
  {
#ifdef DEBUG_DESTRUCTOR
    printf("~LineEntry called for %s\n", LineName.orig_str);
#endif
  }

  /*
   *   This routine is called to actually do something with the LineArgTok
   *   argument.
   *   What's here is just a stub routine. The main work is done by
   *   the member function process_LineEntry of derived classes of this 
   *   class.
   *   However, here we increment the NumProcessedLines variable,
   * , i.e., everything that is common to  all LineEntry objects.
   *   This may be called be every LineEntry child class before it
   *   does its own specific processing.
   */
  void
  LineEntry::process_LineEntry(const TK_TOKEN *lineArgTok)
  {
    for (int i = 0; i < NumEntryDependencies; i++) {
      BI_Dependency *idep = EntryDepList[i];
      if (idep->ResultType() ==  BIDRT_PT_ERROR) {
	if (!idep->checkDependencySatisfied()) {  
	  string ename = EntryName.orig_str;
	  string dname = idep->TargetBaseEntry()->keyname();
	  throw BI_InputError("LineEntry::process_LineEntry on \"" + ename + "\"",
			      "Dependency not satisfied wrt keyname \"" + dname + "\"");
	}
      }

      if (idep->ResultType() == BIDRT_ANTITHETICAL_ERROR) {
	if (idep->checkDependencySatisfied()) {
	  string ename = EntryName.orig_str;
	  string dname = idep->TargetBaseEntry()->keyname();
	  throw BI_InputError("LineEntry::process_LineEntry on \"" + ename + "\"",
			      "Mutual Exclusive Dependency is satisfied wrt keyname \"" + dname + "\"");
	}
      }
    }
    m_numTimesProcessed++;
  }

  /*
   * LineEntry::print_ProcessedLine() (virtual function):
   *
   *   This routine will print out a processed line
   *
   *  The default behavior is to print the original line with a "=>"
   *  prefix to indicate that action has been taken on it.
   */
  void
  LineEntry::print_ProcessedLine(const TK_TOKEN *lineArgTok) const {
    printf(" => %s", EntryName.orig_str);
    if (lineArgTok) {
      if (lineArgTok->orig_str) {
	printf(" = %s\n", lineArgTok->orig_str);
      } else {
	printf("\n");
      }
    } else {
      printf("\n");
    }
  }

}
