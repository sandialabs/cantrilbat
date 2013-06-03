/**
 * @file tok_readFile.cpp
 *
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

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <climits>
#include <cctype>

#include "tok_input_util.h"
#include "tok_readFile.h"
#include <vector>
using std::vector;

using namespace TKInput;

namespace TKInput {
  /**
   * read a list of doubles from a file.
   */
  int tok_read2dbl(const char * const file_name, 
		   vector<double>& a1, vector<double>& a2) {
    TOKEN lineTok;
    int nvsize = 0;
    a1.resize(0);
    a2.resize(0);
    char input[1024];
    int nchar;
    int iline = 0;
    int numExpected = -1;
    double tmp1, tmp2;
    BOOLEAN error = 0;

    static FILE *ifp;
    ifp = fopen(file_name, "r");
    if (!ifp) {
      return -1;
    }

    do {
      /*
       * go get the next line
       */
      nchar = read_line(ifp, input, 0);
      iline++;
      if (nchar > 0) {
	fillTokStruct(&lineTok, input);
	if (lineTok.ntokes != 2) {
	  if (iline == 1 && lineTok.ntokes == 1) {
	    numExpected = str_to_int(lineTok.tok_ptrV[0], INT_MAX, 1, 
				     -1, &error);
	    continue;
	  } else { 
	    printf("ERROR iline = %d, Number tokens != 2: %s\n", iline,
		   lineTok.orig_str);
	    return(-1);
	  }
	}

	tmp1 = str_to_double(lineTok.tok_ptrV[0], DBL_MAX, -DBL_MAX, 
			     NO_DEFAULT_DOUBLE, &error);
	if (error) {
	  printf("ERROR iline = %d, reading first token = %s\n",
		 iline, lineTok.tok_ptrV[0]);
	  return(-1);
	}
	a1.push_back(tmp1);
	tmp2 = str_to_double(lineTok.tok_ptrV[1], DBL_MAX, -DBL_MAX, 
			     NO_DEFAULT_DOUBLE, &error);
	if (error) {
	  printf("ERROR iline = %d, reading second token = %s\n",
		 iline, lineTok.tok_ptrV[1]);
	  return(-1);
	}
	a2.push_back(tmp2);
	nvsize++;
      }
    } while (nchar >= 0);
    fclose(ifp);
    return nvsize;
  }
}
