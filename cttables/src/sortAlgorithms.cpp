/*
 * @file sortAlgorithms.h
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

#include "sortAlgorithms.h"

/**************************************************************/
/**
 * sort_dbl_1():
 *   This algorithm sorts a vector of doubles in increasing
 *   value order using the heapsort algorithm.
 *    (Numerical Recipes in C routine).
 */
void sort_dbl_1(double * const x, const int n) {
    if (!x) return;
    if (n <= 1) return;
    double rra;
    int ll = n/2;
    int iret = n - 1;
    while (1 > 0) {
      if (ll > 0) {
        ll--;
        rra = x[ll];
      } else {
        rra = x[iret];
        x[iret] = x[0];
        iret--;
        if (iret == 0) {
          x[0] = rra;
          return;
        }
      }     
      int i = ll;
      int j = ll + ll + 1;
      while (j <= iret) {
        if (j < iret) {
          if (x[j] < x[j+1])
	      j++;
        }
        if (rra < x[j]) {
          x[i] = x[j];
          i = j;
          j = j + j + 1;
        } else {
          j = iret + 1;
        }
      }
      x[i] = rra;
    }
}
/*****************************************************/
