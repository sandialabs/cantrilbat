/**
 * @file subtractRD.cpp
 *
 * $Author: hkmoffa $
 * $Revision: 5 $
 * $Date: 2012-02-23 14:34:18 -0700 (Thu, 23 Feb 2012) $
 */
/*
 * Copywrite 2005 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include <math.h>
#define MIN(x,y)     (( (x) < (y) ) ? (x) : (y))
namespace mdpUtil {
/*************************************************************************
 *
 * subtractRD():
 *   This routine subtracts 2 numbers. If the difference is less
 *   than 1.0E-14 times the magnitude of the smallest number,
 *   then diff returns an exact zero. 
 *   It also returns an exact zero if the difference is less than
 *   1.0E-300.
 *
 *   returns:  a - b
 *
 *   This routine is used in numerical differencing schemes in order
 *   to avoid roundoff errors resulting in creating Jacobian terms.
 *   Note: This is a slow routine. However, jacobian errors may cause
 *         loss of convergence. Therefore, in practice this routine
 *         has proved cost-effective.
 */
double subtractRD(double a, double b) {
    double diff = a - b;
    double d = MIN(fabs(a), fabs(b));
    d *= 1.0E-14;
    double ad = fabs(diff);
    if (ad < 1.0E-300) {
      diff = 0.0;
    }
    if (ad < d) {
      diff = 0.0;
    }
    return diff;
}

}

/*
 * This is the fortran binding version of this routine
 */
extern "C" double subtractrd_(double *aptr, double *bptr) {
    double diff = *aptr - *bptr;
    double d = MIN(fabs(*aptr), fabs(*bptr));
    d *= 1.0E-14;
    double ad = fabs(diff);
    if (ad < 1.0E-300) {
      diff = 0.0;
    }
    if (ad < d) {
      diff = 0.0;
    }
    return diff;
}
