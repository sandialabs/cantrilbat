/**
 *  @file atofChecker
 *
 *  $Id: atofChecker.cpp 497 2013-01-07 21:17:04Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "zuzax/base/stringUtils.h"
#include "zuzax/base/ctexceptions.h"
#include "zuzax/base/global.h"

#include <iostream>
#include <new>
#include <string>
#include <vector>
#include <typeinfo>

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

using namespace Zuzax;
using namespace std;


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

int main(int argc, char** argv) { 

    char hhhs[80];
    double dval;

    /*
     * General Catch block to trap Zuzax Errors and print them
     */
    try {
      strcpy(hhhs, "34.325d10");
      string hhh(hhhs);
      dval = fpValueCheck(hhh);
      printf("dval = %g\n", dval);
      //strcpy(hhh, "-34.325d-10  3");
      //dval = atofCheck(hhh);
      //printf("dval = %g\n", dval);
      strcpy(hhhs, "-34.2d 01");
      hhh = hhhs;
       
      dval = fpValueCheck(hhh);
      printf("dval = %g\n", dval);

    } catch (ZuzaxError) {
      showErrors(cout);
      return 0;
    }
    return 0;
}
/*******************************************************************/
