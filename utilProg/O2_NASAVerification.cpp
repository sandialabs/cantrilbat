
// Example of using O2 NASA polynomial

/*
 *  Theory
 * -------------------------
 *
 *
 *
 *
 *
 */
#include <stdio.h>
#include <cmath>
#include <cstdlib>

#include "thermoUtilProg.h"

// -------------------------------------------------------
int main () {


  double coeffs_ts_O2_low[7] = {   3.782456360E+00,  -2.996734160E-03,   9.847302010E-06,  -9.681295090E-09,
             3.243728370E-12,  -1.063943560E+03,   3.657675730E+00  };

  Thermo_NASA ts_O2_low(250., 1000., 100000., coeffs_ts_O2_low) ;

   double hkJ = ts_O2_low.h(298.15);
   printf("hkJ = %g kJ/gmol\n", hkJ);

   ts_O2_low.printThermoBlock(5);

   ts_O2_low.adjustH(298.15, 1.0);

    printf(" Hf should be 1 -> %20.11E \n", ts_O2_low.h(298.15)); 

   ts_O2_low.adjustH(298.15, 0.0);

    ts_O2_low.printThermoBlock(5);


  return 0;

}


