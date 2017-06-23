
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

     double S800 = ts_O2_low.s(800.);
     printf(" S(800) = %.11g J/gmol/K\n", S800);

     double S1000 = ts_O2_low.s(1000.);
     printf(" S(1000) = %.11g J/gmol/K\n", S1000);


     double Cp800 = ts_O2_low.cp(800.);
     printf(" Cp(800) = %.11g J/gmol/K\n", Cp800);

     double Cp1000 = ts_O2_low.cp(1000.);
     printf(" Cp(1000) = %.11g J/gmol/K\n", Cp1000);

     double coeffs_ts_O2_high[7] =  {     3.282537840E+00,   1.483087540E-03,  -7.579666690E-07,   2.094705550E-10,
             -2.167177940E-14,  -1.088457720E+03,   5.453231290E+00 };

     Thermo_NASA ts_O2_high(1000., 3500., 100000., coeffs_ts_O2_high) ;

     Regions_NASA O2_reg;

     O2_reg.addRegionPoly(&ts_O2_low);
     O2_reg.addRegionPoly(&ts_O2_high);

     double H1000_low = ts_O2_low.h(1000.0);
     double S1000_low = ts_O2_low.s(1000.0);
     double Cp1000_low = ts_O2_low.cp(1000.0);
     double H1000_high = ts_O2_high.h(1000.0);
     double S1000_high = ts_O2_high.s(1000.0);
     double Cp1000_high = ts_O2_high.cp(1000.0);

     printf(" Before: \n");
     printf(" H:    % -20.15E     % -20.15E    |  %-14.4E\n", H1000_low, H1000_high,   H1000_high  - H1000_low);
     printf(" S:    % -20.15E     % -20.15E    |  %-14.4E\n", S1000_low, S1000_high,   S1000_high  - S1000_low);
     printf("Cp:    % -20.15E     % -20.15E    |  %-14.4E\n", Cp1000_low, Cp1000_high, Cp1000_high - Cp1000_low);

     O2_reg.adjust_high();

     H1000_low = ts_O2_low.h(1000.0);
     S1000_low = ts_O2_low.s(1000.0);
     Cp1000_low = ts_O2_low.cp(1000.0);
     H1000_high = ts_O2_high.h(1000.0);
     S1000_high = ts_O2_high.s(1000.0);
     Cp1000_high = ts_O2_high.cp(1000.0);

     printf(" After: \n");
     printf(" H:    % -20.15E     % -20.15E    |  %-14.4E\n", H1000_low, H1000_high,   H1000_high  - H1000_low);
     printf(" S:    % -20.15E     % -20.15E    |  %-14.4E\n", S1000_low, S1000_high,   S1000_high  - S1000_low);
     printf("Cp:    % -20.15E     % -20.15E    |  %-14.4E\n", Cp1000_low, Cp1000_high, Cp1000_high - Cp1000_low);

     O2_reg.printThermoBlock(5, 15);

  return 0;

}


