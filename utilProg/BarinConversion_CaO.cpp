
// Example of taking a Barin-Knack Table entry and creating a shomate polynomial for use in Zuzax

/*
 *  Theory
 * -------------------------
 *
 *  Barin expresses his units for heat capacity in cal/gmol/K
 *
 *  Cp (cal/gmol/K) = A + B 1.0E-3 T + C 1.0E5 T^-2 + D 1.0E-6 T^2
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


  Cp_Barin Cp_CaO_Sol(11.86, 1.08, -1.66, 0.0);   //298 - 2888

  Thermo_Shomate ts_CaO_Sol;

  ts_CaO_Sol.convertCpBarinToShomate(Cp_CaO_Sol);

  ts_CaO_Sol.convert_Hf_Barin_ToShomate(-151.6);
  ts_CaO_Sol.convert_S298_Barin_ToShomate(9.5);


  double dens = 3.34;  // gm / cm3;
  double mw = 40.078 + 15.9994; 
  double mv = mw / dens * 1.0E-3;
  
  printf("      <!-- species CaO_Sol   -->\n");
 
  printf("      <species name=\"CaO_Sol\">\n");
  printf("        <atomArray> Ca:1 O:1 </atomArray>\n");
  printf("        <!--  Melting point = 2888.0 (Barin & Knacke 1973) -->");
  ts_CaO_Sol.printThermoBlock(8, 3500. , 200.);
  printf("        <standardState  model=\"constant_incompressible\">\n");
  printf("          <!-- molar volume from CaO(s) based on density of %g gm/cm3 --> \n", dens);
  printf("          <molarVolume units=\"m3/kmol\">  %-11.8g </molarVolume>\n", mv);
  printf("        </standardState>\n");
  printf("      </species>\n");
  printf("                \n");

  //---------------------------------------------------------------------------------------------------------------


  return 0;

}


