
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


  Cp_Jgmol Cp_U3O8_Sol(276.747750, 27.328834, -40.733348, 0.0);   //298 - 2888

  Thermo_Shomate ts_U3O8_Sol;

  ts_U3O8_Sol.convertCpJgmolToShomate(Cp_U3O8_Sol);

  ts_U3O8_Sol.convert_Hf_Jgmol_ToShomate(-3576818.45);
  ts_U3O8_Sol.set_S298(285.022053);


//  double dens = 3.34;  // gm / cm3;
//  double mw = 40.078 + 15.9994; 
  double mv = 24.615;  // took UO2 value
  
  printf("      <!-- species U3O8(S)  -->\n");
 
  printf("      <species name=\"U3O8(S)\">\n");
  printf("        <atomArray> U:3 O:8 </atomArray>\n");
  printf("        <!--  Melting point = 1420 k (Barin & Knacke 1973) -->");
  ts_U3O8_Sol.printThermoBlock(8, 3500. , 200.);
  printf("        <standardState  model=\"constant_incompressible\">\n");
//  printf("          <!-- molar volume from CaO(s) based on density of %g gm/cm3 --> \n", dens);
  printf("          <molarVolume units=\"m3/kmol\">  %-11.8g </molarVolume>\n", mv);
  printf("        </standardState>\n");
  printf("      </species>\n");
  printf("                \n");

  //---------------------------------------------------------------------------------------------------------------


  return 0;

}


