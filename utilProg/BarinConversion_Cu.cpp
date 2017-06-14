
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


  Cp_Barin Cp_Cu_SolA(5.41, 1.50, 0.0, 0.0);   //298 - 1357
  Cp_Barin Cp_Cu_liq(7.5, 0.0, 0.0, 0.0);

  Thermo_Shomate ts_Cu_SolA;

  ts_Cu_SolA.convertCpBarinToShomate(Cp_Cu_SolA);

  ts_Cu_SolA.convert_Hf_Barin_ToShomate(0.0);

  //ts_Ca_SolA.convert_S298_Barin_ToShomate(9.95);
  // Set it to the exact value that we are using for the element Ca entropy contribution
  //ts_Cu_SolA.convert_S298_Barin_ToShomate(7.97);
  ts_Cu_SolA.set_S298(33.164);  // Set to the value given in elements.xml

  double dens = 8.96;  // gm / cm3;
  double mw = 63.546; 
  double mv = mw / dens * 1.0E-3;
  
  printf("      <!-- species Cu_Sol   -->\n");
 
  printf("      <species name=\"Cu_Sol\">\n");
  printf("        <atomArray> Cu:1 </atomArray>\n");
  ts_Cu_SolA.printThermoBlock(8, 2000. , 200.);
  printf("        <standardState  model=\"constant_incompressible\">\n");
  printf("          <!-- molar volume from Cu(s) based on density of %g gm/cm3 --> \n", dens);
  printf("          <molarVolume units=\"m3/kmol\">  %-11.8g </molarVolume>\n", mv);
  printf("        </standardState>\n");
  printf("      </species>\n");
  printf("                \n");

  //---------------------------------------------------------------------------------------------------------------

  Thermo_Shomate ts_Cu_liq;

  ts_Cu_liq.convertCpBarinToShomate(Cp_Cu_liq);
  ts_Cu_liq.convert_PhaseChange_Barin_ToShomate(ts_Cu_SolA, 1357., 3.12, 2.299);
  dens = 8.02;  // gm / cm3;
  mv = mw / dens * 1.0E-3;

  printf("      <!-- species Cu_liq  -->\n");

  printf("      <species name=\"Cu_liq\">\n");
  printf("       <atomArray> Cu:1  </atomArray>\n");
  ts_Cu_liq.printThermoBlock(8, 3000., 500.);
  printf("       <standardState  model=\"constant_incompressible\">\n");
  printf("          <!-- molar volume from Cu(s) --> \n");
  printf("          <molarVolume units=\"m3/kmol\">  %-11.8g </molarVolume>\n", mv);
  printf("       </standardState>\n");
  printf("      </species>\n");
  printf("                \n");

  //---------------------------------------------------------------------------------------------------------------

  Thermo_Shomate ts_Cu_gas;
  Cp_Barin Cp_Cu_gas(5.369, -0.720, -0.205, 0.307);
   
  ts_Cu_gas.convertCpBarinToShomate(Cp_Cu_gas);
  ts_Cu_gas.convert_PhaseChange_Barin_ToShomate(ts_Cu_liq, 2846., 72.6, 25.509);

  printf("      <!-- species Cu_gas  -->\n");
  printf("      <species name=\"Cu_gas\">\n");
  printf("       <atomArray> Cu:1  </atomArray>\n");
  ts_Cu_gas.printThermoBlock(8, 3500., 200.);
  printf("       <standardState  model=\"IdealGas\">\n");
  printf("       </standardState>\n");
  printf("      </species>\n");
  printf("                \n");



  return 0;

}


