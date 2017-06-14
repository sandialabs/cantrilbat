
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


  Cp_Barin Cp_Ca_SolA(5.24, 3.50, 0.0, 0.0);   //298 - 737
  Cp_Barin Cp_Ca_SolB(2.59, 6.66, 0.0, 0.0);

  Thermo_Shomate ts_Ca_SolA;

  ts_Ca_SolA.convertCpBarinToShomate(Cp_Ca_SolA);

  ts_Ca_SolA.convert_Hf_Barin_ToShomate(0.0);

  ts_Ca_SolA.convert_S298_Barin_ToShomate(9.95);


  double dens = 1.55;  // gm / cm3;
  double mw = 40.078; 
  double mv = mw / dens * 1.0E-3;
  
  printf("      <!-- species Ca_SolA   -->\n");
 
  printf("      <species name=\"Ca_SolA\">\n");
  printf("        <atomArray> Ca:1 </atomArray>\n");
  ts_Ca_SolA.printThermoBlock(8, 737. , 200.);
  printf("        <standardState  model=\"constant_incompressible\">\n");
  printf("          <!-- molar volume from Ca(s) based on density of %g gm/cm3 --> \n", dens);
  printf("          <molarVolume units=\"m3/kmol\">  %-11.8g </molarVolume>\n", mv);
  printf("        </standardState>\n");
  printf("      </species>\n");
  printf("                \n");

  //---------------------------------------------------------------------------------------------------------------

  Thermo_Shomate ts_Ca_SolB;

  ts_Ca_SolB.convertCpBarinToShomate(Cp_Ca_SolB);
  ts_Ca_SolB.convert_PhaseChange_Barin_ToShomate(ts_Ca_SolA, 737, 0.06, 0.081);

  printf("      <!-- species Ca_SolB  -->\n");

  printf("      <species name=\"Ca_SolB\">\n");
  printf("       <atomArray> Ca:1  </atomArray>\n");
  ts_Ca_SolB.printThermoBlock(8, 1150., 500.);
  printf("       <standardState  model=\"constant_incompressible\">\n");
  printf("          <!-- molar volume from Ca(s) --> \n");
  printf("          <molarVolume units=\"m3/kmol\">  %-11.8g </molarVolume>\n", mv);
  printf("       </standardState>\n");
  printf("      </species>\n");
  printf("                \n");

  //---------------------------------------------------------------------------------------------------------------

  Thermo_Shomate ts_Ca_liq;
  Cp_Barin Cp_Ca_liq(7.4, 0.0, 0.0, 0.0);
   
  ts_Ca_liq.convertCpBarinToShomate(Cp_Ca_liq);
  ts_Ca_liq.convert_PhaseChange_Barin_ToShomate(ts_Ca_SolB, 1123., 2.0, 1.781);

  printf("      <!-- species Ca_liq  -->\n");

  printf("      <species name=\"Ca_liq\">\n");
  printf("       <atomArray> Ca:1  </atomArray>\n");
  ts_Ca_liq.printThermoBlock(8, 2700., 600.);
  printf("       <standardState  model=\"constant_incompressible\">\n");
  printf("          <!-- molar volume from Ca_liq --> \n");
  printf("          <molarVolume units=\"m3/kmol\">  %-11.8g </molarVolume>\n", mv);
  printf("       </standardState>\n");
  printf("      </species>\n");
  printf("                \n");


  //---------------------------------------------------------------------------------------------------------------

  Thermo_Shomate ts_Ca_gas;
  Cp_Barin Cp_Ca_gas(4.97, 0.0, 0.0, 0.0);
   
  ts_Ca_gas.convertCpBarinToShomate(Cp_Ca_gas);
  ts_Ca_gas.convert_PhaseChange_Barin_ToShomate(ts_Ca_liq, 1762., 36.5, 20.715);

  printf("      <!-- species Ca_gas  -->\n");

  printf("      <species name=\"Ca_gas\">\n");
  printf("       <atomArray> Ca:1  </atomArray>\n");
  ts_Ca_gas.printThermoBlock(8, 2700., 298.);
  printf("       <standardState  model=\"idealGas\">\n");
  printf("       </standardState>\n");
  printf("      </species>\n");
  printf("                \n");


  return 0;

}


