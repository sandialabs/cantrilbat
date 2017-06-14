
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


  Cp_Barin Cp_Cgraph_SolA(0.026, 9.307, -0.354, -4.155);   //298 - 737
  Cp_Barin Cp_Cgraph_SolB(5.841, 0.104, -7.559, 0.0);

  Thermo_Shomate ts_Cgraph_SolA;

  ts_Cgraph_SolA.convertCpBarinToShomate(Cp_Cgraph_SolA);

  ts_Cgraph_SolA.convert_Hf_Barin_ToShomate(0.0);

  //ts_Cgraph_SolA.convert_S298_Barin_ToShomate(1.372);
  // Set it to the exact value that we are using for the element Ca entropy contribution
  ts_Cgraph_SolA.set_S298(5.740);


  double dens = 2.15;  // gm / cm3;
  double mw = 12.011; 
  double mv = mw / dens * 1.0E-3;
  
  printf("      <!-- species Cgraph_SolA   -->\n");
 
  printf("      <species name=\"Cgraph_SolA\">\n");
  printf("        <atomArray> C:1 </atomArray>\n");
  ts_Cgraph_SolA.printThermoBlock(8, 1100 , 200.);
  printf("        <standardState  model=\"constant_incompressible\">\n");
  printf("          <!-- molar volume from C(s) based on density of %g gm/cm3 --> \n", dens);
  printf("          <molarVolume units=\"m3/kmol\">  %-11.8g </molarVolume>\n", mv);
  printf("        </standardState>\n");
  printf("      </species>\n");
  printf("                \n");

  //---------------------------------------------------------------------------------------------------------------

  Thermo_Shomate ts_Cgraph_SolB;

  ts_Cgraph_SolB.convertCpBarinToShomate(Cp_Cgraph_SolB);
  ts_Cgraph_SolB.convert_PhaseChange_Barin_ToShomate(ts_Cgraph_SolA, 1100., 0.00, 0.0);

  printf("      <!-- species Cgraph_SolB  -->\n");

  printf("      <species name=\"Cgraph_SolB\">\n");
  printf("       <atomArray> C:1  </atomArray>\n");
  ts_Cgraph_SolB.printThermoBlock(8, 4073, 1100.);
  printf("       <standardState  model=\"constant_incompressible\">\n");
  printf("          <!-- molar volume from C(s) --> \n");
  printf("          <molarVolume units=\"m3/kmol\">  %-11.8g </molarVolume>\n", mv);
  printf("       </standardState>\n");
  printf("      </species>\n");
  printf("                \n");

  //---------------------------------------------------------------------------------------------------------------


  return 0;

}


