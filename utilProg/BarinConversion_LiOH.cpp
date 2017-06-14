
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

//==================================================================================================================================
// -------------------------------------------------------
int main () {


  Cp_Barin Cp_LiOH(11.99, 8.24, -2.27, 0.0);

  Thermo_Shomate ts_LiOH;

  ts_LiOH.convertCpBarinToShomate(Cp_LiOH);

  ts_LiOH.convert_Hf_Barin_ToShomate(-115.84);

  ts_LiOH.convert_S298_Barin_ToShomate(10.2);


  double dens = 1.46;  // gm / cm3;
  double mw = 2.0 * 6.941 + 32.066; 
  double mv = mw / dens * 1.0E-3;
  
  printf("      <!-- species LiOH(S)   -->\n");
 
  printf("      <species name=\"LiOH(S)\">\n");
  printf("        <atomArray> Li:1 O:1 H:1 </atomArray>\n");
  ts_LiOH.printThermoBlock(8, 1700., 250.);
  printf("        <standardState  model=\"constant_incompressible\">\n");
  printf("          <!-- molar volume from Li2O(s) based on density of %g gm/cm3 --> \n", dens);
  printf("          <molarVolume units=\"m3/kmol\">  %-11.8g </molarVolume>\n", mv);
  printf("        </standardState>\n");
  printf("      </species>\n");
  printf("                \n");

  //---------------------------------------------------------------------------------------------------------------
  // do LiOH(l)

  Thermo_Shomate LiOH_l;

  Cp_Barin Cp_LiOH_l(20.74, 0.0, 0.0, 0.0);
  Thermo_Shomate ts_LiOH_l;
  ts_LiOH_l.convertCpBarinToShomate(Cp_LiOH_l);
  

  ts_LiOH_l.convert_DeltaH_Barin_ToShomate(ts_LiOH, 744.3, 5.01);
  ts_LiOH_l.convert_DeltaS_Barin_ToShomate(ts_LiOH, 744.3, 6.731);
  ts_LiOH_l.convert_PhaseChange_Barin_ToShomate(ts_LiOH, 744.3, 5.01, 6.731);

  printf("      <!-- species LiOH(L)   -->\n");
  printf("      <species name=\"LiOH(L)\">\n");
  printf("       <atomArray> Li:1 O:1 H:1 </atomArray>\n");
  ts_LiOH_l.printThermoBlock(8, 2500., 1500.);
  printf("       <standardState  model=\"constant_incompressible\">\n");
  printf("          <!-- molar volume from Li2O(s) --> \n");
  printf("          <molarVolume units=\"m3/kmol\">  %-11.8g </molarVolume>\n", mv);
  printf("       </standardState>\n");
  printf("      </species>\n");
  printf("                \n");

  return 0;
}
