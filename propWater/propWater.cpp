/* ======================================================================= */
/* $RCSfile: propWater.cpp,v $ */
/* $Author: hkmoffa $ */
/* $Date: 2012-02-23 14:34:18 -0700 (Thu, 23 Feb 2012) $ */
/* $Revision: 5 $ */
/* ======================================================================= */

#include <iostream>

#include "propWater.h" 

#include "Water.h"

using namespace tpx;


void Ephase(water &w) {
  if (w.TwoPhase()) {
    printf("we are in two phase zone, T = %g, P = %g, x = %g\n",
	   w.Temp(), w.P(), w.x());
  } else {
    printf("we are not in two phase zone, T = %g, P = %g, x = %g\n",
	   w.Temp(), w.P(), w.x());
  }
}

int main()
{
  

  water w;
  //w.setStdState();

  //w.set_TPp(300., 1.013E5);
  w.Set(TP, 500., 1.013E7);
  double result = w.Ps();
  printf("result = %g\n", result);

  printf("Critical Temperature = %g Kelvin -> CRC value = %g \n",
	 w.Tcrit(), 373.99 + 273.15);
  printf("Critical Pressure    = %g Pa -> CRC value = %g\n",  
	 w.Pcrit(), 22.064 * 1.0E6);
  printf("Critical Volume = %g m^3/kg -> CRC value = %g\n",
	 w.Vcrit(), 1.0/ (0.322) / 1.0E6 * 1.0E3);
  printf("Molecular Weight = %g gm/mol > CRC value = %g\n",
	 w.MolWt(), 18.01528);
  printf("Tmin = %g Kelvin\n", w.Tmin());
  printf("Tmax = %g Kelvin\n", w.Tmax());

  w.set_TPp(273.16, 1.013E5);
  printf("PSat(0 C) = %g Pa -> CRC value = %g\n",
	 w.Psat(), 0.6113 * 1.0E3);

  w.set_TPp(273.15 + 20, 1.013E5);
  printf("PSat(20C) = %g Pa -> CRC value = %g\n",
	 w.Psat(), 2.3388 * 1.0E3);

  w.set_TPp(273.15 + 25, 1.013E5);
  printf("PSat(25C) = %g Pa\n", w.Psat());

  w.set_TPp(273.15 + 40, 1.013E5);
  printf("PSat(40C) = %g Pa -> CRC value = %g\n",
	 w.Psat(), 7.3814 * 1.0E3);

  w.set_TPp(273.15 + 60, 1.013E5);
  printf("PSat(60C) = %g Pa -> CRC value = %g\n",
	 w.Psat(), 19.932 * 1.0E3);

  w.set_TPp(373.15, 1.013E5);
  printf("PSat(100C) = %g Pa -> CRC value = %g\n",
	 w.Psat(), 101.325 * 1.0E3);

  printf("Name = %s\n", w.name());
  printf("Formula = %s\n", w.formula());
  printf("------------------------------------\n");
  printf("Triple Point:\n");
  w.Set(TP, 273.17, 1.0E5);
  Ephase(w);
  printf("v = %g, rho = %g\n", w.v(), 1.0/w.v());
  double ptrip = w.Psat();
  printf("ptrip = %g\n", ptrip);
  /* Keep this line -> This is the triple point */
  w.Set(TV, 273.16, 1.0/999.771605);

  Ephase(w);
  printf("Liquid Entropy and Enthalpy are set to zero at triple point:\n");
  printf("s = %g J/(kg K)\n", w.s());
  printf("h = %g J/kg\n", w.h());
  printf("T = %g, P = %g\n", w.Temp(), w.P());
  printf("g_liq = %g J/kg", w.g());
  w.Set(TX, 273.16, 1.0);
  Ephase(w);
  printf("s = %g J/(kg K)\n", w.s());
  printf("h = %g J/kg\n", w.h());
  printf("T = %g, P = %g\n", w.Temp(), w.P());
  printf("g_gas = %g J/kg\n", w.g());
  printf("Psat = %g\n", w.Psat());
  printf("------------------------------------\n");



  w.Set(TP, 298.15, 1.00E5);
  Ephase(w);
  printf("s = %g J/(kg K)\n", w.s());
  printf("h = %g J/kg\n", w.h());
  printf("s = %g J/(K gmol)\n", w.s() * w.MolWt() * 1.0E-3);
  printf("h = %g kJ/gmol\n", w.h() * w.MolWt() * 1.0E-6);

  double h0 = -285.83 * 1.0E6 / w.MolWt();
  double s0 = 69.95 * 1.0E3 / w.MolWt();
  w.setStdState(h0, s0, 298.15, 1.0E5) ;

  printf("h = %g kJ/gmol\n", w.h() * w.MolWt() * 1.0E-6);
  printf("s = %g J/(K gmol)\n", w.s() * w.MolWt() * 1.0E-3);

  w.Set(TP, 300., 1.00E5);
  Ephase(w);
  printf("h - h0 = %g kJ/gmol\n", (w.h() - h0) * w.MolWt() * 1.0E-6);
  printf("s = %g J/(K gmol)\n", w.s() * w.MolWt() * 1.0E-3);

  w.Set(TP, 373.21, 1.11E5);
  Ephase(w);
  printf("h - ho = %g kJ/gmol\n", (w.h()-h0) * w.MolWt() * 1.0E-6);
  printf("s = %g J/(K gmol)\n", w.s() * w.MolWt() * 1.0E-3);

  w.Set(TP, 400., 1.00E5);
  Ephase(w);
  printf("h = %g kJ/gmol\n", (w.h()) * w.MolWt() * 1.0E-6);
  printf("s = %g J/(K gmol)\n", w.s() * w.MolWt() * 1.0E-3);

  w.Set(TP, 298.15, 1.00);
  Ephase(w);
  printf("h = %g kJ/gmol\n", (w.h()) * w.MolWt() * 1.0E-6);
  printf("s = %g J/(K gmol)\n", w.s() * w.MolWt() * 1.0E-3);
  double stmp =  w.s() * w.MolWt() * 1.0E-3;
  stmp -= 8.31451 * log(1.0E5);
  printf("stmp - corrected = %g J/(K gmol)\n",stmp);
} 
