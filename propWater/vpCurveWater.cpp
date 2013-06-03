/* ======================================================================= */
/* $RCSfile: vpCurveWater.cpp,v $ */
/* $Author: hkmoffa $ */
/* $Date: 2012-02-23 14:34:18 -0700 (Thu, 23 Feb 2012) $ */
/* $Revision: 5 $ */
/* ======================================================================= */

#include <iostream>
 

#include "Water.h"

#include <vector>


using namespace tpx;
using std::vector;

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
  w.setStdState();

 

  vector<double> Ttable;
  vector<double> DelHLit;

  Ttable.push_back(0.1);
  DelHLit.push_back(45.054);

  Ttable.push_back(25.);
  DelHLit.push_back(43.990);

  Ttable.push_back(40.);
  DelHLit.push_back(43.35);

  Ttable.push_back(60.);
  DelHLit.push_back(42.482);

  Ttable.push_back(80.);
  DelHLit.push_back(41.585);

  Ttable.push_back(100.);
  DelHLit.push_back(40.657);

 Ttable.push_back(120.);
  DelHLit.push_back(39.684);
 Ttable.push_back(140.);
  DelHLit.push_back(38.643);
 Ttable.push_back(160.);
  DelHLit.push_back(37.518);
 Ttable.push_back(180.);
  DelHLit.push_back(36.304);
 Ttable.push_back(200.);
  DelHLit.push_back(34.962);
 Ttable.push_back(220.);
  DelHLit.push_back(33.468);
 Ttable.push_back(240.);
  DelHLit.push_back(31.809);
 Ttable.push_back(260.);
  DelHLit.push_back(29.930);
 Ttable.push_back(280.);
  DelHLit.push_back(27.795);
 Ttable.push_back(300.);
  DelHLit.push_back(25.3);
 Ttable.push_back(320.);
  DelHLit.push_back(22.297);
Ttable.push_back(340.);
  DelHLit.push_back(18.502);
Ttable.push_back(360.);
  DelHLit.push_back(12.966);
Ttable.push_back(374.);
  DelHLit.push_back(2.066);





  printf ("Temperature  Pressure  H_liq   S_liq    H_gas   S_gas   ");
  printf("DelH (DelHlit)  DelS  DelG  Fug  FugCoeff\n");
  printf("------------------------------------------------------------");
  printf("------------------------------------------------------------");
  printf("\n");



  for (int i = 0; i < Ttable.size(); i++) {
    double T = Ttable[i] + 273.15;
    double delH_book = DelHLit[i];
        
    w.Set(TP, T, 1.0E3);
    double pressure = w.Psat();
    w.Set(TP, T, pressure);
    printf(" %12g   %12g", T, pressure);

    /*
     * Liquid
     */
    w.Set(TX, T, 0.0);
    double h_liq = w.h() * w.MolWt() * 1.0E-6;
    double s_liq = w.s() * w.MolWt() * 1.0E-6;
    double g_liq = w.g() * w.MolWt() * 1.0E-6;
    printf(" %12g %12g", h_liq, s_liq);

    /*
     * gas
     */
    w.Set(TX, T, 1.0);
    double s_gas = w.s() * w.MolWt() * 1.0E-6;
    double h_gas = w.h() * w.MolWt() * 1.0E-6;
    printf(" %12g %12g", h_gas, s_gas);
    double g_gas = w.g() * w.MolWt() * 1.0E-6;
  
    double Del_H = (h_gas - h_liq);
    double Del_S = (s_gas - s_liq);
    printf(" %g (%g)", Del_H, delH_book);
    printf(" %g ", Del_S * 1.0e3);
    printf(" %g", Del_H - T * Del_S);

    w.Set(TP, T, 1.0);
    double Rgas = 8.31451E3;
    double g_low = w.g() * w.MolWt() * 1.0E-6;
    double fugac = exp((g_gas - g_low) * 1.0E6 / (Rgas * T));
    printf(" %g  %g\n", fugac, fugac / pressure);
  }
  printf("------------------------------------------------------------");
  printf("------------------------------------------------------------");
  printf("\n");
} 
