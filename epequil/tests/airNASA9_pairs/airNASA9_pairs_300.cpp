/*
 *  $Author: hkmoffa $
 *  $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
 *  $Revision: 508 $
 *
 *  Copyright 2002 California Institute of Technology
 *
 */

#include "zuzax/IdealGasMix.h"
#include "zuzax/equilibrium.h"

#include <cstdio>

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

int main(int argc, char **argv) {
  try {
    double vol, inte;
    IdealGasMix g("airNASA9.xml", "airNASA9");
    double pres = 1.0E5;
    g.setState_TPX(300.0, pres, "N2:0.7, O2:0.3");
    // equilibrate(g, "TP", -1);
    //cout << g;

    double dens = g.density();
    vol = 1.0/dens; 

#ifdef DEBUG_CHEMEQUIL
    Zuzax::ChemEquil_print_lvl = 3;
#endif

    inte = g.intEnergy_mass();
    printf(" inte = %g\n", inte);
    inte -= 3.0E3;
    printf("attempted equil at (U,V,dens) = %10.5g, %10.5g, %10.5g\n", inte, vol, 1.0/vol);
    g.setState_UV(inte, vol);
    equilibrate(g, "UV", -1);
    cout << g;



    double enth = g.enthalpy_mass();
    printf(" enth = %g\n", enth);
    enth += 1.0E5;
    printf("attempted equil at (H,P) = %10.5g, %10.5g\n", enth, pres);
    g.setState_HP(enth, pres);
    equilibrate(g, "HP", -1);
    cout << g;

    double entrop = g.entropy_mass();
    printf(" entropy = %g\n", entrop);
    entrop += 1.0E2;
    printf("attempted equil at (S,P) = %10.5g, %10.5g\n", entrop, pres);
    g.setState_SP(entrop, pres);
    equilibrate(g, "SP", -1);
    cout << g;

    dens = g.density();
    printf(" dens = %g\n", dens);
    dens *= 0.9;
    vol = 1.0/dens; 
    printf("attempted equil at (S,V,dens) = %10.5g, %10.5g, %10.5g\n", entrop, vol, 1.0/vol);
    g.setState_SV(entrop, vol);
    equilibrate(g, "SV", -1);
    cout << g;

    double temp = 300.;
    printf("attempted equil at (T,V, dens) = %10.5g, %10.5g, %10.5g\n", temp, vol, 1.0/vol);
    g.setTemperature(temp);
    equilibrate(g, "TV", -1);
    cout << g;

    inte = g.intEnergy_mass();
    printf(" inte = %g\n", inte);
    inte -= 1.0E5;
    printf("attempted equil at (U,V,dens) = %10.5g, %10.5g, %10.5g\n", inte, vol, 1.0/vol);
    g.setState_UV(inte, vol);
    equilibrate(g, "UV", -1);
    cout << g;



    return 0;
  }
  catch (ZuzaxError) {
    showErrors(cerr);
    cerr << "program terminating." << endl;
    return -1;
  }
}
