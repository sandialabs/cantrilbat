/* ======================================================================= */
/* $RCSfile: FeS2_OpenCircuitVoltage.cpp,v $ */
/* $Author: hkmoffa $ */
/* $Date: 2013-01-07 15:05:01 -0700 (Mon, 07 Jan 2013) $ */
/* $Revision: 498 $ */
/* ======================================================================= */



#include "cantera/thermo.h"
#include "cantera/thermo/MargulesVPSSTP.h"
#include "cantera/thermo/FixedChemPotSSTP.h"

#include "cantera/equilibrium.h"

#include "cantera/thermo.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/base/logger.h"
#include <cmath>


#include <cstdio>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace Cantera;

int CHECK_DEBUG_MODE = 0;

void printUsage() {
    cout << "usage: fixedElemPot [-h] [-d #] " <<  endl;
    cout <<"                -> Everything is hardwired" << endl;
}


int main(int argc, char **argv)
{
  VCSnonideal::vcs_timing_print_lvl = 0;
  int printLvl = 2;
  int retn = 0;
  
  try {

    /*
     * Process the command line arguments
     */
    if (argc > 1) {
      string tok;
      for (int j = 1; j < argc; j++) {
	tok = string(argv[j]);
	if (tok[0] == '-') {
	  int nopt = static_cast<int>(tok.size());
	  for (int n = 1; n < nopt; n++) {
	    if (!strcmp(tok.c_str() + 1, "help_cmdfile")) {
	    } else if (tok[n] == 'h') {
	      printUsage();
	      exit(1);
	    } else if (tok[n] == 'd') {
	      int lvl = 2;
	      if (j < (argc - 1)) {
		string tokla = string(argv[j+1]);
		if (strlen(tokla.c_str()) > 0) {
		  lvl = atoi(tokla.c_str());
		  n = nopt - 1;
		  j += 1;
		  if (lvl >= 0 && lvl <= 1000) {
		    printLvl = lvl;
		  }
		}
	      }
	    } else {
	      printUsage();
	      exit(1);
	    }
	  }
	} else {
	  printUsage();
	  exit(1);
	}
      }
    }

    ofstream outfile;
    outfile.open("ocvOutput.csv",ios::out);
    outfile << "Open Circuit Voltage Curve" << endl;
    outfile << "\n";


    int tab = 15;
    outfile << setw(tab+1) << "x_Fe," << setw(tab+1) << "x_S," << setw(tab+1) << "x_Li," 
            << setw(tab+1) << "Ah/kg FeS2," << setw(tab+1) << "Volts," << setw(tab+1) << "DeltaGibbs\n";


    double um[20];
    
    double Temp = 723.15;

    FixedChemPotSSTP *LiFixed = new FixedChemPotSSTP("Li", -2.3E7);

    LiFixed->setState_TP(Temp, OneAtm);

    LiFixed->getChemPotentials(um);
    printf(" chem pot = %g\n", um[0]);

    ThermoPhase *FeS2 = newPhase("FeS2.xml","FeS2(S)");
    ThermoPhase *Li3Fe2S4 = newPhase("Li3Fe2S4.xml","Li3Fe2S4(S)");
    ThermoPhase *Li_liq = newPhase("Li_Liquid.xml","Li(L)");

    MargulesVPSSTP *Li2FeS2 = new MargulesVPSSTP("Li2_xFe1_xS2.xml", "Li[2+x]Fe[1-x]S2(S)");
 
    MargulesVPSSTP *FeS = new MargulesVPSSTP("Fe1_xS.xml", "Fe[1-x]S(S)");

    FeS2->setState_TP(Temp, OneAtm);
    Li3Fe2S4->setState_TP(Temp, OneAtm);
    Li_liq->setState_TP(Temp, OneAtm);
    Li2FeS2->setState_TP(Temp, OneAtm);
    FeS->setState_TP(Temp, OneAtm);


    FeS2->getChemPotentials(um);
    double um_FeS2 = um[0];

    Li3Fe2S4->getChemPotentials(um);
    double um_Li3Fe2S4 = um[0];

    Li_liq->getChemPotentials(um);
    double um_Li_liq = um[0];

    double deltaG = 1.0/3.0 * um_Li3Fe2S4 - 2.0/3.0 * um_FeS2 - um_Li_liq;

    printf("delta G = %g J/kmol \n",  deltaG);
    double voltsg = - deltaG / Faraday;
    printf("               in volts = %g\n", voltsg);
    

    printf("um_li = %g\n", um_Li_liq);


    double X_feS[3];
    double X_Li2FeS2[3];
    int i_FeS = FeS->speciesIndex("FeS(B)");
    int i_VaS = FeS->speciesIndex("VaS(B)");

    double z = 0.1;
    double y = (0.2 / (1.0 - 2 * z));
  

    X_feS[i_VaS] = z;
    X_feS[i_FeS] = 1 - z;
    FeS->setState_TPX(Temp, OneAtm, X_feS);
    
    
    int i_Li2FeS2 = Li2FeS2->speciesIndex("Li2FeS2(A)");
    int i_Li73 = Li2FeS2->speciesIndex("Li[7/3]Fe[2/3]S2(A)");

    X_Li2FeS2[i_Li73] = 0.3;
    X_Li2FeS2[i_Li2FeS2] = 1.0 - 0.3;
    Li2FeS2->setState_TPX(Temp, OneAtm, X_Li2FeS2);



    ThermoPhase *Fe = newPhase("Fe_Solid.xml","Fe(S)");
    ThermoPhase *Li2S = newPhase("Li2S.xml","Li2S(S)");
    Fe->setState_TP(Temp, OneAtm);
    Li2S->setState_TP(Temp, OneAtm);

    Fe->getChemPotentials(um);
    double um_Fe = um[0];


    Li2S->getChemPotentials(um);
    double um_Li2S = um[0];

    Li2FeS2->getStandardChemPotentials(um);
    double um_Li2FeS2 = um[i_Li2FeS2];

    deltaG = um_Li2S + 0.5 * um_Fe - 0.5 * um_Li2FeS2 - um_Li_liq;

    printf("delta G = %g J/kmol \n",  deltaG);
    voltsg = - deltaG / Faraday;
    printf("               in volts = %g\n", voltsg);


    //  double  volts = 1.700;  // Linear Region
    //  double  volts = 1.68;  // Linear Region
    //  double  volts = 1.66;  // Linear Region
    //   double  volts = 1.655;  // Linear Region
    // double  volts = 1.654;  //  Linear Region
    //  double  volts = 1.653;  //  has some Fe in it
    // double  volts = 1.651;  // has some Fe in it
    //     double  volts = 1.650;  // has some Fe in it
    double  volts = 1.635;     // has some Fe in it // test suite
    // double volts = 1.63;  // has some Fe in it
    // volts = 1.62953         // nominal boundary
    //double volts = 1.629; 
    // double  volts = 1.625;
    // double  volts = 1.61;


    double dg_corr =  - volts * Faraday; 
    printf("dg_corr = %g\n", dg_corr);

    double um_li_chempot = um_Li_liq + dg_corr;

    printf("um_li_chempot = %g\n", um_li_chempot);
    LiFixed->setChemicalPotential(um_li_chempot);
    
    //int numVoltages = 39;
    //double voltages[] = {2.0, 1.99, 1.98, 1.97, 1.96, 1.95, 1.94, 1.93, 1.92, 1.91, 1.9000,    1.8900,    1.8800,    1.8700,    1.8600,    1.8500,    1.8400,    1.8300,    1.8200,    1.8100,    1.8000,    1.7900,    1.7800,    1.7700,    1.7600,    1.7500,    1.7400,    1.7300,    1.7200,    1.7100,    1.7000,    1.6900,    1.6800,    1.6700,    1.6600,    1.6500,    1.6400,    1.6300,    1.6200};

    // do it this way...wrong transition point
    //int numVoltages = 3;
    //double voltages[] = {1.97, 1.96, 1.95};

    // do it this way...right transition point
    //int numVoltages = 2;
    //double voltages[] = {1.96, 1.95};

    double dV = 0.01;
    double voltageMin = 1.5;
    double voltageMax = 2.07;
    int numVoltages;
    numVoltages = 1 + (voltageMax-voltageMin)/dV;

    double elemN[3];
    for (int i = 0; i < numVoltages; i++) {

      //volts = voltages[i];
      volts = voltageMin + dV*i;

      Cantera::MultiPhase mmm;

      mmm.addPhase(FeS2, 1.);
      mmm.addPhase(Li3Fe2S4, 1.);
      mmm.addPhase(LiFixed, 18.);
      mmm.addPhase(Li2FeS2, 1.);
      mmm.addPhase(FeS, y);

      mmm.addPhase(Fe, 1.);
      mmm.addPhase(Li2S, 2.);
      mmm.init();

      int iLiFixed = 2;
      int iFe = mmm.elementIndex("Fe");
      int iS = mmm.elementIndex("S");
      int iLi = mmm.elementIndex("Li");
      int   solver = 2;
      int estimateEquil = 0;

      mmm.getElemAbundances(elemN);
      //double ratio = elemN[iS]/elemN[iFe];
      if (elemN[iS] > 2*elemN[iFe]) {
	mmm.setPhaseMoles(5,1.+elemN[iS]/2-elemN[iFe]);
      } else {
	mmm.setPhaseMoles(6,2.+2*elemN[iFe]-elemN[iS]);
      }
      
      dg_corr =  - volts * Faraday; 
      um_li_chempot = um_Li_liq + dg_corr;
      LiFixed->setChemicalPotential(um_li_chempot);
      
      vcs_equilibrate(mmm, "TP", estimateEquil, printLvl, solver);
      //cout << mmm << endl;
      
      mmm.getElemAbundances(elemN);
      elemN[iLi] -= mmm.phaseMoles(iLiFixed);
      double mass = elemN[iFe] * 55.845 + elemN[iS] * 32.065;
      double gibbs = 0;
      for (size_t j = 0; j <  mmm.nPhases(); j++) {
	gibbs += mmm.phaseMoles(j) * mmm.phase(j).gibbs_mole();
      }
      gibbs -= mmm.phaseMoles(iLiFixed) * mmm.phase(iLiFixed).gibbs_mole();
      
      outfile << setw(tab) << elemN[iFe] << "," << setw(tab) << elemN[iS] << "," 
              << setw(tab) << elemN[iLi] << "," << setw(tab) << Faraday * elemN[iLi] / mass  
              << "," << setw(tab) << volts << "," << setw(tab) << gibbs << endl;
      
    }

    outfile.close();
   
    Cantera::appdelete();

    return retn;
  

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
    
}
