/*
 *  $Author: hkmoffa $
 *  $Date: 2010/10/20 18:22:53 $
 *  $Revision: 1.1 $
 *
 *  Copyright 2002 California Institute of Technology
 *
 */

#include <stdio.h>
#include "IdealGasMix.h"
#include "equilibrium.h"
#include "PureFluid.h"

#include "thermo.h"
#include "zuzax/equil/vcs_internal.h"
#include "kernel/logger.h"
#include "kernel/HMWSoln.h"
#include "kernel/IdealSolidSolnPhase.h"


using namespace Zuzax;
using namespace std;

void printUsage() {
  cout << "usage: TBP_H2O_titr [-h] [-help_cmdfile] [-d #] "
       <<  endl;
  cout << "    -h           help" << endl;
  cout << "    -d           #   : level of debug printing" << endl;
  cout << endl;
}

int main(int argc, char **argv) {
  try {

    bool printInputFormat = false; // print cmdfile.txt format 
    bool printedUsage = false;


    int solver = 2;
 //   double pres = OneAtm;
    double pres = 1.0E+05;  // Pascal units
    
#ifdef USE_VCSNONIDEAL
    solver = 2;
    VCSnonideal::vcs_timing_print_lvl = 0;
#endif
    int printLvl = 2;
    int estimateEquil = 0;

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
	      printInputFormat = true;
	    } else if (tok[n] == 'h') {
	      printUsage();
	      printedUsage = true;
	      exit(1);
	    } else if (tok[n] == 'd') {
	      printLvl = 2;
	      int lvl = 2;
	      if (j < (argc - 1)) {
		string tokla = string(argv[j+1]);
		if (strlen(tokla.c_str()) > 0) {
		  lvl = atoi(tokla.c_str());
		  n = nopt - 1;
		  j += 1;
		  if (lvl >= 0 && lvl <= 1000) {
		    if (lvl == 0) printLvl = 0;
		    else          printLvl = lvl;
		  }
		}
	      }
	    } else {
	      printUsage();
	      printedUsage = true;
	      exit(1);
	    }
	  }
	} else {
	  printUsage();
	  printedUsage = true;
	  exit(1);
	}
      }
    }
    
    double C[30];
    int c = 0;
    // HMWSoln *HMW = new HMWSoln("HMW_UO2_HNO3.xml", "HNO3_electrolyte");
    HMWSoln *HMW = new HMWSoln("HMW_TBP_UO2_HNO3.xml", "HNO3_electrolyte");
    IdealSolidSolnPhase *tbpLiquid = new IdealSolidSolnPhase("TBP_liquid.xml", "TBP_liquid");
//    IdealSolidSolnPhase *tbpLiquid = new IdealSolidSolnPhase("TBP_liquid_TBPorg_only.xml", "TBP_liquid");
 //   Zuzax::ThermoPhase *tbpLiquid = Zuzax::newPhase("TBP_liquid_StoichSubs.xml", "TBP_liquid");
//     Zuzax::ThermoPhase *tbpLiquid = Zuzax::newPhase("TBP_liquid_StoichSubs.xml", "TBP(org)");
//    Zuzax::ThermoPhase *tbpLiquid = Zuzax::newPhase("TBP_liquid_TBPorg_only.xml", "TBP_liquid");
    HMW->setState_TP(298.15, pres);
    // HMW->setState_TPX(298.15, pres, "H2O(L):0.9339, TBP(org):0.0661");
 //   tbpLiquid->setState_TPX(298.15, pres, "TBP(org):0.75");
 //   tbpLiquid->setState_TPX(298.15, pres, "dodecane:0.99, TBP(org):0.01");
 //   tbpLiquid->setState_TPX(298.15, pres, "AMSCO:0.3857, TBP(org):0.6143");   // for moles of org liquid equiv to [TBP] of 2.334 M (65 Vol.%) for TBP-AMSCO addPhase only 
//    tbpLiquid->setState_TPX(298.15, pres, "AMSCO:0.403, TBP(org):0.597");   // for moles of org liquid equiv to [TBP] of 2.334 M (65 Vol.%) for TBP-AMSCO addPhase only for UO2(NO3)2
 //   tbpLiquid->setState_TPX(298.15, pres, "AMSCO:0.7228, TBP(org):0.2772");   // for moles of org liquid equiv to [TBP] of 1.092 M (30 Vol.%) for TBP-AMSCO addPhase only 
    tbpLiquid->setState_TPX(298.15, pres, "AMSCO:0.725, TBP(org):0.275");   // for moles of org liquid equiv to [TBP] of 1.092 M (30 Vol.%) for TBP-AMSCO addPhase only for UO2(NO3)2
 //   tbpLiquid->setState_TPX(298.15, pres, "AMSCO:0.8694, TBP(org):0.1306");   // for moles of org liquid equiv to [TBP] of 0.522 M (15 Vol.%) for TBP-AMSCO addPhase only 
 //   tbpLiquid->setState_TPX(298.15, pres, "AMSCO:0.8697, TBP(org):0.1303");   // for moles of org liquid equiv to [TBP] of 0.515 M (15 Vol.%) for TBP-AMSCO addPhase only for UO2(NO3)2
//    tbpLiquid->setState_TPX(298.15, pres, "AMSCO:0.9124, TBP(org):0.0876");   // for moles of org liquid equiv to [TBP] of 0.3518 M (10 Vol.%) for TBP-AMSCO addPhase only 
//    tbpLiquid->setState_TPX(298.15, pres, "AMSCO:0.7559, TBP(org):0.2431");   // for 2:1 or 1:1 moles of addPhase for org liquid & electrolyte
//     tbpLiquid->setState_TPX(298.15, pres, "AMSCO:0.73, TBP(org):0.27");   // for 1:1 moles of addPhase for org liquid & electrolyte
 //   tbpLiquid->setState_TP(298.15, pres);

    Zuzax::MultiPhase mmm;
    
//    mmm.addPhase(HMW, 8.0);
//    mmm.addPhase(HMW, 5.0);
//    mmm.addPhase(HMW, 2.0);
    mmm.addPhase(tbpLiquid, 1.0);
    mmm.addPhase(HMW, 1.0);
//    mmm.addPhase(tbpLiquid, 0.5);

 // The following 4 lines below yield the same answer as setting mole fractions by setState_TPX(298.15, pres, "dodecane:0.7384, TBP(org):0.2616") 
 /*  mmm.addPhase(tbpLiquid, 0.2616);
    mmm.init();
    int iDodec = mmm.speciesIndex("dodecane", "TBP_liquid");
    mmm.addSpeciesMoles(iDodec, 0.7384);
*/

    vcs_equilibrate(mmm, "TP", estimateEquil, printLvl, solver);
    // equilibrate(mmm, "TP", solver);
    cout << mmm << endl;

//     C[c++] = 0.0;
    C[c++] = 0.01;
    C[c++] = 0.015;
    C[c++] = 0.15;
    C[c++] = 0.15;
    C[c++] = 0.15;
    C[c++] = 0.15;
    C[c++] = 0.10;
    C[c++] = 0.10;
    C[c++] = 0.10;
    C[c++] = 0.10;
    C[c++] = 0.10;
    C[c++] = 0.10; // CUT HERE FOR 15 VOL.% TBP
    C[c++] = 0.10;
    C[c++] = 0.10;
    C[c++] = 0.10;
    C[c++] = 0.10;
    C[c++] = 0.10;
    C[c++] = 0.10;
    C[c++] = 0.10;
/*    C[c++] = 0.25;  
    C[c++] = 0.25;
    C[c++] = 0.25;
    C[c++] = 0.15; // CUT HERE FOR 30 VOL.% TBP 
    C[c++] = 0.15; 
    C[c++] = 0.15; 
    C[c++] = 0.15; 
    C[c++] = 0.15; 
    C[c++] = 0.15;
    C[c++] = 0.148;
    C[c++] = 0.25; // CUT HERE FOR 15 and 30 VOL.% TBP
    C[c++] = 0.3; 
    C[c++] = 0.4;
    C[c++] = 0.6;
    C[c++] = 1.0;
    C[c++] = 1.2;
    C[c++] = 1.5;
    C[c++] = 1.8;
    C[c++] = 0.3;
    C[c++] = 0.5;
    C[c++] = 0.5;
    C[c++] = 0.5;
    C[c++] = 0.5;
    C[c++] = 1.0;
    C[c++] = 1.0;
    C[c++] = 1.0;
    C[c++] = 0.5;
    C[c++] = 1.2;
    C[c++] = 2.0;
    C[c++] = 2.5;
    C[c++] = 1.0;
*/

 //   int iTBPorg = mmm.speciesIndex("TBP(org)", "TBP_liquid"); 
    int iAMSCO = mmm.speciesIndex("AMSCO", "TBP_liquid"); 
 //   int idodec = mmm.speciesIndex("dodecane", "TBP_liquid"); 
//    int iTBPorg = mmm.speciesIndex("TBP(org)", "TBP(org)"); 
  //   int iH2O = mmm.speciesIndex("H2O(L)", "HNO3_electrolyte"); 
  //  int iNa = mmm.speciesIndex("Na+", "HNO3_electrolyte"); 
  //   int iCl = mmm.speciesIndex("Cl-", "HNO3_electrolyte"); 
  //  int iHaq = mmm.speciesIndex("H+", "HNO3_electrolyte"); 
    int iUO2aq = mmm.speciesIndex("UO2++", "HNO3_electrolyte"); 
    int iNO3aq = mmm.speciesIndex("NO3-", "HNO3_electrolyte"); 
    
    // for (int it = 0; it < 2; it++) {
    for (int i = 0; i < c; i++) {
     //   mmm.addSpeciesMoles(iNa, C[i] * 0.1);
   //     mmm.addSpeciesMoles(iCl, C[i] * 0.2);
    //  mmm.addSpeciesMoles(iTBPorg, C[i]);  // Increased the amount of titrated species to be used with "0" starting moles of HMW
   //   mmm.addSpeciesMoles(idodec, C[i]);  
   //    mmm.addSpeciesMoles(iH2O, C[i] * 70.);  // Increased the amount of titrated H2O(l) to maintain the overall H2O activity to unity
   //    mmm.addSpeciesMoles(iUO2aq, C[i] * 0.1);
   //    mmm.addSpeciesMoles(iHaq, C[i] * 0.1);
       mmm.addSpeciesMoles(iUO2aq, C[i] * 0.1);
       mmm.addSpeciesMoles(iNO3aq, C[i] * 2 * 0.1);
  //     mmm.addSpeciesMoles(iAMSCO, C[i] * 0.022); // For UO2(NO3)2 titration w/ prog. TBP dilution at 65 Vol. % TBP
   //    mmm.addSpeciesMoles(iAMSCO, C[i] * 0.013);  //For UO2(NO3)2 titration w/ prog. TBP dilution at 30 Vol. % TBP
       mmm.addSpeciesMoles(iAMSCO, C[i] * 0.024);  //For UO2(NO3)2 titration w/ prog. TBP dilution at 15 Vol. % TBP
      vcs_equilibrate(mmm, "TP", estimateEquil, printLvl, solver);
      cout << mmm << endl;
    }

    return 0;
  }
  
  catch (ZuzaxError) {
    showErrors(cerr);
    cerr << "program terminating." << endl;
    return -1;
  }
}
