/*
 *  $Author: hkmoffa $
 *  $Date: 2013-01-07 15:07:12 -0700 (Mon, 07 Jan 2013) $
 *  $Revision: 500 $
 *
 *  Copyright 2002 California Institute of Technology
 *
 */

#include "cantera/IdealGasMix.h"
#include "cantera/equilibrium.h"

#include "cantera/thermo.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/base/logger.h"
#include "cantera/thermo/HMWSoln.h"

#include <cstring>

#ifdef useZuzaxNamespace
using namespace Zuzax;
#define ZZCantera Zuzax
#else
using namespace Cantera;
#define ZZCantera Cantera
#endif

using namespace std;

void printUsage() {
    cout << "usage: nacl_dessicate [-h] [-help_cmdfile] [-d #] "
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
    double pres = OneAtm;
    
#ifdef USE_VCSNONIDEAL
    solver = 2;
    //VCSnonideal::vcs_timing_print_lvl = 0;
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


    HMWSoln *HMW = new HMWSoln("HMW_BrineDatabase_CO2.xml", "carbonate_electrolyte");
    HMW->setState_TPX(300.0, pres, "H2O(L):0.8, Na+:0.1 Cl-:0.1 CO2(aq):1.0E-10");

    string nacl_s = "NaCl_Solid.xml";
    string id = "NaCl(S)";
    ZZCantera::ThermoPhase *solid = ZZCantera::newPhase(nacl_s, id);
    solid->setState_TP(300.0, pres);

    
    IdealGasMix gas("gas.xml", "air");
    gas.setState_TPX(300.0, pres, "H2O:0.12, CO2:0.88");



    ZZCantera::MultiPhase mmm;
    
    mmm.addPhase(HMW, 10.);
    mmm.addPhase(solid, 0.001);
    mmm.addPhase(&gas, 0.0000000000);

    //equilibrate(mmm, "TP", solver);

    vcs_equilibrate(mmm, "TP", estimateEquil, printLvl, solver);
    cout << mmm << endl;
    int iCO2g = mmm.speciesIndex("CO2", "air"); 
    int iH2OL = mmm.speciesIndex("H2O(L)", "carbonate_electrolyte"); 
    
    for (int it = 0; it < 10; it++) {
      mmm.addSpeciesMoles(iCO2g, 1.0);
      mmm.addSpeciesMoles(iH2OL, -1.0);
      vcs_equilibrate(mmm, "TP", estimateEquil, printLvl, solver);
      cout << mmm << endl;
    }
   
    return 0;
  }
  catch (CanteraError) {
    showErrors(cerr);
    cerr << "program terminating." << endl;
    return -1;
  }
}
