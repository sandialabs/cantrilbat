/*
 *  $Author: hkmoffa $
 *  $Date: 2013-01-07 15:24:38 -0700 (Mon, 07 Jan 2013) $
 *  $Revision: 501 $
 *
 *
 */

#include "cantera/IdealGasMix.h"
#include "cantera/equilibrium.h"
#include "cantera/equil/vcs_MultiPhaseEquil.h"
#include "cantera/equil/vcs_internal.h"

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
    cout << "usage: silane_equil [-h] [-help_cmdfile] [-d #] "
         <<  endl;
    cout << "    -h           help" << endl;
    cout << "    -d           #   : level of debug printing" << endl;
    cout << "  vcs_Cantera.inp    : command file" << endl;
    cout << "                     : (if missing, assume vcs_Cantera.inp)"
         << endl;
    cout << endl;
}

int main(int argc, char **argv) {
  int numSucc = 0;
  int numFail = 0;
  int printLvl = 1;
  int estimateEquil = -1;
  string inputFile = "gri30.xml";
  bool printInputFormat = false; // print cmdfile.txt format 
  bool printedUsage = false;

  //VCSnonideal::vcs_timing_print_lvl = 0;

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
      } else if (inputFile == "gri30.xml") {
	inputFile = tok;
      } else {
	printUsage();
	printedUsage = true;
	exit(1);
      }
    }
  }



  try {
    int retnSub;
    double T = 2000;
    int solver = 2;


    IdealGasMix g(inputFile.c_str(), "");
    double pres = OneAtm;

    int kk = g.nSpecies();
    vector_fp Xmol(kk, 0.0);
    int iC2H4 = g.speciesIndex("C2H4");
    int iO2 = g.speciesIndex("O2");
    int iN2 = g.speciesIndex("N2");

    /*
     * Do an initial calculation that can be debugged
     * easily
     */
    Xmol[iC2H4] = 0.4;
    Xmol[iO2] = 0.4;
    Xmol[iN2] = 0.2;

    g.setState_TPX(T, pres, DATA_PTR(Xmol));
    try {
      retnSub = ZZCantera::vcs_equilibrate(g, "TP", estimateEquil, printLvl, solver, 1.0e-9, 1000, 100, -99);

      cout << g;
      if (retnSub != 1) {
	cerr << "ERROR: MultiEquil equilibration step failed at " 
	     << " T    = " << T 
	     << " Pres = " << pres 
	     << endl;
	cout << "ERROR: MultiEquil equilibration step failed at " 
	     << " T    = " << T 
	     << " Pres = " << pres
	     << endl;
	exit(-1);
      }
      numSucc++;
    } catch (CanteraError) {
      cout << g;
      showErrors(cerr);
      cerr << "ERROR: MultiEquil equilibration step failed at " 
	   << " T    = " << T 
	   << " Pres = " << pres 
	   << endl;
      cout << "ERROR: MultiEqiul equilibration step failed at " 
	   << " T    = " << T 
	   << " Pres = " << pres
	   << endl;
      exit(-1);
    }

    cout << "NUMBER OF SUCCESSES =  " << numSucc << endl;
    cout << "NUMBER OF FAILURES  =  " << numFail << endl;

    return numFail;
  }
  catch (CanteraError) {
    showErrors(cerr);
    cerr << "ERROR: program terminating due to unforeseen circumstances." << endl;
    return -1;
  }
}
