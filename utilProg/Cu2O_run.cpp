/*
 *
 */

#include "cantera/base/ct_defs.h"
#include "cantera/equilibrium.h"
#include "cantera/thermo.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/thermo/HMWSoln.h"
#include "cantera/IdealGasMix.h"

#include <cstdio>
#include <cstring>

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
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



    double pres = OneAtm;
    
#ifdef USE_VCSNONIDEAL
    ZZVCSnonideal::vcs_timing_print_lvl = 0;
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
	  } else if (tok[n] == 'h') {
	    printUsage();
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


//    HMWSoln *HMW = new HMWSoln("HMW_NaCl.xml", "NaCl_electrolyte");
    HMWSoln *HMW = new HMWSoln("HMW_CuNaCl.xml", "brine_electrolyte");
    HMW->setState_TPX(298.15, pres, "H2O(L):0.85 Cu++:0.05 Cl-:0.1");

    double totalElementEntrop = 0.0;

    size_t m = HMW->elementIndex("Cu");
    double se_Cu = HMW->entropyElement298(m);
    totalElementEntrop += se_Cu * 1;

    m = HMW->elementIndex("O");
    double se_O = HMW->entropyElement298(m);
    totalElementEntrop += se_O * 2;

    m = HMW->elementIndex("H");
    double se_H = HMW->entropyElement298(m);
    totalElementEntrop += se_H * 2;


    m = HMW->elementIndex("Cl");
    double se_Cl = HMW->entropyElement298(m);
 

    printf("totalElementEntrop = %g\n", totalElementEntrop);
    double s0_sp = 87.00E3;
    double deltaGf = -359.92E6;

    double Hf = deltaGf  - 298.15 * totalElementEntrop + 298.15 * s0_sp;

    printf ("Hf = %g kJ / gmol \n", Hf / 1.0E6);

    double a_shomate = 48.597;
    double b_shomate = 7.427;
    double c_shomate = 0.0;
    double d_shomate = 0.0;
    double e_shomate = -0.761;

    double t = 298.15 / 1000.0;
     //  {h}^0(T) = A t + \frac{B t^2}{2} + \frac{C t^3}{3} + \frac{D t^4}{4}  - \frac{E}{t}  + F. 

    double f_shomate  = Hf / 1.0E6  - a_shomate * t - b_shomate * t * t / 2.0  + e_shomate / t;
    
    printf("f_shomate = %g\n", f_shomate);

    //  {s}^0(T) = A\ln t + B t + \frac{C t^2}{2} + \frac{D t^3}{3} - \frac{E}{2t^2}  + G.

    double g_shomate = s0_sp / 1.0E3 - a_shomate * log(t) - b_shomate * t + e_shomate / ( 2.0 * t * t);

    printf("g_shomate = %g\n", g_shomate);

    printf("Shomate Polynomials:\n");
    printf("    %14.6E , %14.6E , %14.6E , %14.6E , %14.6E , %14.6E, %14.6E \n", a_shomate,  b_shomate ,  c_shomate , d_shomate , e_shomate ,  f_shomate ,  g_shomate );


    
    
    return 0;
  }
  catch (CanteraError) {
    showErrors(cerr);
    cerr << "program terminating." << endl;
    return -1;
  }
}
