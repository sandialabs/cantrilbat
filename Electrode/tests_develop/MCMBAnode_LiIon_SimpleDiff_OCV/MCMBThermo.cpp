/*
 * $Id: MCMBAnode_SimpleDiff.cpp 496 2013-01-07 21:15:37Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "cantera/equilibrium.h"
#include "cantera/thermo/MolalityVPSSTP.h"

#include "cantera/equil/vcs_prob.h"
#include "cantera/numerics/NonlinearSolver.h"

#include "Electrode_SimpleDiff.h"
#include "Electrode_RadialDiffRegions.h"  
#include "ExtraGlobalRxn.h"
#include "RxnMolChange.h"
#include "ReactingSurDomain.h"

#include <sstream>
#include <iomanip>

using namespace std;
using namespace Cantera;
using namespace VCSnonideal;

// a lvl of one prints out the .csv file
int mpequil_debug_print_lvl = 1;
int VCS_Debug_Print_Lvl = 3;

void printUsage() {
    cout << "usage: MCMBAnode_SimpleDiff [-h] [-help_cmdfile] [-d #] [anode.inp]"
         <<  endl;
    cout << "    -h               : Prints this help" << endl;
    cout << "    -help_cmdfile    : Prints a list of block commands understood by this parser - add anode.inp for more information" << endl;
    cout << "   -d #              : Level of debug printing" << endl;
    cout << "   anode.inp         : Command file (if missing, assume anode.inp)" << endl;
    cout << endl;
}

//=====================================================================================================


int main(int argc, char **argv)
{
 



  VCSnonideal::vcs_timing_print_lvl = 0;
  NonlinearSolver::s_TurnOffTiming = true;
  NonlinearSolver::s_print_NumJac = true;

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
		  mpequil_debug_print_lvl = lvl;
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

  try {

     PhaseList*  pl = new PhaseList();
     pl->addVolPhase("metal_Li_LiIon_electrons.xml");
     pl->addVolPhase ("MCMB_RedlichKister.xml");
     pl->addVolPhase("ECsoln_ion.xml");
     pl->addVolPhase("Li_Metal.xml");
     pl->addSurPhase("MCMBAnode_electrode.xml");

     pl->setState_TP(300.0, OneAtm);

     ReactingSurDomain* rsd = new ReactingSurDomain();
     rsd->importFromPL(pl, 0);

     double dg[5];
     rsd->getDeltaGibbs(dg);


     printf("dg[0] = %g\n", dg[0]);
     printf("dg[1] = %g\n", dg[1]);

     


     double nstoic = 1.0;
     double ocv = dg[1] / Faraday / nstoic;

     printf("ocv = %g\n", ocv);

     double xmol[10];

     ThermoPhase& mcmb = pl->thermo(1);

     ThermoPhase& ecsoln = pl->thermo(2);
     xmol[1] = 0.1;
     xmol[2] = xmol[1];
     xmol[0] = 1.0 - 2.0 * xmol[1]; 
     ecsoln.setState_TPX(300., OneAtm, xmol);
     double aa[10];
     ecsoln.getActivities(aa);
     printf("activity Li+ = %g\n", aa[1]);

     int kLiC = mcmb.speciesIndex("Li_C6-bulk");
     int KC = mcmb.speciesIndex("V_C6-bulk");
     int numP = 51;
     printf ("Fig 3 Karthikeyan Sikha and White\n");
     printf("        xLi              xV           OCV \n") ;
     for (int i =0; i < numP; i++) {
	 double xKC = 0.0 + (double) i / (numP - 1.0);
         xmol[KC] = xKC;
         xmol[kLiC] = 1.0 - xmol[KC];
         mcmb.setState_TPX(300., OneAtm, xmol);

	 rsd->getDeltaGibbs(dg);
	 double ocv = dg[1] / Faraday / nstoic;
         dg[0] -= GasConstant * 300. * std::log(0.1);
	 double ocvE = dg[0] / Faraday / nstoic;

	 printf(" %12.6f   %12.6f   %12.6f   %12.6f\n",  xKC, 1.0 - xKC,  ocv, ocvE);
     }
     


    Cantera::appdelete();

    return 0;

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
} 
