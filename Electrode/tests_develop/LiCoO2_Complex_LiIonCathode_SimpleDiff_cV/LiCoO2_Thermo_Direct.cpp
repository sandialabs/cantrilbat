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
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif
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

double voltsDirect(double x_li) 
{
    double x_alpha = 1.0 - x_li;
    if (x_alpha == 0.5) {
       x_alpha = 0.5000001;
    }
    if (x_alpha <= 0.0) {
      x_alpha = 1.0E-13;
    }
    if (x_alpha >= 1.0) {
      x_alpha = 0.9999999;
    }
    double A[10];
    A[0] =   6.4832E09;
    A[1] =        -6.5173E09;
    A[2] =             6.5664E09;
    A[3] =            -6.5787E09;
    A[4] =            6.3021E09;
    A[5] =           -5.0465E09;
    A[6] =            2.7113E09;
    A[7] =            -6.9045E08;
    double u0 = -29.614;
    
    double rt = GasConstant * 298.15;
    double rtf = rt / Faraday;
    int N = 8;
    double fac1 = 2.0 * x_alpha - 1.0;
    double fac2 = 1.0;
    double fac3 = fac1 * fac1;
    double lterm = log( (1.0 - x_alpha) / x_alpha);
    double ocv = u0 + lterm * rtf;
    for (int  k = 0; k < N; k++) {
	fac2 *= fac1;
	double firstT = fac2;
	fac3 = fac3 / fac1;
	double secT = 2.0 *  x_alpha * k * (1.0 -  x_alpha) / fac3;

	 ocv += A[k] * (firstT - secT) / Faraday;
    }
    return ocv;
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
     pl->addVolPhase("LiCoO2_RedlichKister.xml");
     pl->addVolPhase("ECsoln_ion.xml");
     pl->addVolPhase("Li_Metal.xml");
     pl->addSurPhase("LiCoO2Cathode_electrode_extra.xml");

     pl->setState_TP(300.0, OneAtm);

     ReactingSurDomain* rsd = new ReactingSurDomain();
     rsd->importFromPL(pl, 0);

     double dg[5];
     rsd->getDeltaGibbs(dg);


     //printf("dg[0] = %g\n", dg[0]);
     //printf("dg[1] = %g\n", dg[1]);

     


     double nstoic = 1.0;
     double ocv = dg[1] / Faraday / nstoic;

     //printf("ocv = %g\n", ocv);

     double xmol[10];

     ThermoPhase& LiCoO2phase = pl->thermo(1);

     ThermoPhase& ecsoln = pl->thermo(2);
     xmol[1] = 0.1;
     xmol[2] = xmol[1];
     xmol[0] = 1.0 - 2.0 * xmol[1]; 
     ecsoln.setState_TPX(300., OneAtm, xmol);
     double aa[10];
     ecsoln.getActivities(aa);
    // printf("activity Li+ = %g\n", aa[1]);

     int kLiC = LiCoO2phase.speciesIndex("LiCoO2");
     int KC = LiCoO2phase.speciesIndex("CoO2");
     int numP = 151;
     printf ("Fig 3 Karthikeyan Sikha and White\n");
     printf("        xLi              xV           OCV \n") ;
     
     for (int i =0; i < numP; i++) {
	 double xKC = 0.0 + (double) i / (numP - 1.0);
         xmol[KC] = xKC;
         xmol[kLiC] = 1.0 - xmol[KC];
	 ocv = voltsDirect( xmol[kLiC]);
	 printf(" %12.6f   %12.6f   %12.6f \n",  xKC, 1.0 - xKC,  ocv);
     }
 

    ZZCantera::appdelete();

    return 0;

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
} 
