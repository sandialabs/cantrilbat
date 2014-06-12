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

#include <fenv.h>

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

     fexcept_t ff;
     fegetexceptflag(&ff, FE_ALL_EXCEPT);
     printf("fexcept = %d\n", ff);
     PhaseList*  pl = new PhaseList();
     fegetexceptflag(&ff, FE_ALL_EXCEPT);
     printf("fexcept 2 = %d\n", ff);
     pl->addVolPhase("metal_Li_LiIon_electrons.xml");
     pl->addVolPhase("LiCoO2_RedlichKister.xml");
     pl->addVolPhase("ECsoln_ion.xml");
     pl->addVolPhase("Li_Metal.xml");
     pl->addSurPhase("LiCoO2Cathode_electrode_extra.xml");
     fegetexceptflag(&ff, FE_ALL_EXCEPT);
     printf("fexcept 2 = %d\n", ff);

     pl->setState_TP(300.0, OneAtm);

     ReactingSurDomain* rsd = new ReactingSurDomain();
     rsd->importFromPL(pl, 0);

     double dg[5];
     rsd->getDeltaGibbs(dg);
     fegetexceptflag(&ff, FE_ALL_EXCEPT);
     printf("fexcept = %d\n", ff);


     printf("dg[0] = %g\n", dg[0]);
     printf("dg[1] = %g\n", dg[1]);

     int res = feclearexcept(FE_ALL_EXCEPT);
     res = fegetexceptflag(&ff, FE_ALL_EXCEPT);

     
     double val = 2000.;
     double eval = std::exp(val);
     printf("eval = %g\n", eval);

     res = fegetexceptflag(&ff, FE_ALL_EXCEPT);

     double nstoic = 1.0;
     double ocv = dg[1] / Faraday / nstoic;

     printf("ocv = %g\n", ocv);

     double xmol[10];

     ThermoPhase& LiCoO2phase = pl->thermo(1);

     ThermoPhase& ecsoln = pl->thermo(2);
     xmol[1] = 0.1;
     xmol[2] = xmol[1];
     xmol[0] = 1.0 - 2.0 * xmol[1]; 
     ecsoln.setState_TPX(300., OneAtm, xmol);
     double aa[10];
     double ac[10];
     double lnac[10];
     ecsoln.getActivities(aa);
     printf("activity Li+ = %g\n", aa[1]);

     int kLiC = LiCoO2phase.speciesIndex("LiCoO2");
     int KC = LiCoO2phase.speciesIndex("CoO2");
     int numP = 51;
     printf ("Fig 3 Karthikeyan Sikha and White\n");
     printf("        xCoO2            x_LiCoO2           OCV_1     OCV_2  ac(CoO2)  ac(LiCoO2)  lnac(CoO2) lnac(LiCoO2) \n") ;
     double xStart = 1.0;
     for (int i =0; i < numP; i++) {
	 double xKC = xStart * (1.0 - (double) i / (numP ));
         xmol[KC] = xKC;
         xmol[kLiC] = 1.0 - xmol[KC];
         LiCoO2phase.setState_TPX(300., OneAtm, xmol);

	 rsd->getDeltaGibbs(dg);
	 double ocv = dg[1] / Faraday / nstoic;
         dg[0] -= GasConstant * 300. * std::log(0.1);
	 double ocvE = dg[0] / Faraday / nstoic;

         LiCoO2phase.getActivityCoefficients(ac);

       LiCoO2phase.getLnActivityCoefficients(lnac);



	 printf(" %12.6f   %12.6f   %12.6f   %12.6f %12.6E %12.6E %12.6E %12.6E \n",  xKC, 1.0 - xKC,  ocv, ocvE, ac[KC], ac[kLiC], lnac[KC], lnac[kLiC] );
     }
     

    Cantera::appdelete();

    return 0;

  } catch (CanteraError) {

    showErrors();
    return -1;
  } catch (std::exception &ee) {
    printf("standard exception caught\n");
    return -1;

  }
} 
