/*
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

//#include "Electrode_SimpleDiff.h"
//#include "Electrode_RadialDiffRegions.h"  
#include "PhaseList.h"
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
    cout << "    -help_cmdfile    : Prints a list of block commands understood by this parser "
      "- add anode.inp for more information" << endl;
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

     double xmol[10];
     PhaseList*  pl = new PhaseList();
     pl->addVolPhase("metal_Li_LiIon_electrons.xml");
     pl->addVolPhase ("MCMB_RedlichKister.xml");
     pl->addVolPhase("ECsoln_ion.xml");
     pl->addVolPhase("Li_Metal.xml");
     pl->addSurPhase("MCMBAnode_electrode_extra.xml");

     pl->setState_TP(300.0, OneAtm);

     ReactingSurDomain* rsd = new ReactingSurDomain();
     rsd->importFromPL(pl, 0);

     double dg[5], dh[5], ds[5], hh[5], ss[5], gg[5], hh0[5], ss0[5], gg0[5];
     rsd->getDeltaGibbs(dg);

     ThermoPhase& mcmb = pl->thermo(1);
     string nn = mcmb.id();
     printf("mcmb name = %s\n", nn.c_str());
     mcmb.getMoleFractions(xmol);


     printf("dg[0] = %g for xLi = %g\n", dg[0], xmol[0]);
     printf("dg[1] = %g for xLi = %g\n", dg[1], xmol[0]);


     int ii = pl->globalPhaseIndex("Li(S)");
     ThermoPhase& Limetal = pl->thermo(ii);
     


     double nstoic = 1.0;
     double ocv = dg[1] / Faraday / nstoic;

     printf("ocv = %g\n", ocv);



     ThermoPhase& ecsoln = pl->thermo(2);
     xmol[1] = 0.1;
     xmol[2] = xmol[1];
     xmol[0] = 1.0 - 2.0 * xmol[1]; 
     ecsoln.setState_TPX(300., OneAtm, xmol);
     double aa[10];
     ecsoln.getActivities(aa);
     printf("activity Li+ = %g\n", aa[1]);

     int kLiC = mcmb.speciesIndex("Li_C6-bulk");
     int kC = mcmb.speciesIndex("V_C6-bulk");
     int numP = 51;
     printf ("Fig 3 Karthikeyan Sikha and White\n");
     printf("        xV              xLi          OCV           OCV[1]           DeltaH         DeltaS       DeltaG\n") ;
     for (int i =0; i < numP; i++) {
	 double xKC = 0.0 + (double) i / (numP - 1.0);
         xmol[kC] = xKC;
         xmol[kLiC] = 1.0 - xmol[kC];
         mcmb.setState_TPX(300., OneAtm, xmol);

	 rsd->getDeltaGibbs(dg);
	 double ocv = dg[1] / Faraday / nstoic;
         dg[0] -= GasConstant * 300. * std::log(0.1);
	 double ocvE = dg[0] / Faraday / nstoic;

         rsd->getDeltaEnthalpy(dh);
         rsd->getDeltaEntropy(ds);

	 printf(" %12.6f   %12.6f   %12.6f   %12.6f, %12.3E, %12.3E, %12.3E\n", 
		xKC, 1.0 - xKC,  ocv, ocvE, dh[1], ds[1], dg[1]);
     }

     printf("  Analysis of the reaction as experienced at the start of the calculation\n");
     xmol[1] = 0.0780266;
     xmol[2] = xmol[1];
     xmol[0] = 1.0 - 2.0 * xmol[1]; 
     ecsoln.setState_TPX(300., OneAtm, xmol);
     printf("  Mole fraction of Li+ = %g\n", xmol[1]);
     ecsoln.getActivities(aa);
     printf("activity Li+ = %g\n", aa[1]);
     double liqS[10];
     ecsoln.getPartialMolarEntropies(liqS);
     
     printf("        xV              xLi          OCV           OCV[1]           DeltaH         DeltaS       DeltaG\n") ;
     for (int i =0; i < numP; i++) {
	 double xKC = 0.0 + (double) i / (numP - 1.0);
         xmol[kC] = xKC;
         xmol[kLiC] = 1.0 - xmol[kC];
         mcmb.setState_TPX(300., OneAtm, xmol);

	 rsd->getDeltaGibbs(dg);
	 double ocv = dg[1] / Faraday / nstoic;
         dg[0] -= GasConstant * 300. * std::log(0.1);
	 double ocvE = dg[0] / Faraday / nstoic;

         rsd->getDeltaEnthalpy(dh);
         rsd->getDeltaEntropy(ds);

	 printf(" %12.6f   %12.6f   %12.6f   %12.6f, %12.3E, %12.3E, %12.3E\n", 
		xKC, 1.0 - xKC,  ocv, ocvE, dh[0], ds[0], dg[0]);
     }


     printf("\n\n SETTING xLi = 0.68 in the tables below\n\n");
     
     double RT = GasConstant * 300.; 
     xmol[kLiC] = 0.68;
     xmol[kC] = 1.0 -  xmol[kLiC];
     mcmb.setState_TPX(300., OneAtm, xmol);

     mcmb.getPartialMolarEnthalpies(hh);
     mcmb.getEnthalpy_RT(hh0);
     hh0[0] *= RT;
     hh0[1] *= RT;
     mcmb.getPartialMolarEntropies(ss);
     mcmb.getEntropy_R(ss0);
     ss0[0] *= GasConstant; 
     ss0[1] *= GasConstant;
     mcmb.getChemPotentials(gg);
     mcmb.getGibbs_RT(gg0);
     gg0[0] *= RT;
     gg0[1] *= RT;


     
     printf(" Species    h(kJ/gmol)     h0(kJ/gmol)          s        s0        g          \n");
     printf("Li-C6    ");
     printf("%12.3E %12.3E %12.3E %12.3E %12.3E %12.3E  \n", hh[kLiC] * 1.0E-6, hh0[kLiC] * 1.0E-6, 
            ss[kLiC] * 1.0E-3, ss0[kLiC] * 1.0E-3, gg[kLiC] * 1.0E-6, gg0[kLiC] * 1.0E-6 );
     printf(" V-C6    ");
     printf("%12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n", hh[kC] * 1.0E-6,   hh0[kC] * 1.0E-6,
            ss[kC] * 1.0E-3, ss0[kC] * 1.0E-3, gg[kC] * 1.0E-6, gg0[kC] * 1.0E-6  );
 
     double h_Li = Limetal.enthalpy_mole();
     double s_Li = Limetal.entropy_mole();
     double g_Li = Limetal.gibbs_mole();
     

     printf("LiMetal    %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n", h_Li * 1.0E-6,  h_Li * 1.0E-6, 
            s_Li * 1.0E-3,  s_Li * 1.0E-3 ,  g_Li * 1.0E-6,  g_Li * 1.0E-6  );

     printf("Total   ");


    Cantera::appdelete();

    return 0;

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
} 
