/*
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "zuzax/equilibrium.h"
#include "zuzax/thermo/MolalityVPSSTP.h"

#include "zuzax/equil/vcs_prob.h"
#include "zuzax/numerics/NonlinearSolver_JAC.h"

//#include "Electrode_SimpleDiff.h"
//#include "Electrode_RadialDiffRegions.h"  
#include "zuzax/multiphase/PhaseList.h"
#include "ReactingSurDomain.h"

#include <sstream>
#include <iomanip>
#include <cstring>

using namespace Zuzax;
using namespace std;

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

  NonlinearSolver_JAC::s_TurnOffTiming = true;
  NonlinearSolver_JAC::s_print_NumJac = true;

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
     double Temp = 298.15;
     PhaseList*  pl = new PhaseList();
     pl->addVolPhase("metal_Li_LiIon_electrons.xml");
     pl->addVolPhase ("MCMB_RedlichKister.xml");
     pl->addVolPhase("ECsoln_ion.xml");
     pl->addVolPhase("Li_Metal.xml");
     pl->addSurPhase("MCMBAnode_electrode_extra.xml");

     pl->setState_TP(Temp, OneAtm);

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

     int iPhMetal =  pl->globalPhaseIndex("metal_Li_LiIon_electrons");
     ThermoPhase& metal = pl->thermo(iPhMetal);


     double nstoic = 1.0;
     double ocv = dg[1] / Faraday / nstoic;

     printf("ocv = %g\n", ocv);


     //
     //  We wil assume that the Li+ mole fraction in the electrolyte is this value
     //
     double xMoleFractionLip = 0.0780266;

     ThermoPhase& ecsoln = pl->thermo(2);
     xmol[1] = xMoleFractionLip;
     xmol[2] = xmol[1];
     xmol[0] = 1.0 - 2.0 * xmol[1]; 
     ecsoln.setState_TPX(Temp, OneAtm, xmol);
     double aa[10];
     ecsoln.getActivities(aa);
     printf("activity Li+ = %g\n", aa[1]);

     const string& rs =rsd->reactionString(1);

     int kLiC = mcmb.speciesIndex("Li_C6-bulk");
     int kC = mcmb.speciesIndex("V_C6-bulk");

     printf("                                 Thermo For Full-Cell Rxn, 1 : %s\n", rs.c_str());
     int numP = 51;
     printf ("Fig 3 Karthikeyan Sikha and White\n");
     printf("        xV              xLi OCV_HalfCellPlusLiP OCV_FullCell[1] DeltaH[1]     DeltaS[1]     DeltaG[1]   (DH - T DS)\n") ;
     for (int i =0; i < numP; i++) {
	 double xKC = 0.0 + (double) i / (numP - 1.0);
         xmol[kC] = xKC;
         xmol[kLiC] = 1.0 - xmol[kC];
         mcmb.setState_TPX(Temp, OneAtm, xmol);

	 rsd->getDeltaGibbs(dg);
	 double ocv = dg[1] / Faraday / nstoic;
         dg[0] -= GasConstant * Temp * std::log(xMoleFractionLip);
	 double ocvE = dg[0] / Faraday / nstoic;

         rsd->getDeltaEnthalpy(dh);
         rsd->getDeltaEntropy(ds);

	 printf(" %12.6f   %12.6f   %12.6f   %12.6f, %12.3E, %12.3E, %12.3E", 
		xKC, 1.0 - xKC, ocvE, ocv, dh[1], ds[1], dg[1]);
	 double dgCheck = dh[1] - Temp * ds[1];
	 printf("  %12.3E ", dgCheck);
	 printf("\n");
     }

     printf("  Analysis of the reaction as experienced at the start of the calculation\n");
     const string& rs0 =rsd->reactionString(0);

     xmol[1] = xMoleFractionLip;
     xmol[2] = xmol[1];
     xmol[0] = 1.0 - 2.0 * xmol[1]; 
     ecsoln.setState_TPX(Temp, OneAtm, xmol);
     printf("  Mole fraction of Li+ = %g\n", xmol[1]);
     ecsoln.getActivities(aa);
     printf("activity Li+ = %g\n", aa[1]);
     double liqS[10];
     ecsoln.getPartialMolarEntropies(liqS);
     double liqH[10];
     ecsoln.getPartialMolarEnthalpies(liqH);
     double liqG[10];
     ecsoln.getChemPotentials(liqG);
     printf("                                 Thermo For Half-Cell Rxn, Rxn 0: %s\n", rs0.c_str());
     
     printf("        xV              xLi      OCV_FullCell OCV_HalfCell[0]  DeltaH[0]    DeltaS[0]     DeltaG[0]     dgCheck\n") ;
     for (int i =0; i < numP; i++) {
	 double xKC = 0.0 + (double) i / (numP - 1.0);
         xmol[kC] = xKC;
         xmol[kLiC] = 1.0 - xmol[kC];
         mcmb.setState_TPX(Temp, OneAtm, xmol);

	 rsd->getDeltaGibbs(dg);
	 double ocv = dg[1] / Faraday / nstoic;
	 double ocvE = dg[0] / Faraday / nstoic;

         rsd->getDeltaEnthalpy(dh);
         rsd->getDeltaEntropy(ds);
	 double dgCheck = dh[0] - Temp * ds[0];

	 printf(" %12.6f   %12.6f   %12.6f   %12.6f, %12.3E, %12.3E, %12.3E %12.3E\n", 
		xKC, 1.0 - xKC,  ocv, ocvE, dh[0], ds[0], dg[0], dgCheck);
     }


     printf("\n\n SETTING xLi = 0.80 in the tables below\n\n");
     
     double RT = GasConstant * Temp; 
     xmol[kLiC] = 0.80;
     xmol[kC] = 1.0 -  xmol[kLiC];
     mcmb.setState_TPX(Temp, OneAtm, xmol);

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



     printf(" Reaction 1\n"); 
     printf(" Species    h(J/kmol)     h0(kJ/gmol)          s        s0        g          \n");
     printf("Li-C6    ");
     printf("% 12.3E % 12.3E % 12.3E % 12.3E % 12.3E % 12.3E  \n", hh[kLiC] , hh0[kLiC],  
            ss[kLiC] , ss0[kLiC] , gg[kLiC] , gg0[kLiC]  );
     printf(" V-C6    ");
     printf("% 12.3E % 12.3E % 12.3E % 12.3E % 12.3E % 12.3E \n", hh[kC] ,   hh0[kC] ,
            ss[kC] , ss0[kC] , gg[kC] , gg0[kC]  );
 
     double h_Li = Limetal.enthalpy_mole();
     double s_Li = Limetal.entropy_mole();
     double g_Li = Limetal.gibbs_mole();
     

     printf("LiMetal  % 12.3E % 12.3E % 12.3E % 12.3E % 12.3E % 12.3E \n", h_Li ,  h_Li , 
            s_Li ,  s_Li  ,  g_Li ,  g_Li );

     double deltaH = - hh[kLiC] +  hh[kC] + h_Li;
     double deltaS = - ss[kLiC] +  ss[kC] + s_Li;
     double deltaG = - gg[kLiC] +  gg[kC] + g_Li;
     printf("Total    %12.3E              % 12.3E              % 12.3E\n", deltaH, deltaS, deltaG);

     printf("\n\n");
     printf(" Reaction 0\n");
     printf(" Species    h(J/kmol)     h0(kJ/gmol)          s        s0        g       g0   \n");
     printf("Li-C6    ");
     printf("% 12.3E % 12.3E % 12.3E % 12.3E % 12.3E % 12.3E  \n", hh[kLiC] , hh0[kLiC] ,
            ss[kLiC] , ss0[kLiC] , gg[kLiC] , gg0[kLiC]);
     printf(" V-C6    ");
     printf("% 12.3E % 12.3E % 12.3E % 12.3E % 12.3E % 12.3E \n", hh[kC] ,   hh0[kC] ,
            ss[kC] , ss0[kC] , gg[kC] , gg0[kC]   );

     double h_e = metal.enthalpy_mole();
     double s_e = metal.entropy_mole();
     double g_e = metal.gibbs_mole();

     double h_Lip = liqH[1];
     double s_Lip = liqS[1];
     double g_Lip = liqG[1];


     printf("Li+(l)   % 12.3E % 12.3E % 12.3E % 12.3E % 12.3E % 12.3E \n", h_Lip, h_Lip ,
            s_Lip , s_Lip , g_Lip,  g_Lip );
     printf("Electron % 12.3E % 12.3E % 12.3E % 12.3E % 12.3E % 12.3E \n", h_e ,  h_e ,
            s_e ,  s_e  ,  g_e ,  g_e   );

     deltaH = -hh[kLiC] +  hh[kC] + h_Lip + h_e;
     deltaS = -ss[kLiC] +  ss[kC] + s_Lip + s_e;
     deltaG = -gg[kLiC] +  gg[kC] + g_Lip + g_e;
     printf("Total    % 12.3E              % 12.3E              % 12.3E \n", deltaH, deltaS ,   deltaG);


     double nTdeltaS =  - Temp * deltaS;
     
     double phiSoln = 1.6075E-2;
     double phiMetal = -5.4418E-04;
     ecsoln.setElectricPotential(phiSoln);
     metal.setElectricPotential(phiMetal);

     double icurD = rsd->getCurrentDensityRxn();

     double netROP[10];
     rsd->getNetRatesOfProgress(netROP);
     
     int irxn = 0;
     double nStoich;
     double OCV;
     double io;
     double nu;
     double beta, resist;
     //double icurD2 = rsd->getExchangeCurrentDensityFormulation(irxn, &nStoich, &OCV, &io, &nu, &beta, &resist);
     double icurD2;
#ifdef DONOTREMOVE
     icurD2 = rsd->getExchangeCurrentDensityFormulation(irxn, &nStoich, &OCV, &io, &nu, &beta, &resist);
#else
     bool okk = rsd->getExchangeCurrentDensityFormulation(irxn, nStoich, OCV, io, nu, beta, resist);
     if (okk) {
                icurD2 = rsd->calcCurrentDensity(nu, nStoich, io, beta, Temp, resist);
     } else {
	 icurD2 = 0.0;
     }
#endif
     printf("\nSample Rate Calculation Results:\n");
     printf("          T = %g \n", Temp);
     printf("          X_Li+ = %g\n", xMoleFractionLip);
     mcmb.getMoleFractions(xmol);
     printf("          X_LiC6 = %g\n", xmol[kLiC]);
     printf(" icurr  = %g (from getCurrentDensityRxn())\n", icurD);
     printf(" icurr2 = %g (from getExchangeCurrentDensityFormulation())\n", icurD2);
     printf("      Calculated exchange current density formulation:\n");
     printf("            nu = %g\n", nu);
     printf("            OCV = %g         -> volts = %g\n", OCV, nu + OCV);
     printf("            io =  %g\n", io);
     printf("            beta = %g\n", beta);
     printf("            nStoich = %g\n", nStoich);
     printf("            resist = %g\n", resist);
     double icalc =  rsd->calcCurrentDensity(nu, nStoich, io, beta, Temp);
     printf(" icurr3 = %g (from calcCurrentDensity()\n\n", icalc);


     double delx = (7.4667E-05 - 5.3333E-05) * 0.5;
     double sad = 2.25E5;  // surface area density
     double icurrRxnCell = sad * delx * 1.0E-4 * icurD  / 1.0E-4;

     printf("icurrRxnCell = %g amps / m2\n", icurrRxnCell);
  
     printf("-TdeltaS = %g Joule / kmol \n", nTdeltaS);

     double deltaT = 1.0E-8;

     double iTdelS = netROP[0] * nTdeltaS * deltaT * sad * delx;

     printf("Cell source term = %g J / m2\n", iTdelS);

     double volts = phiMetal - phiSoln;
     printf("volts = %g \n", volts);
     double iVTerm =  icurrRxnCell * volts * deltaT;
     double ndeltaHterm = - netROP[0] * deltaH * deltaT * sad * delx;
     double ndeltaHphiterm = iVTerm +  ndeltaHterm;
     printf("iVterm = %g\n", iVTerm);
     printf("ndeltaHterm = %g\n", ndeltaHterm);
     printf("ndeltaHphiterm = %g\n",  ndeltaHphiterm);

     double ndeltaGterm = - netROP[0] * deltaG * deltaT * sad * delx;
     printf("ndeltaGterm = %g\n", ndeltaGterm);
     double ndeltaGphiterm = iVTerm +  ndeltaGterm;
     printf("ndeltaGphiterm = %g\n",  ndeltaGphiterm);

     appdelete();

    return 0;

  } catch (ZuzaxError) {

    showErrors();
    return -1;
  }
} 
