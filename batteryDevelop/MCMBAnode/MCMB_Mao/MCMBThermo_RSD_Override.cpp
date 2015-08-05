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

#include "Electrode_input.h"
#include "Electrode_Exception.h"

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

     // FILL THIS IN
     OCV_Override_input ocvInput;

/*
 *  Set up the override by hand.
 *  the usual procedure is to do via an input deck like this
 *
start block Open Circuit Potential Override for interface  anode_surface
   Open Circuit Voltage Model = MCMB2528
   Replaced Species = Li_C6-bulk
   Identify Reaction for OCV Model = 0
   Temperature Derivative =  MODEL
   Temperature for OCV = 298.15
   Open Circuit Voltage Temperature Derivative Model = MCMB2528
end block Open Circuit Potential Override for interface anode_surface
*/
     ocvInput.replacedSpeciesName = "Li_C6-bulk";
     ocvInput.rxnID = 0;
     ocvInput.rxnID_deltaS = 1;
     ocvInput.OCVModel= "MCMB2528_Dualfoil";
     ocvInput.temperatureDerivType = 2; //"Model"
     ocvInput.temperatureBase = 298.15;
     ocvInput.OCVTempDerivModel = "MCMB2528_Dualfoil";
     ocvInput.DoDSurrogateSpeciesName = "V_C6-bulk";

     /*
      *  Post process the results
      */
     ocvInput.numTimes = 1;
     int kg = ocvInput.replacedGlobalSpeciesID = pl->globalSpeciesIndex(ocvInput.replacedSpeciesName);
     if (kg < 0) {
	 throw Electrode_Error("main", "Species not found in phaselist : "
			       + ocvInput.replacedSpeciesName);
     } 
     int phaseID;
     int localSpeciesIndex;
     pl->getLocalIndecisesFromGlobalSpeciesIndex(kg, phaseID, localSpeciesIndex);

     //
     // Store the phase index and local species index of the replaced species
     // 
     ocvInput.replacedSpeciesPhaseID = phaseID;
     ocvInput.replacedLocalSpeciesID = localSpeciesIndex;

     if (ocvInput.DoDSurrogateSpeciesName != "") {
	 size_t k =  pl->thermo(phaseID).speciesIndex(ocvInput.DoDSurrogateSpeciesName);
	 if (k != npos) {
	     ocvInput.MF_DoD_LocalSpeciesID = k;
	 } else {
	     exit(-1);
	 }
     }
     
     /*
      *  Post the override information - this changes the kinetics of the solid phase
      */
     rsd->addOCVoverride(&ocvInput);

     double dg[5], dh[5], ds[5], hh[5], ss[5], gg[5], hh0[5], ss0[5], gg0[5];

     ThermoPhase& mcmb = pl->thermo(1);
   
     int ii = pl->globalPhaseIndex("Li(S)");
     ThermoPhase& Limetal = pl->thermo(ii);

     int iPhMetal =  pl->globalPhaseIndex("metal_Li_LiIon_electrons");
     ThermoPhase& metal = pl->thermo(iPhMetal);


     double nstoic = 1.0;



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
     
     /*
      *  There are two reactions in the mechanism 0 - half cell reaction
      *                                           1 - full cell reaction with zero rate constant.
      */

     printf("                                 Thermo For Full-Cell Rxn, 1 : %s\n", rs.c_str());
     int numP = 51;
     printf("                      - notes = half cell reaction with LiP contribution is equal to the full cell reaction\n");
     printf("                             \n");
     printf("        xV              xLi OCV_HalfCellPlusLiP OCV_FullCell[1] DeltaH[1]     DeltaS[1]     DeltaG[1]   (DH - T DS)\n") ;
     for (int i = -1; i < numP; i++) {
         double xKC;
         if (i == -1) {
             xKC = 0.3; 
         } else {
	     xKC = 0.0 + (double) i / (numP - 1.0);
         }
         xmol[kC] = xKC;
         xmol[kLiC] = 1.0 - xmol[kC];
         mcmb.setState_TPX(Temp, OneAtm, xmol);

	 rsd->getDeltaGibbs(dg);
	 double ocv = dg[1] / Faraday / nstoic;
         dg[0] -= GasConstant * Temp * std::log(xMoleFractionLip);
	 double ocvE = dg[0] / Faraday / nstoic;

         rsd->getDeltaEnthalpy(dh);
         rsd->getDeltaEntropy(ds);

	 printf(" %12.6f   %12.6f   %12.7f   %12.7f, %12.5E, %12.5E, %12.5E", 
		xKC, 1.0 - xKC, ocvE, ocv, dh[1], ds[1], dg[1]);
	 double dgCheck = dh[1] - Temp * ds[1];
	 printf("  %12.5E ", dgCheck);
	 printf("\n");
     }

     printf("    What's suppose to occur:   The last two columns of the table are suppuse to be equal and \n");
     printf("    The OCV value is suppose to agree with Fig. 8 Mao et al.\n");

     printf("\n\n");
     printf("         xV              xLi OCV[1] dOCVdT[1] \n");
       
      for (int i =0; i < numP; i++) {
         
         double xKC = 0.0 + (double) i / (numP - 1.0);
         xmol[kC] = xKC;
         xmol[kLiC] = 1.0 - xmol[kC];
         mcmb.setState_TPX(Temp, OneAtm, xmol);

         rsd->getDeltaGibbs(dg);
         double ocv300 = dg[1] / Faraday / nstoic;
         dg[0] -= GasConstant * Temp * std::log(xMoleFractionLip);
         //double ocvE = dg[0] / Faraday / nstoic;

         rsd->getDeltaEnthalpy(dh);
         rsd->getDeltaEntropy(ds);

         double dOCVdt = - ds[1] / Faraday;  

         printf(" %12.6f   %12.6f   %12.7f   %12.7f",
                xKC, 1.0 - xKC, ocv300, dOCVdt*1.0E3);
         printf("\n");
     }



     const string& rs0 =rsd->reactionString(0);

     xmol[1] = xMoleFractionLip;
     xmol[2] = xmol[1];
     xmol[0] = 1.0 - 2.0 * xmol[1]; 
     ecsoln.setState_TPX(Temp, OneAtm, xmol);
     printf("  Mole fraction of Li+ = %g\n", xmol[1]);
     ecsoln.getActivities(aa);
     //printf("activity Li+ = %g\n", aa[1]);
     double liqS[10];
     ecsoln.getPartialMolarEntropies(liqS);
     double liqH[10];
     ecsoln.getPartialMolarEnthalpies(liqH);
     double liqG[10];
     ecsoln.getChemPotentials(liqG);
     printf("                                 Thermo For Half-Cell Rxn, Rxn 0: %s\n", rs0.c_str());
     
     printf("        xV              xLi      OCV_FullCell OCV_HalfCell[0]  DeltaH[0]    DeltaS[0]     DeltaG[0]     dgCheck\n") ;
     printf("                                            (with Lip additions)\n");
     for (int i =0; i < numP; i++) {
	 double xKC = 0.0 + (double) i / (numP - 1.0);
         xmol[kC] = xKC;
         xmol[kLiC] = 1.0 - xmol[kC];
         mcmb.setState_TPX(Temp, OneAtm, xmol);

	 rsd->getDeltaGibbs(dg);
	 double ocv = dg[1] / Faraday / nstoic;
         //double dg0 = dg[0] - GasConstant * Temp * std::log(xMoleFractionLip);
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

     double deltaG_species,deltaH_species, deltaS_species;
     rsd->getOCVThermoOffsets_ReplacedSpecies(deltaG_species, deltaH_species, deltaS_species);

     printf("\n");
     printf("           DETAILED LOOK AT EACH REACTION \n");
     printf("\n");

     printf(" Reaction 1\n"); 
     printf(" Species    h(J/kmol)     h0(kJ/gmol)          s        s0        g         g0 \n");
     printf("Li-C6    ");
     printf("% 12.3E % 12.3E % 12.3E % 12.3E % 12.3E % 12.3E  \n", hh[kLiC] , hh0[kLiC], ss[kLiC] , ss0[kLiC] , 
            gg[kLiC] , gg0[kLiC]  );
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
     printf("Dsp      %12.3E              % 12.3E              % 12.3E\n", deltaH_species,  deltaS_species,  deltaG_species);
     deltaH -= deltaH_species;
     deltaS -= deltaS_species;
     deltaG -= deltaG_species;
     printf("revTotal %12.3E              % 12.3E              % 12.3E\n", deltaH, deltaS, deltaG);

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
     printf("Dsp      %12.3E              % 12.3E              % 12.3E\n", deltaH_species,  deltaS_species,  deltaG_species);
     deltaH -= deltaH_species;
     deltaS -= deltaS_species;
     deltaG -= deltaG_species;
     printf("revTotal %12.3E              % 12.3E              % 12.3E\n", deltaH, deltaS, deltaG);

     printf("\n");
     printf("Now let's look at a single reaction point:\n");
     double nTdeltaS =  - Temp * deltaS;
     
     double phiSoln = -0.1;
     double phiMetal = 0.109659;
     double volts = phiMetal - phiSoln;
     printf(" phiMetal = %g, phiSoln = %g , Voltage = %g\n", phiMetal, phiSoln, volts);
     ecsoln.setElectricPotential(phiSoln);
     metal.setElectricPotential(phiMetal);

     double icurD = rsd->getCurrentDensityRxn();

     const std::vector<double>& netRate = rsd->calcNetSurfaceProductionRateDensities();
     double netROP[10];
     rsd->getNetRatesOfProgress(netROP);
     
     int irxn = 0;
     double nStoich;
     double OCV;
     double io;
     double nu;
     double beta, resist;
     double icurD2 = rsd->getExchangeCurrentDensityFormulation(irxn, &nStoich, &OCV, &io, &nu, &beta, &resist);

     printf(" icurr = %g   icurr2 = %g \n", icurD, icurD2);
     printf("      Calculated exchange current density formulation:\n");
     printf("            nu = %g\n", nu);
     printf("            OCV = %g         -> volts = %g\n", OCV, nu + OCV);
     printf("            io =  %g\n", io);
     printf("            beta = %g\n", beta);
     printf("            nStoich = %g\n", nStoich);
     printf("            resist = %g\n", resist);
     double icalc =  rsd->calcCurrentDensity(nu, nStoich, io, beta, Temp);
     printf("            i = %g\n", icalc);

     Temp = 298.15 + 100;
     rsd->setState_TP(Temp, OneAtm);


     printf("\n\n");
     printf("         xV              xLi OCV[1] dOCVdT[1] \n");
       
     for (int i =0; i < numP; i++) {
         
         double xKC = 0.0 + (double) i / (numP - 1.0);
         xmol[kC] = xKC;
         xmol[kLiC] = 1.0 - xmol[kC];
         mcmb.setState_TPX(Temp, OneAtm, xmol);

         rsd->getDeltaGibbs(dg);
         double ocv300 = dg[1] / Faraday / nstoic;
         dg[0] -= GasConstant * Temp * std::log(xMoleFractionLip);
         //double ocvE = dg[0] / Faraday / nstoic;

         rsd->getDeltaEnthalpy(dh);
         rsd->getDeltaEntropy(ds);

         double dOCVdt = - ds[1] / Faraday;  

         printf(" %12.6f   %12.6f   %12.7f   %12.7f",
                xKC, 1.0 - xKC, ocv300, dOCVdt*1.0E3);
         printf("\n");
     }

     /*
      *  There are two reactions in the mechanism 0 - half cell reaction
      *                                           1 - full cell reaction with zero rate constant.
      */

     printf("                                 Thermo For Full-Cell Rxn, 1 (Redo at 398.15 K: %s\n", rs.c_str());
     numP = 51;
     printf("                      - notes = half cell reaction with LiP contribution is equal to the full cell reaction\n");
     printf("                             \n");
     printf("        xV              xLi OCV_HalfCellPlusLiP OCV_FullCell[1] DeltaH[1]     DeltaS[1]     DeltaG[1]   (DH - T DS)\n") ;
     for (int i = -1; i < numP; i++) {
         double xKC;
         if (i == -1) {
             xKC = 0.3; 
         } else {
	     xKC = 0.0 + (double) i / (numP - 1.0);
         }
         xmol[kC] = xKC;
         xmol[kLiC] = 1.0 - xmol[kC];
         mcmb.setState_TPX(Temp, OneAtm, xmol);

	 rsd->getDeltaGibbs(dg);
	 double ocv = dg[1] / Faraday / nstoic;
         dg[0] -= GasConstant * Temp * std::log(xMoleFractionLip);
	 double ocvE = dg[0] / Faraday / nstoic;

         rsd->getDeltaEnthalpy(dh);
         rsd->getDeltaEntropy(ds);

	 printf(" %12.6f   %12.6f   %12.7f   %12.7f, %12.5E, %12.5E, %12.5E", 
		xKC, 1.0 - xKC, ocvE, ocv, dh[1], ds[1], dg[1]);
	 double dgCheck = dh[1] - Temp * ds[1];
	 printf("  %12.5E ", dgCheck);
	 printf("\n");
     }




     Cantera::appdelete();

    return 0;

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
} 
