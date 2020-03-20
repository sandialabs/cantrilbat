/*
 * $Id: mpequil_main.cpp 502 2013-01-07 22:25:47Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */


#include "zuzax/equilibrium.h"
#include "zuzax/base/stringUtils.h"
#include "mpequil_input.h"
#include "mpequil_prep.h"

#ifndef USE_VCSNONIDEAL
#define USE_VCSNONIDEAL
#endif

#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y)    (( (x) < (y) ) ? (x) : (y))
#endif

#include "zuzax/equil/vcs_MultiPhaseEquil.h"
#include "zuzax/equil/vcs_prob.h"
#include "zuzax/equil/vcs_solve.h"
#include "zuzax/equil/vcs_VolPhase.h"
#include "zuzax/equil/vcs_internal.h"

#include <cstring>

using namespace std;
using namespace Zuzax;

using namespace vcs_nonideal;

// a lvl of one prints out the .csv file
int mpequil_debug_print_lvl = 1;
int VCS_Debug_Print_Lvl = 3;

void printUsage() {
    cout << "usage: mpequil [-h] [-help_cmdfile] [-d #] [mpequil.inp]"
         <<  endl;
    cout << "    -h           help" << endl;
    cout << "    -d           #   : level of debug printing" << endl;
    cout << "  mpequil.inp    : command file" << endl;
    cout << "                     : (if missing, assume mpequil.inp)" 
	 << endl;
    cout << endl;
}


//! This routine converts the MPEQUIL_INPUT problem into a VCS_PROB
//! structure
#ifdef USE_VCSNONIDEAL
int mpequil_convert(MPEQUIL_INPUT *prob_input, vcs_nonideal::VCS_PROB *vprob) {

  Zuzax::MP_EquilStatic *mix = prob_input->m_mp;
  /*
   *     Extract the current state information
   *     from the MultiPhase object and
   *     Transfer it to VCS_PROB object.
   */
  int res = vcs_Cantera_to_vprob(mix, vprob);
  if (res != 0) {
    printf("problems\n");
  }

  // Set the estimation technique
  if (prob_input->iest) {
    vprob->iest = 0;
  } else {
    vprob->iest = -1;
  }

  // Figure out where the element abundances are
  // coming from. Is it from the multiphase object are are
  // they specified?
  // mix->getElementAbundances()
  // mix->nElements()
  // number of element constraints in vprob may be different than the
  // number of elements in the mixture. How do we do the mapping?
  // In general the number of element constraints in VCS_PROB is larger
  // than the number of elements, because there may be multiple oxidation
  // 
  if (prob_input->specifiedElementAbundances) {
    for (size_t e = 0; e < vprob->ne; e++) {
      // gather some useful info about the element constraint
      //int elType = vprob->m_elType[e];
      //int elActive = vprob->ElActive[e];
      std::string elName = vprob->ElName[e];
      // Zero out the element abundance for this element.
      vprob->m_elemAbund[e] = 0.0;
      // Now find the equivalent element in the multiphase object by comparing names
      int mixNE = mix->nElements();
      for (int ie = 0; ie < mixNE; ie++) {
	std::string mname = mix->elementName(ie);
	if (elName == mname) {
	  //! override the element abundances
	  vprob->m_elemAbund[e] = prob_input->elementMoles[ie];
	}
      }
    }
  }

  return 0;
}
#endif

//! This routine takes the problem as specified in the MPEQUIL_INPUT structure
//! and solves it to equilibrium
/*!
 *
 */
#ifdef USE_VCSNONIDEAL
int  mpequil_equilibrate(MPEQUIL_INPUT *prob_input, int estimateInit, int printFlag) {

  MP_EquilStatic *mix = prob_input->m_mp;
  /*
   * Malloc a new vcs problem strcutre
   */
  vcs_nonideal::VCS_PROB *vprob =
    new vcs_nonideal::VCS_PROB(prob_input->nspecies, prob_input->ne, prob_input->nphase);
#ifdef DEBUG
  vprob->setDebugPrintLvl(VCS_Debug_Print_Lvl);
#endif
  /*
   * vcs problem input to Zuzax conversion
   */
  int res = mpequil_convert(prob_input, vprob);
  if (res != 0) {
    printf("problems\n");
  }

  /*
   * Print out the problem specification from the point of
   * view of the vprob object.
   */
  vprob->prob_report(2);

  /*
   * Call the thermo Program
   */
  int ip1 = 0;
   if (printFlag >= 3) {
     ip1 = printFlag - 2;
   } 
  int ipr = MAX(0, printFlag - 1);
  int maxit = 1000;
  vcs_nonideal::VCS_SOLVE *vsolvePtr = new vcs_nonideal::VCS_SOLVE();
  int iconv = vsolvePtr->vcs(vprob, 0, ipr, ip1, maxit);


  /*
   * Transfer the information back to the MultiPhase object.
   * Note we don't just call setMoles, because some multispecies
   * solution phases may be zeroed out, and that would cause a problem
   * for that routine. Also, the mole fractions of such zereod out
   * phases actually contain information about likely reemergent
   * states.
   */
  mix->uploadMoleFractionsFromPhases();
  int kGlob = 0;
  for (size_t ip = 0; ip < vprob->NPhase; ip++) {
    double phaseMole = 0.0;
    ThermoPhase &tref = mix->phase(ip);
    int nspPhase = tref.nSpecies();
    for (int k = 0; k < nspPhase; k++, kGlob++) {
      phaseMole += vprob->m_spMoles[kGlob];
    }
    //phaseMole *= 1.0E-3;
    mix->setPhaseMoles(ip, phaseMole);
  }
  
  //   double te = vcs_second();
  if (printFlag > 0) {
    printf("\n Results from vcs:\n");
    if (iconv != 0) {
      printf("\nVCS FAILED TO CONVERGE!\n");
    }
    printf("\n");
    printf("Temperature = %g Kelvin\n",  vprob->temperature());
    printf("Pressure    = %g Pa\n", vprob->PresPA);
    printf("\n");
    printf("----------------------------------------"
	   "---------------------\n");
    printf(" Name            KMole_Number");
    printf("(kmol)");
    printf("  Mole_Fraction     Chem_Potential");
    if (vprob->m_VCS_UnitsFormat == VCS_UNITS_KCALMOL) 
      printf(" (kcal/mol)\n");
    else if (vprob->m_VCS_UnitsFormat == VCS_UNITS_UNITLESS) 
      printf(" (Dimensionless)\n");
    else if (vprob->m_VCS_UnitsFormat == VCS_UNITS_KJMOL) 
      printf(" (kJ/mol)\n");
    else if (vprob->m_VCS_UnitsFormat == VCS_UNITS_KELVIN) 
      printf(" (Kelvin)\n");
    else if (vprob->m_VCS_UnitsFormat == VCS_UNITS_MKS) 
      printf(" (J/kmol)\n");
    printf("--------------------------------------------------"
	   "-----------\n");
    for (size_t i = 0; i < vprob->nSpecies(); i++) {
      printf("%-12s", vprob->SpName[i].c_str());
      if (vprob->SpeciesUnknownType[i] == 
	  VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	printf("  %15.3e %15.3e  ", 0.0, vprob->mf[i]);
	printf("%15.3e\n", vprob->m_gibbsSpecies[i]);
      } else {
	printf("  %15.3e   %15.3e  ", vprob->m_spMoles[i], vprob->mf[i]);
	if (vprob->m_spMoles[i] <= 0.0) {
	  int iph = vprob->PhaseID[i];
	  vcs_VolPhase *VPhase = vprob->VPhaseList[iph];
	  if (VPhase->nSpecies() > 1) {
	    printf("     -1.000e+300\n");
	  } else {
	    printf("%15.3e\n", vprob->m_gibbsSpecies[i]);
	  }
	} else {
	  printf("%15.3e\n", vprob->m_gibbsSpecies[i]);
	}
      }
    }
    printf("------------------------------------------"
	   "-------------------\n"); 

    //     printf("Total time = %12.6e seconds\n", te - ts);
  }
  // hard code a csv output file.
  if (printFlag > 0) {
    static int counter = 0;
    string reportFile = "vcs_equilibrate_res.csv";
    if (counter > 0) {
      reportFile = "vcs_equilibrate_res_" + int2str(counter) + ".csv";
    }
    vprob->reportCSV(reportFile);
    counter++;
  }
  delete vsolvePtr;
  delete vprob;
  return iconv;
}
#endif

int main(int argc, char **argv)
{
  int i;
  FILE *inputFP = stdin;
  MPEQUIL_INPUT *prob_input = new MPEQUIL_INPUT();
  int retn;
  string commandFile = "mpequil.inp";
  // printed usage

  vcs_nonideal::vcs_timing_print_lvl = 0;

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
      } else if (commandFile == "" || commandFile == "mpequil.inp") {
	commandFile = tok;
      } else {
	printUsage();
	exit(1);
      }
    }
  }

  if (commandFile != "") {
    inputFP = fopen(commandFile.c_str(), "r");
    if (!inputFP) {
      printf("Can't open file, %s, bailing\n", commandFile.c_str());
      return -1;
    }
  }


  try {
    /*
     * Go get the problem description from the input file
     */
    retn = mpequil_input(prob_input, commandFile);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }
 
    /*
     * prep the input
     */
    int res = mpequil_prep(prob_input);
    if (res != 0) {
      printf("problems\n");
    }

    /*
     * Print out the problem specification
     */
    //vprob->prob_report(2);

    /*
     * Call the thermo Program
     */
    MP_EquilStatic *mp =  prob_input->m_mp;

#ifdef USE_VCSNONIDEAL
    int printFlag = mpequil_debug_print_lvl;
    int estimateInit = 0;
    //vcs_MultiPhaseEquil mpe(*mp);
    //(void) mpe.equilibrate(prob_input->prob_type, estimateInit, printFlag);
    // Call this one because we output a csv file. This may be used 
    // in a test suite
    if (prob_input->specifiedElementAbundances) {

      mpequil_equilibrate(prob_input, estimateInit, printFlag);
   
    } else {
      vcs_equilibrate_1(*mp, prob_input->prob_type, estimateInit, printFlag);
    }
#else
    MultiPhaseEquil mpe(mp);
	 
    (void) mpe.equilibrate(prob_input->prob_type);
#endif
     
    mpequil_query(prob_input);

    printf("\n Results from mpequil:\n");
    printf("\n");
    printf("Temperature = %g Kelvin\n",  prob_input->T);
    printf("Pressure    = %g ", prob_input->Pres);
    printf("Pa\n");
 
    printf("\n");
    printf("----------------------------------------"
	   "---------------------\n");
    printf(" Name             Mole_Number     Mole_Fraction     Chem_Potential");
    printf(" (J/kmol)\n");
    printf("--------------------------------------------------"
	   "-----------\n");
    string spname;
    for (i = 0; i < prob_input->nspecies; i++) {
      spname = mp->speciesName(i);
      printf("%-12s", spname.c_str());
      printf("  %15.6e %15.6e  %15.6e\n", prob_input->spMoles[i],
	     prob_input->spMf[i], prob_input->spChemPot[i]);
    }
    printf("------------------------------------------"
	   "-------------------\n");   

    int nph = mp->nPhases();
    for (int iph = 0; iph < nph; iph++) {
      ThermoPhase & tp = mp->phase(iph);
      cout << tp.report() << endl;
    }
   
    delete prob_input;


    Zuzax::appdelete();

    return retn;

  } catch (ZuzaxError) {

    showErrors();
    return -1;
  }
} 
