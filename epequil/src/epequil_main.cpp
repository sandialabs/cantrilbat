/*
 * $Id: epequil_main.cpp 508 2013-01-07 22:54:04Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */


#include "zuzax/equilibrium.h"
#include "epequil_input.h"
#include "epequil_prep.h"

#include <cstdio>
#include <cstring>

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif
using namespace std;

#ifdef DEBUG
int epequil_debug_print_lvl = 2;
#endif

#ifdef DEBUG_HKM_EPEQUIL
extern int debug_prnt_lvl;
#endif

void printUsage() {
    cout << "usage: epequil [-h] [-help_cmdfile] [-d #] [epequil.inp]"
         <<  endl;
    cout << "    -h           help" << endl;
    cout << "    -d           #   : level of debug printing" << endl;
    cout << "  epequil.inp    : command file" << endl;
    cout << "                     : (if missing, assume epequil.inp)" 
	 << endl;
    cout << endl;
}



void epequil_updateOne(ThermoPhase &p_ref, EPEQUIL_INPUT *pi) {
  int i, k;
  int nsp = p_ref.nSpecies();
  /*
   * The correct mole fractions are located in p_ref. Transfer these
   * to the EPEQUIL_INPUT object
   */
  p_ref.getMoleFractions(pi->spMf);
  /*
   * The correct temperature and pressure are in p_ref. Transfer these to the
   * EPEQUIL_INPUT object
   */
  pi->T = p_ref.temperature();
  pi->Pres = p_ref.pressure();
  /*
   * Search for the element with the most abundance
   */
  int eMost = 0;
  double eVal = 0.0;
  double sVal = 0.0;
 
  for (i = 0; i < pi->ne; i++) {
    sVal = 0.0;
    for (k = 0; k < nsp; k++) {
      double natoms = p_ref.nAtoms(k, i);
      if (natoms * pi->spMf[k] > sVal) {
	sVal = natoms * pi->spMf[k];
      }
    }
    if (sVal * pi->elementMoles[i] > eVal) {
      eMost = i;
      eVal = sVal * pi->elementMoles[i];
    }
  }
  /*
   * Deal with the degenerate case of zero moles
   */
  if (eVal <= 0.0) {
    pi->phaseMoles[0] = 0.0;
    for (i = 0; i < nsp; i++) {
      pi->spMoles[i] = 0.0;
    }
    pi->Vol = 0.0;
  } else {
    /*
     * Find the element abundance 
     */
    double baseVal = 0.0;
    for (k = 0; k < nsp; k++) {
      double natoms = p_ref.nAtoms(k, eMost);
      baseVal += natoms * pi->spMf[k];
    }
    if (baseVal <= 0.0) {
      printf("shouldn't be here\n");
      exit(-1);
    }
    double factor = pi->elementMoles[eMost] / baseVal;

    double sum = 0.0;
    for (k = 0; k < nsp; k++) {
      pi->spMoles[k] = factor * pi->spMf[k];
      sum += pi->spMoles[k];
    }
    pi->phaseMoles[0] = sum;
  

    /*
     * Check the element abundance calculation
     */
    double elAbundSum = 0.0;
    for (i = 0; i < pi->ne; i++) {
      elAbundSum += pi->elementMoles[i];
    }
    for (i = 0; i < pi->ne; i++) {
      double elAbundCalc = 0.0;
      for (k = 0; k < nsp; k++) {
	double natoms = p_ref.nAtoms(k, i);
	elAbundCalc += natoms * pi->spMoles[k];
      }
      if (fabs(elAbundCalc - pi->elementMoles[i]) / elAbundSum > 1.0E-6) {
	printf("Error in the calculation!\n");
	exit(-1);
      }
    }
  }
  /*
   * Store the chemical potentials
   */
  p_ref.getChemPotentials(pi->spChemPot);
  /*
   * Store the partial molar volumes
   */
  p_ref.getPartialMolarVolumes(pi->VolPM);
  /*
   * Calculate the total volume
   */
  pi->Vol = 0.0;
  for (k = 0; k < nsp; k++) {
    pi->Vol += pi->spMoles[k] * pi->VolPM[k];
  }

  /*
   * Now update the MultiPhase object
   */
  MultiPhase *mp = pi->m_mp;
  mp->setPhaseMoleFractions(0, pi->spMf);
  mp->setPhaseMoles(0, pi->phaseMoles[0]);
  mp->setTemperature(pi->T);
  mp->setPressure(pi->Pres);

}




int main(int argc, char **argv)
{
   int i;
   FILE *inputFP = stdin;
   EPEQUIL_INPUT *prob_input = new EPEQUIL_INPUT();
   int retn;
   string commandFile = "epequil.inp";

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
#ifdef DEBUG
		   epequil_debug_print_lvl = lvl;
#endif
		 }
	       }  
	     }
	   } else {
	     printUsage();
	     exit(1);
	   }
	 }
       } else if (commandFile == "") {
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
     retn = epequil_input(prob_input, commandFile);
     if (retn == -1) {
       printf("exiting with error\n");
       exit(-1);
     }
 
     /*
      * prep the input
      */
     int res = epequil_prep(prob_input);
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
     MultiPhase *mp =  prob_input->m_mp;

     /*
      * Find the first Thermophase
      */
     if (mp->nPhases() > 1) {
       printf("ERROR: epequil currently limited to 1 phase\n");
       exit(-1);
     }
     ThermoPhase &p_ref = mp->phase(0);

     Zuzax::ChemEquil ep;

     //beginLogGroup("epequil", 10);
#ifdef DEBUG_HKM_EPEQUIL
     debug_prnt_lvl = 0;
#endif

     ep.equilibrate(p_ref, "TP");
     //endLogGroup("epequil");

    
     epequil_updateOne(p_ref, prob_input);
 
     printf("\n Results from epequil:\n");
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
   
     delete prob_input;


     Zuzax::appdelete();

     return retn;

   } catch (ZuzaxError) {

     showErrors();
     return -1;
   }
}
