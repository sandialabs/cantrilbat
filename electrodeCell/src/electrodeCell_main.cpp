/**
 *  @file electrodeCell_main.cpp
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "zuzax/equilibrium.h"
#include "zuzax/thermo/MolalityVPSSTP.h"

#include "zuzax/equil/vcs_MultiPhaseEquil.h"
#include "zuzax/equil/vcs_prob.h"
#include "zuzax/equil/vcs_solve.h"
#include "zuzax/equil/vcs_VolPhase.h"
#include "zuzax/equil/vcs_internal.h"

#include "Electrode_input.h"
#include "Electrode.h"
#include "Electrode_Factory.h"
#include "Electrode_Equilibrium.h"
#include "electrodeCell_prep.h"
#include "electrodeCell_kin.h"
#include "cell_input.h"
#include "zuzax/kinetics/RxnMolChange.h"

//#include <cstdio>

using namespace std;
using namespace Zuzax;
using namespace mdpUtil;


// a lvl of one prints out the .csv file
int mpequil_debug_print_lvl = 1;
int VCS_Debug_Print_Lvl = 3;

void printUsage() {
    cout << "usage: electrodeCell [-h] [-help_cmdfile] [-d #] [mpequil.inp]"
         <<  endl;
    cout << "    -h           help" << endl;
    cout << "    -d           #   : level of debug printing" << endl;
    cout << "  electrodeCell.inp    : command file" << endl;
    cout << "                     : (if missing, assume mpequil.inp)" 
	 << endl;
    cout << endl;
}


//! This routine converts the MPEQUIL_INPUT problem into a VCS_PROB
//! structure

int mpequil_convert(Zuzax::Electrode *electrode, vcs_nonideal::VCS_PROB *vprob,  Zuzax::MP_EquilStatic *mix) {

  /*
   *     Extract the current state information
   *     from the MultiPhase object and
   *     Transfer it to VCS_PROB object.
   */
  //int res = vcs_Cantera_to_vprob(mix, vprob);
  int res = vprob->input(mix);
  if (res != 0) {
    printf("problems\n");
  }

  // Set the estimation technique
  vprob->setInitialEstimateMethod(-1);


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

  return 0;
}


//! This routine takes the problem as specified in the MPEQUIL_INPUT structure
//! and solves it to equilibrium
/*!
 *
 */

int mpequil_equilibrate(Zuzax::Electrode *electrode, int estimateInit, int printFlag) {

  Electrode_Equilibrium *ee_equil = new Electrode_Equilibrium(electrode);

  MP_EquilStatic *mix = ee_equil->MultiPhase_Obj();
  /*
   * Malloc a new vcs problem strcutre
   */
  int ne = electrode->nElements();
  int nVphase = electrode->nVolPhases();
  int nVspecies = electrode->nVolSpecies();
  vcs_nonideal::VCS_PROB *vprob = new vcs_nonideal::VCS_PROB(nVspecies, ne, nVphase);
#ifdef DEBUG
  vprob->setDebugPrintLvl(VCS_Debug_Print_Lvl);
#endif
  /*
   * vcs problem input to Zuzax conversion
   */
  int res = mpequil_convert(electrode, vprob, mix);
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
  int ipr = std::max(0, printFlag - 1);
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
  for (int ip = 0; ip < (int) vprob->nPhases(); ip++) {
    double phaseMole = 0.0;
    ThermoPhase &tref = mix->phase(ip);
    int nspPhase = tref.nSpecies();
    for (int k = 0; k < nspPhase; k++, kGlob++) {
      phaseMole += vprob->speciesMoles(kGlob);
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
    printf("Pressure    = %g Pa\n", vprob->pressure());
    printf("\n");
    printf("----------------------------------------"
	   "---------------------\n");
    printf(" Name            KMole_Number");
    printf("(kmol)");
    printf("  Mole_Fraction     Chem_Potential");
    if (vprob->vcsUnitsFormat() == VCS_UNITS_KCALMOL) 
      printf(" (kcal/mol)\n");
    else if (vprob->vcsUnitsFormat() == VCS_UNITS_UNITLESS) 
      printf(" (Dimensionless)\n");
    else if (vprob->vcsUnitsFormat() == VCS_UNITS_KJMOL) 
      printf(" (kJ/mol)\n");
    else if (vprob->vcsUnitsFormat() == VCS_UNITS_KELVIN) 
      printf(" (Kelvin)\n");
    else if (vprob->vcsUnitsFormat() == VCS_UNITS_MKS) 
      printf(" (J/kmol)\n");
    printf("--------------------------------------------------"
	   "-----------\n");
    const std::vector<doublevalue>& vg = vprob->vec_gibbsSpecies();
    for (size_t i = 0; i < vprob->nSpecies(); i++) {
      printf("%-12s", vprob->speciesName(i).c_str());
      if (vprob->speciesUnknownType(i) == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	printf("  %15.3e %15.3e  ", 0.0, vprob->mfVectorC()[i]);
	printf("%15.3e\n", vg[i]);
      } else {
	printf("  %15.3e   %15.3e  ", vprob->speciesMoles(i), vprob->mfVectorC()[i]);
	if (vprob->speciesMoles(i) <= 0.0) {
	  size_t iph = vprob->globalPhaseIndexFromGlobalSpeciesIndex(i);
	  vcs_nonideal::vcs_VolPhase *VPhase = vprob->phaseController(iph);
	  if (VPhase->nSpecies() > 1) {
	    printf("     -1.000e+300\n");
	  } else {
	    printf("%15.3e\n", vg[i]);
	  }
	} else {
	  printf("%15.3e\n", vg[i]);
	}
      }
    }
    printf("-------------------------------------------------------------\n"); 

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
//======================================================================================================================


int main(int argc, char **argv)
{
  int i;
  FILE* inputFP = stdin;

  int retn;
  //bool doCathode = false;
  std::string commandFileNet = "cell.inp";
  std::string commandFileA = "electrodeAnode.inp";
  std::string commandFileC = "electrodeCathode.inp";
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
      } else if (commandFileNet == "") {
	commandFileNet = tok;
      } else {
	printUsage();
	exit(1);
      }
    }
  }

  if (commandFileNet != "") {
    inputFP = fopen(commandFileNet.c_str(), "r");
    if (!inputFP) {
      printf("Can't open file, %s, bailing\n", commandFileNet.c_str());
      return -1;
    }
  }


  try {

    /*
     * Go get the Cell problem description from the input file
     */
    retn = cell_input(commandFileNet);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

    // ------------------------------------------------------------------------------------

    commandFileA = CellO.anodeElectrodeFileName_;

    // Initialize a block input structure for the command file
    BEInput::BlockEntry *cfA = new BEInput::BlockEntry("command_file");

    //  Go get the name of the Electrode model first, then delete the stub reader
    ELECTRODE_KEY_INPUT *electrodeA_input = new ELECTRODE_KEY_INPUT();
    retn = electrodeA_input->electrode_input(commandFileA, cfA);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }
    std::string aen = electrodeA_input->electrodeModelName;
    delete electrodeA_input;

    //  Use the ElectrodeModelName value as input to the electrode factory to create the electrode
    Zuzax::Electrode* electrodeA = newElectrodeObject(aen);
    if (!electrodeA) {
        throw ZuzaxError("electrodeCell_main()", "newElectrodeObject failed for model");
    }
    electrodeA_input = newElectrodeKeyInputObject(aen);
    if (!electrodeA_input) {
        throw ZuzaxError("electrodeCell_main()", "newElectrodeKeyInputObject() failed for model %s", aen.c_str());
    }
    
    //  Parse the complete child input file
    retn = electrodeA_input->electrode_input_child(commandFileA, cfA);
    if (retn == -1) {
        throw ZuzaxError("electrodeCell_main()", "electrode_input_child() method failed on input file %s", commandFileA.c_str());
    }
    delete cfA;

    // Now create the innards of the Electrode object
    retn = electrodeA->electrode_model_create(electrodeA_input);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

    // Read in the initial conditions
    retn = electrodeA->setInitialConditions(electrodeA_input);
    if (retn == -1) {
        throw ZuzaxError("electrodeCell_main()", "Electrode::setInitialConditions() for anode failed");
    }

    // ------------------------------------------------------------------------------------

    commandFileC = CellO.cathodeElectrodeFileName_;
    
    // Initialize a block input structure for the command file
    BEInput::BlockEntry *cfC = new BEInput::BlockEntry("command_file");
    ELECTRODE_KEY_INPUT *electrodeC_input = new ELECTRODE_KEY_INPUT();
    retn = electrodeC_input->electrode_input(commandFileC, cfC);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }
    /*
     *  Use the ElectrodeModelName value as input to the electrode factory to create the electrode
     */
    std::string en = electrodeC_input->electrodeModelName;
    delete electrodeC_input;
    Zuzax::Electrode* electrodeC = newElectrodeObject(en);
    electrodeC_input = newElectrodeKeyInputObject(en);
    if (!electrodeC_input) {
        throw ZuzaxError("electrodeCell_main()", "newElectrodeKeyInputObject() failed for model %s", en.c_str());
    }

    /*
     *  Parse the complete child input file
     */
    retn = electrodeC_input->electrode_input_child(commandFileC, cfC);
    if (retn == -1) {
        throw ZuzaxError("electrodeCell_main()", "electrode_input_child() method failed while reading file, %s",
                         commandFileC.c_str());
    }

    // Now construct the Electrode object's innards based on the input file
    retn = electrodeC->electrode_model_create(electrodeC_input);
    if (retn == -1) {
        throw ZuzaxError("electrodeCell_main()::electrode_model_create()", 
                         "Error initializing the cathode electrode object");
    }

    // Now read and implement the initial conditions
    retn = electrodeC->setInitialConditions(electrodeC_input);
    if (retn == -1) {
        throw ZuzaxError("electrodeCell_main()", "Electrode::setInitialConditions() for cathode failed");
    }
    delete cfC;

    // ---------------------------------------------------------------------------------------------------------------

    /*
     * prep the input
     */
    int res = electrode_prep(electrodeA);
    if (res != 0) {
      printf("problems\n");
    }

    /*
     * Print out the problem specification
     */
    //vprob->prob_report(2);

    /*
     * Call the thermo Program for the first electrode
     */
    Electrode_Equilibrium* ee_equilA = new Electrode_Equilibrium(electrodeA);
    ee_equilA->setupEquilibriumProblem();
    MP_EquilStatic *mpA =  ee_equilA->MultiPhase_Obj();

    int printFlag = mpequil_debug_print_lvl;
    int estimateInit = 0;
 
    printFlag = 9;
    vcs_equilibrate_1(*mpA, TP, estimateInit, printFlag);

    ee_equilA->uploadMP();


         
    double AopenCircuitVoltageEst = electrodeA->voltage();

  
    printf("\n Results from mpequil:\n");
    printf("\n");
    printf("Temperature = %g Kelvin\n",  electrodeA->temperature());
    printf("Pressure    = %g ", electrodeA->pressure());
    printf("Pa\n");
 
    printf("\n");
    printf("----------------------------------------"
	   "---------------------\n");
    printf(" Name             Mole_Number     Mole_Fraction     Chem_Potential");
    printf(" (J/kmol)\n");
    printf("--------------------------------------------------"
	   "-----------\n");
    string spname;
   
    for (i = 0; i < (int) mpA->nSpecies(); i++) {
      spname = mpA->speciesName(i);
      printf("%-12s", spname.c_str());
      printf("  %15.6e %15.6e  %15.6e\n", electrodeA->moleNumSpecies(i),
	     electrodeA->moleFraction(i), electrodeA->speciesElectrochemPotential(i));
    }

 
    printf("------------------------------------------"
	   "-------------------\n");   
   
    for (int iph = 0; iph < (int) mpA->nPhases(); iph++) {
      ThermoPhase *tp = &(mpA->phase(iph));
      string phname = tp->id();
      printf("%-12s", phname.c_str());
      printf("  %15.6e %15.6e\n", electrodeA->phaseMoles(iph),
	     electrodeA->phaseElectricPotential(iph));
    }

    double IcurrNet = 0.0;

    // Set the Interfacial Voltage
    // InterfaceKinetics *iKA = electrodeA->m_rSurDomain;
    InterfaceKinetics *iKA = electrodeA->currOuterReactingSurface();
    //int nReactionsA = iKA->nReactions();
    //
    // Read in the number of extra global pathways
    //
    int numExtraGlobalRxnsA = electrodeA->numExtraGlobalRxnPathways();

    for (int iextra = 0; iextra <  numExtraGlobalRxnsA; iextra++) {
      //struct EGRInput * egr_ptr = electrodeA->m_EGRList[iextra];

     // electrodeA->addExtraGlobalRxn(egr_ptr);
      /*
      Zuzax::ExtraGlobalRxn *egr = new ExtraGlobalRxn(iKA);
      double *RxnVector  = new double[nReactionsA];
      for (int i = 0; i < nReactionsA; i++) {
	RxnVector[i] = 0.0;
      }
      for (int iErxn = 0; iErxn < egr_ptr->m_numElemReactions; iErxn++) {
	struct ERSSpec * ers_ptr = egr_ptr->m_ERSList[iErxn];
	RxnVector[ers_ptr->m_reactionIndex] = ers_ptr->m_reactionMultiplier;
      }
      egr->setupElemRxnVector(RxnVector, egr_ptr->m_SS_KinSpeciesKindex);

      RxnMolChange *rmcEGR = new RxnMolChange(iKA, egr);

      electrodeA->m_egr.push_back(egr);
      electrodeA->m_rmcEGR.push_back(rmcEGR);
      */
      Zuzax::RxnMolChange *rmcEGR = electrodeA->rxnMolChangesEGR(iextra); 
      Zuzax::ExtraGlobalRxn *egr = electrodeA->extraGlobalRxnPathway(iextra);
      if ((rmcEGR)->m_ChargeTransferInRxn != 0.0) {
	//processGERCurrentVsPotTable(rmcEGR, pl, 0, TT,
	//		    *iK, *egr, *rts);

	IcurrNet = processGERCurrent(rmcEGR, electrodeA, iextra, *iKA, *egr);
      }

 
      //delete [] RxnVector;

    }


    // Equilibrate the cathode




    /*
     * Call the thermo Program for the first electrode
     */
    Electrode_Equilibrium* ee_equilC = new Electrode_Equilibrium(electrodeC);
    ee_equilC->setupEquilibriumProblem();
    MP_EquilStatic *mpC =  ee_equilC->MultiPhase_Obj();

    printFlag = mpequil_debug_print_lvl;
    estimateInit = 0;
 
    printFlag = 9;
    vcs_equilibrate_1(*mpC, TP, estimateInit, printFlag);

    ee_equilC->uploadMP();

    double CopenCircuiteVoltageEst = electrodeC->voltage();


   // Process the ExtraGlobalRxns in the Cathode
    InterfaceKinetics *iKC = electrodeC->currOuterReactingSurface();
 
    int numExtraGlobalRxnsC = electrodeC->numExtraGlobalRxnPathways();

    for (int iextra = 0; iextra < numExtraGlobalRxnsC; iextra++) {
   //   struct EGRInput * egr_ptr = electrodeC->m_EGRList[iextra];

    //  electrodeC->addExtraGlobalRxn(egr_ptr);
      /*
      Zuzax::ExtraGlobalRxn *egr = new ExtraGlobalRxn(iKC);
      double *RxnVector  = new double[nReactionsC];
      for (int i = 0; i < nReactionsC; i++) {
	RxnVector[i] = 0.0;
      }
      for (int iErxn = 0; iErxn < egr_ptr->m_numElemReactions; iErxn++) {
	struct ERSSpec * ers_ptr = egr_ptr->m_ERSList[iErxn];
	RxnVector[ers_ptr->m_reactionIndex] = ers_ptr->m_reactionMultiplier;
      }
      egr->setupElemRxnVector(RxnVector, egr_ptr->m_SS_KinSpeciesKindex);

      RxnMolChange *rmcEGR = new RxnMolChange(iKC, egr);

      electrodeC->m_egr.push_back(egr);
      electrodeC->m_rmcEGR.push_back(rmcEGR);
      */
      Zuzax::ExtraGlobalRxn *egr = electrodeC->extraGlobalRxnPathway(iextra);
      Zuzax::RxnMolChange *rmcEGR = electrodeC->rxnMolChangesEGR(iextra); 

      if ((rmcEGR)->m_ChargeTransferInRxn != 0.0) {
	//processGERCurrentVsPotTable(rmcEGR, pl, 0, TT,
	//		    *iK, *egr, *rts);

	IcurrNet = processGERCurrent(rmcEGR, electrodeC, iextra, *iKC, *egr);
      }

 
      //delete [] RxnVector;

    }

    printf("IcurrNet = %g\n", IcurrNet);


    for (int i = 0; i < 11; i++) {
      double deltaV = AopenCircuitVoltageEst  + i * 0.02;
      //electrodeA->setDeltaVoltage(deltaV); 
      electrodeA->setVoltages(deltaV, 0.0);
      Zuzax::RxnMolChange * rmcEGR = electrodeA->rxnMolChangesEGR(0);
      Zuzax::ExtraGlobalRxn * egr = electrodeA->extraGlobalRxnPathway(0);
      IcurrNet = processGERCurrent(rmcEGR, electrodeA, 0, *iKA, 
				   *(egr));

      printf("Voltage = %g, Icurr = %g\n", deltaV, IcurrNet);
    }
    double Vlow  = AopenCircuitVoltageEst  - 1.0 ;
    double Vhigh = AopenCircuitVoltageEst  + 1.0;
    double deltaV = findV(electrodeA, 0.01,  Vlow, Vhigh, 9, 1.0E-7, 1000); 
					
    printf("delta V = %g for I = 0.01\n", deltaV);
    double totalV;
    double deltaVA = AopenCircuitVoltageEst;
    double deltaVC = CopenCircuiteVoltageEst;
    for (int i = 0; i < 11; i++) {
      double Itarget = 0.0 + i * 0.001;
      Vlow = deltaVA - 1.0;
      Vhigh = deltaVA + 1.0;
      deltaVA = findV(electrodeA, Itarget, Vlow, Vhigh, 9, 1.0E-7, 1000); 
      Vlow = deltaVC - 1.0;
      Vhigh = deltaVC + 1.0;
      deltaVC = findV(electrodeC, -Itarget, Vlow, Vhigh, 9, 1.0E-7, 1000); 
      totalV = deltaVC - deltaVA;
      
      printf("Itarget = %g, deltaVA = %g, deltaVC = %g, totalV = %g\n",
	     Itarget, deltaVA, deltaVC, totalV);
    }
    
    

    delete electrodeA_input;


    Zuzax::appdelete();

    return retn;

  } catch (ZuzaxError) {

    showErrors();
    return -1;
  }
} 
