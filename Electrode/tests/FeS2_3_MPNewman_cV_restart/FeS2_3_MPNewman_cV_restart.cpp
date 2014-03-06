/*
 * $Id: FeS2_3_MPNewman_cV_restart.cpp 496 2013-01-07 21:15:37Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "cantera/equilibrium.h"
#include "cantera/thermo/MolalityVPSSTP.h"

#include "cantera/equil/vcs_MultiPhaseEquil.h"
#include "cantera/equil/vcs_prob.h"
#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/thermo/IonsFromNeutralVPSSTP.h"
#include "cantera/numerics/ResidEval.h"
#include "cantera/numerics/RootFind.h"
#include "cantera/numerics/NonlinearSolver.h"

#include "Electrode_input.h"
#include "Electrode.h"
#include "Electrode_MP_RxnExtent.h"
#include "ExtraGlobalRxn.h"
#include "RxnMolChange.h"
#include "BlockEntry.h"

#include <stdio.h>
#include <sstream>
#include <iomanip>

using namespace std;
using namespace Cantera;
using namespace VCSnonideal;

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

//=====================================================================================================
// void setPlateau(Electrode_MP_RxnExtent *eC) {

//   int crap = eC->plateauList_.size();
//   for (int i = 0; i < (int) crap; i++) {
//     delete eC->plateauList_[i];
//   }
//   eC->plateauList_.clear();

//   int isurf = 0;
//   int numInnerPhases = 1;
//   int numOuterPhases = 1;
//   Cantera::Plateau *ppo = new Plateau(isurf,numInnerPhases,numOuterPhases);
//   eC->plateauList_.push_back(ppo);
  
//   // set phase indices

//   ppo->phaseNameInnerSolidPhaseList_[0] = "FeS2(S)";
//   ppo->phaseNameOuterSolidPhaseList_[0] = "Li3Fe2S4(S)";
//   ppo->phaseIndexInnerSolidPhaseList_[0] = 2;
//   ppo->phaseIndexOuterSolidPhaseList_[0] = 3;

//   isurf = 1;
//   numInnerPhases = 1;
//   numOuterPhases = 2;
//   ppo = new Plateau(isurf, numInnerPhases, numOuterPhases);
//   eC->plateauList_.push_back(ppo);
  
//   // set phase indices
//   ppo->phaseIndexInnerSolidPhaseList_[0] = 3;
//   ppo->phaseIndexOuterSolidPhaseList_[0] = 4;
//   ppo->phaseIndexOuterSolidPhaseList_[1] = 5;

//   ppo->phaseNameInnerSolidPhaseList_[0] = "Li3Fe2S4(S)";
//   ppo->phaseNameOuterSolidPhaseList_[0] = "Li[2+x]Fe[1-x]S2(S)";
//   ppo->phaseNameOuterSolidPhaseList_[1] = "Fe[1-x]S(S)";


//   isurf = 2;
//   numInnerPhases = 2;
//   numOuterPhases = 2;
//   ppo = new Plateau(isurf, numInnerPhases, numOuterPhases);
//   eC->plateauList_.push_back(ppo);
  
//   // set phase indices
//   ppo->phaseIndexInnerSolidPhaseList_[0] = 4;
//   ppo->phaseIndexInnerSolidPhaseList_[1] = 5;
//   ppo->phaseIndexOuterSolidPhaseList_[0] = 6;
//   ppo->phaseIndexOuterSolidPhaseList_[1] = 7;

//   ppo->phaseNameInnerSolidPhaseList_[0] = "Li[2+x]Fe[1-x]S2(S)";
//   ppo->phaseNameInnerSolidPhaseList_[1] = "Fe[1-x]S(S)";
//   ppo->phaseNameOuterSolidPhaseList_[0] = "Li2S(S)";
//   ppo->phaseNameOuterSolidPhaseList_[1] = "Fe(S)";

// }

//=====================================================================================================


int main(int argc, char **argv)
{
  int retn = 0;
  //bool doCathode = false;
  string commandFileNet = "cell.inp";
  string commandFileA = "electrodeAnode.inp";
  //  string commandFileC = "electrodeCathode.inp";
  string commandFileC = "cathode.inp";
  // printed usage

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
      } else if (commandFileNet == "") {
	commandFileNet = tok;
      } else {
	printUsage();
	exit(1);
      }
    }
  }


  try {

    /*
     * Go get the Cell problem description from the input file
     */
    //  retn = cell_input(commandFileNet);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }
  
    Cantera::Electrode_MP_RxnExtent *electrodeC  = new Cantera::Electrode_MP_RxnExtent();

    ELECTRODE_KEY_INPUT *electrodeC_input = new ELECTRODE_KEY_INPUT();
    
    std::string commandFileC = "cathode.inp";


	
    // Initialize a block input structure for the command file

    BEInput::BlockEntry *cfC = new BEInput::BlockEntry("command_file");

 
    // Go get the problem description from the input file
    electrodeC_input->printLvl_ = 5;
    retn = electrodeC_input->electrode_input(commandFileC, cfC);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

    retn =  electrodeC->electrode_input_child(&electrodeC_input);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

  
    retn = electrodeC->electrode_model_create(electrodeC_input);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }  
    retn = electrodeC->electrode_stateSave_create();
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }
    retn = 0;



    double deltaT = 0.1;
    double Tinitial = 0.0;
    double Tfinal = 0.0;

    electrodeC->setVoltages(1.95, 0.0);
    double molNum[10];

    electrodeC->setPhaseExistenceForReactingSurfaces(true);
    electrodeC->setVoltages(1.60, 0.0);

    double oc = electrodeC->openCircuitVoltageSSRxn(1);
    oc = electrodeC->openCircuitVoltage(1);
    printf("oc[0] = %g\n", oc);


    int nT = 50;
    deltaT = 2.0E-1;

    electrodeC->printElectrode();
    electrodeC->setPrintLevel(4);
    electrodeC->setDeltaTSubcycle(0.01);
    electrodeC->detailedResidPrintFlag_ = 10;
    electrodeC->enableExtraPrinting_ = 0;
  
    /*
     * load the saved solution
     */
    XML_Node* xSavedSoln = get_XML_File("soln_0_0_orig.xml");

    XML_Node *xGTSI = electrodeC->selectGlobalTimeStepIncrement(xSavedSoln, 1);

    electrodeC->loadGlobalTimeStepTFinalState(xGTSI);
    Tfinal = electrodeC->getFinalTime();

    electrodeC->printElectrode();
    electrodeC->setPrintLevel(4);
    electrodeC->detailedResidPrintFlag_ = 10;
    electrodeC->enableExtraPrinting_ = 0;
    electrodeC->printCSVLvl_ = 9;

    for (int itimes = 0; itimes < nT; itimes++) {
      Tinitial = Tfinal;

      electrodeC->resetStartingCondition(Tinitial);
      
      Tfinal = Tinitial + deltaT;

      electrodeC->integrate(deltaT);
      electrodeC->getMoleNumSpecies(molNum);
      doublereal net[12];
      double amps = electrodeC->getIntegratedProductionRatesCurrent(net);
 
      cout << setw(15) << Tfinal << setw(15) << amps << endl;
      electrodeC->printElectrode();
      electrodeC->writeSolutionTimeIncrement();
  
    }


    delete cfC;
    delete electrodeC_input;
    delete electrodeC;
    Cantera::appdelete();

    return retn;

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
} 
