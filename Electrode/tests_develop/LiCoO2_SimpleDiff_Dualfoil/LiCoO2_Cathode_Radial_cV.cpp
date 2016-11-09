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

#include "cantera/equil/vcs_MultiPhaseEquil.h"
//#include "cantera/equil/vcs_prob.h"
//#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/thermo/IonsFromNeutralVPSSTP.h"
#include "cantera/numerics/ResidEval.h"
#include "cantera/numerics/NonlinearSolver.h"

#include "Electrode_input.h"
#include "Electrode_SimpleDiff.h"
#include "Electrode_RadialDiffRegions.h"  

#include <sstream>
#include <iomanip>

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

// a lvl of one prints out the .csv file
int mpequil_debug_print_lvl = 1;
int VCS_Debug_Print_Lvl = 3;

void printUsage() {
    cout << "usage: LiCoO2_Cathode_Radial_cV [-h] [-help_cmdfile] [-d #] [anode.inp]"
         <<  endl;
    cout << "    -h               : Prints this help" << endl;
    cout << "    -help_cmdfile    : Prints a list of block commands understood by this parser - add cathode.inp for more information" << endl;
    cout << "   -d #              : Level of debug printing" << endl;
    cout << "   cathode.inp         : Command file (if missing, assume cathode.inp)" << endl;
    cout << endl;
}

//=====================================================================================================


int main(int argc, char **argv)
{
  int retn = 0;
  //bool doCathode = false;
  std::string commandFileNet = "cathode.inp";
  std::string commandFileC = "cathode.inp";

  bool printInputFormat = false; // print cmdfile.txt format
  // printed usage

  //VCSnonideal::vcs_timing_print_lvl = 0;
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
	    printInputFormat = true;
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
      } else if (commandFileNet == "anode.inp") {
	commandFileNet = tok;
      } else {
	printUsage();
	exit(1);
      }
    }
  }

  try {

    /*
     *  Use this opportunity to set environment to get consistent printing for test cases
     *  Note, this has moved to electrode_model_create routine, generally.
     */ 
    Electrode::readEnvironmentalVariables();
  
    ZZCantera::Electrode_SimpleDiff *electrodeC  = new ZZCantera::Electrode_SimpleDiff();

    ELECTRODE_RadialDiffRegions_KEY_INPUT *electrodeC_input = new ELECTRODE_RadialDiffRegions_KEY_INPUT();
    
    std::string commandFileC = commandFileNet;
	
    // Initialize a block input structure for the command file

    BEInput::BlockEntry *cfC = new BEInput::BlockEntry("command_file");


    
 
    // Go get the problem description from the input file

    retn = electrodeC_input->electrode_input_child(commandFileC, cfC);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

    if (printInputFormat) {
      cfC->print_usage();
      exit(-1);
    }

    retn = electrodeC->electrode_model_create(electrodeC_input);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

    retn = electrodeC->setInitialConditions(electrodeC_input);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

    retn = electrodeC->electrode_stateSave_create();
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

    double Tinitial = 0.0;
    double Tfinal = 0.0;
    double molNum[10];
    double vvolts = electrodeC->openCircuitVoltageRxn(0, 0);
    printf(" volts = %g\n", vvolts);

    electrodeC->setPhaseExistenceForReactingSurfaces(true);
    electrodeC->setVoltages(3.7, 0.0);
    double oc = electrodeC->openCircuitVoltageSSRxn(0);
    oc = electrodeC->openCircuitVoltage(0);
    printf("oc[0] = %g\n", oc);
    int nT = 30;
    double deltaT = 1.0E-1;
    electrodeC->printCSVLvl_ = 3;

    double pmv[10];

    ThermoPhase *th = electrodeC->phasePtr("LiCoO2_Interstitials_cathode");
    th->getPartialMolarVolumes(pmv);
    
    electrodeC->setDeltaTSubcycle(2.5E-4);
    electrodeC->printElectrode();

    remove("soln.xml");

    electrodeC->setPrintLevel(1);
    //electrodeC->enableExtraPrinting_ = 10;
    //electrodeC->detailedResidPrintFlag_ = 10;
    //electrodeC->setMaxNumberSubCycles(40);


    nT = 50; 
    for (int itimes = 0; itimes < nT; itimes++) {
      Tinitial = Tfinal; 
      electrodeC->resetStartingCondition(Tinitial);
      int numSubIntegrations = electrodeC->integrate(deltaT);
      Tfinal = electrodeC->timeFinalFinal();
      electrodeC->getMoleNumSpecies(molNum);
      doublereal net[12];
      double amps = electrodeC->getIntegratedProductionRatesCurrent(net);
      cout << setw(15) << Tfinal << setw(15) << amps << numSubIntegrations << endl;
      electrodeC->printElectrode();
      electrodeC->writeSolutionTimeIncrement();
    }
    delete cfC;
    delete electrodeC_input;
    delete electrodeC;
    ZZCantera::appdelete();

    return 0;

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
} 
