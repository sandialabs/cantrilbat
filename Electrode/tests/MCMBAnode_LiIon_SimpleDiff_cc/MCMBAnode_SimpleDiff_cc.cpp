/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "zuzax/equilibrium.h"
#include "zuzax/thermo/MolalityVPSSTP.h"

#include "zuzax/equil/vcs_MultiPhaseEquil.h"
#include "zuzax/equil/vcs_solve.h"
#include "zuzax/equil/vcs_VolPhase.h"
#include "zuzax/numerics/ResidEval.h"
#include "zuzax/numerics/NonlinearSolver_JAC.h"

#include "Electrode_input.h"
#include "Electrode_SimpleDiff.h"
#include "Electrode_RadialDiffRegions.h"  

#include <sstream>
#include <iomanip>

using namespace std;
using namespace Zuzax;

// a lvl of one prints out the .csv file
int mpequil_debug_print_lvl = 1;

void printUsage() {
    cout << "usage: MCMBAnode_SimpleDiff [-h] [-help_cmdfile] [-d #] [anode.inp]"
         <<  endl;
    cout << "    -h               : Prints this help" << endl;
    cout << "    -help_cmdfile    : Prints a list of block commands understood by this parser - add anode.inp for more information" << endl;
    cout << "   -d #              : Level of debug printing" << endl;
    cout << "   anode.inp         : Command file (if missing, assume anode.inp)" << endl;
    cout << endl;
}

//=====================================================================================================


int main(int argc, char **argv)
{
  int retn = 0;
  //bool doCathode = false;
  std::string commandFileNet = "anode.inp";
  std::string commandFileA = "anode.inp";

  bool printInputFormat = false; // print cmdfile.txt format
  // printed usage

  //VCSnonideal::vcs_timing_print_lvl = 0;
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
  
    Electrode_SimpleDiff *electrodeA  = new Electrode_SimpleDiff();

    ELECTRODE_RadialDiffRegions_KEY_INPUT *electrodeA_input = new ELECTRODE_RadialDiffRegions_KEY_INPUT();
    
    std::string commandFileA = commandFileNet;
	
    // Initialize a block input structure for the command file

    BEInput::BlockEntry *cfC = new BEInput::BlockEntry("command_file");

 
    // Go get the problem description from the input file

    retn = electrodeA_input->electrode_input_child(commandFileA, cfC);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

    if (printInputFormat) {
      cfC->print_usage();
      exit(-1);
    }

    retn = electrodeA->electrode_model_create(electrodeA_input);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }
    retn = electrodeA->setInitialConditions(electrodeA_input);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }
    retn = electrodeA->electrode_stateSave_create();
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }
    double deltaT = 1.0E-2;
    double Tinitial = 0.0;
    double Tfinal = 0.0;

    electrodeA->setVoltages(0.0, -0.2);

    electrodeA->setPhaseExistenceForReactingSurfaces(true);
    electrodeA->setVoltages(0.0, -0.2);
    double oc = electrodeA->openCircuitVoltageSSRxn(0);
    oc = electrodeA->openCircuitVoltage(0);
    printf("oc[0] = %g\n", oc);
    int nT = 100;
    electrodeA->printCSVLvl_ = 1;

    electrodeA->printElectrode();
    electrodeA->setDeltaTSubcycle(0.001);

    remove("soln.xml");

    // Extra printing for debugging
    electrodeA->enableExtraPrinting_ = true;
    electrodeA->detailedResidPrintFlag_ = 12;
    // Normal printing for testing
    electrodeA->enableExtraPrinting_ = false; 
    electrodeA->detailedResidPrintFlag_ = 0;
    electrodeA->printLvl_ = 4;

    electrodeA->fixCapacityBalances_final();

    double net[50];
    double ampsGoal = 5.0;
    double deltaTgoal = deltaT;
    double coul = 0.0;

    FILE * fp = fopen("outputTable.csv","w");
    fprintf(fp, "Constant Current Curves \n");
    fprintf(fp, "amps = %g  \n", ampsGoal);
    fprintf(fp, "deltaT = %g  \n", deltaT);
    fprintf(fp, "\n");
    fprintf(fp, "  Tfinal      ,      Coul     ,     Ah   , Volts   , amps ,   globalError , numSubcyles \n");

    int maxSV;  
    for (int itimes = 0; itimes < nT; itimes++) {
      Tinitial = Tfinal;

      electrodeA->resetStartingCondition(Tinitial);
      double amps = ampsGoal;
      deltaT = deltaTgoal;

      double volts = electrodeA->integrateConstantCurrent(amps, deltaT, 1.0, -1.0, 1.0E-4);
      Tfinal = Tinitial + deltaT;

      amps = electrodeA->getIntegratedProductionRatesCurrent(net);
      double errGlob = electrodeA->reportStateVariableIntegrationError(maxSV);
      int numSubcycles = electrodeA->numSubcycles();

      coul  += amps * deltaT;
 
      fprintf(fp, " %12.6E ,  %12.6E , %12.6E , %12.6E , %12.6E , %10.2E , %5d \n", Tfinal, coul, coul/3600. , volts, amps , errGlob, numSubcycles );

      electrodeA->printElectrode();
      electrodeA->writeSolutionTimeIncrement();
    }
    delete cfC;
    delete electrodeA_input;
    delete electrodeA;
    appdelete();

    return 0;

  } catch (ZuzaxError) {

    showErrors();
    return -1;
  }
} 
