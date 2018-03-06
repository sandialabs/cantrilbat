/*
 * $Id: MCMBAnode_CSTR_2.cpp 274 2013-03-26 17:04:11Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "cantera/equilibrium.h"
#include "cantera/thermo/MolalityVPSSTP.h"

#include "cantera/numerics/RootFind.h"
#include "cantera/numerics/NonlinearSolver.h"

#include "Electrode_input.h"
#include "Electrode_CSTR_MCMBAnode.h"

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


int main(int argc, char **argv)
{
  int retn = 0;
  //bool doCathode = false;
  string commandFileNet = "cell.inp";

  string commandFileA = "anode.inp";
  // printed usage

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
  
    Electrode_CSTR_MCMBAnode *electrodeA  = new Electrode_CSTR_MCMBAnode();

    ELECTRODE_KEY_INPUT *electrodeA_input = new ELECTRODE_KEY_INPUT();
    
    std::string commandFileA = "anode.inp";

	
    // Initialize a block input structure for the command file

    BEInput::BlockEntry *cfC = new BEInput::BlockEntry("command_file");

 
    // Go get the problem description from the input file
    electrodeA_input->printLvl_ = 5;
    retn = electrodeA_input->electrode_input(commandFileA, cfC);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

    retn = electrodeA->electrode_input_child(&electrodeA_input);
    if (retn == -1) {
      printf("exiting with error\n");
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
    double deltaT = 0.1;
    double Tinitial = 0.0;
    double Tfinal = 0.0;

    electrodeA->setVoltages(0.0, -0.2);
    electrodeA->setPrintLevel(4);
    double molNum[10];

    electrodeA->setPhaseExistenceForReactingSurfaces(true);
    electrodeA->setVoltages(0.0, -0.2);
    double oc = electrodeA->openCircuitVoltageSSRxn(0);
    oc = electrodeA->openCircuitVoltage(0);
    printf("oc[0] = %g\n", oc);
    int nT = 50;
    deltaT = 2.0E-1;
    electrodeA->printCSVLvl_ = 3;

    electrodeA->printElectrode();
    electrodeA->setDeltaTSubcycle(0.01);

    (void) remove("soln.xml");
  
    for (int itimes = 0; itimes < nT; itimes++) {
      Tinitial = Tfinal;

      electrodeA->resetStartingCondition(Tinitial);
      Tfinal = Tinitial + deltaT;
      electrodeA->integrate(deltaT);
      electrodeA->getMoleNumSpecies(molNum);
      double net[12];
      double amps = electrodeA->getIntegratedProductionRatesCurrent(net);
 
      cout << setw(15) << Tfinal << setw(15) << amps << endl;
      electrodeA->printElectrode();
      electrodeA->writeSolutionTimeIncrement();
    }
    delete cfC;
    delete electrodeA_input;
    delete electrodeA;
    appdelete();

    return 0;

  } catch (CanteraError& ce) {

    showErrors();
    return -1;
  }
} 
