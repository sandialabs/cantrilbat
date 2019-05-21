/*
 * $Id: LiCoO2_Cathode_3.cpp 202 2013-01-09 00:10:07Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "zuzax/equilibrium.h"

#include "zuzax/equil/vcs_VolPhase.h"
#include "zuzax/equil/vcs_internal.h"
#include "zuzax/thermo/IonsFromNeutralVPSSTP.h"
#include "zuzax/numerics/ResidEval.h"
#include "zuzax/numerics/NonlinearSolver.h"
#include "zuzax/numerics/DAE_Solver.h"

#include "Electrode_input.h"
#include "Electrode_CSTR_LiCoO2Cathode.h"
#include "GFCEO_Electrode.h"

#include <stdio.h>
#include <sstream>
#include <iomanip>

using namespace std;
using namespace Zuzax;
using namespace vcs_nonideal;
//using namespace mdpUtil;

// a lvl of one prints out the .csv file
int debug_print_lvl = 16;
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
  std::string commandFileA = "cathode.inp";

  vcs_nonideal::vcs_timing_print_lvl = 0;
  NonlinearSolver::s_TurnOffTiming = true;
  NonlinearSolver::s_print_NumJac = true;

  /*
   * Process the command line arguments
   */ 
  if (argc > 1) {
    std::string tok;
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
		  debug_print_lvl = lvl;
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
  
    Electrode_CSTR_LiCoO2Cathode *electrodeC = new Electrode_CSTR_LiCoO2Cathode();
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

    retn = electrodeC->electrode_input_child(&electrodeC_input);
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

    retn = electrodeC->setInitialConditions(electrodeC_input);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

    double deltaT = 0.1;
    double Tinitial = 0.0;
    double Tfinal = 0.0;

    double molNum[10];

    electrodeC->setPhaseExistenceForReactingSurfaces(true);
    electrodeC->setVoltages(3.2, 0.0);

    double oc = electrodeC->openCircuitVoltageSSRxn(0, 0);
    oc = electrodeC->openCircuitVoltage(0);
    printf("oc[0] = %g\n", oc);

    int nT = 50;
    deltaT = 2.0E-2;

    electrodeC->printElectrode();
    electrodeC->setPrintLevel(debug_print_lvl);
    //electrodeC->setPrintLevel(1);
    electrodeC->setDeltaTSubcycle(0.01);
    //electrodeC->detailedResidPrintFlag_ = 10;
    //electrodeC->enableExtraPrinting_ = 10;
    electrodeC->printCSVLvl_ = 3;

    int iOwn = 0;
    double atol = 1.0E-13;
    GFCEO_Electrode gee(electrodeC, atol, iOwn);
    DAE_Solver* integDAE = Zuzax::newDAE_Solver("IDA", gee);
    gee.set_DAE_Integrator(integDAE);
    gee.setTolerances(1.0E-6, 1.0E-13);
  
    for (int itimes = 0; itimes < nT; itimes++) {
      Tinitial = Tfinal;

      electrodeC->resetStartingCondition(Tinitial);
      
      Tfinal = Tinitial + deltaT;

      gee.integrate(deltaT);
      electrodeC->getMoleNumSpecies(molNum);
      double net[12];
      double amps = electrodeC->getIntegratedProductionRatesCurrent(net);
 
      cout << setw(15) << Tfinal << setw(15) << amps << endl;
      electrodeC->printElectrode();
      electrodeC->writeSolutionTimeIncrement();
  
    }

    delete cfC;
    delete electrodeC_input;
    delete electrodeC;
    appdelete();

    return retn;

  } catch (ZuzaxError) {
    showErrors();
    return -1;
  } catch (ZuzaxError) {
    showErrors();
    return -1;
  }
} 
