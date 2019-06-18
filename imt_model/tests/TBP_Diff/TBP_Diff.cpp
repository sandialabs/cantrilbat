/*
 * $Id: AtoB_1.cpp 222 2012-06-26 21:55:32Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */



#include "zuzax/equil/vcs_MultiPhaseEquil.h"
#include "zuzax/equil/vcs_solve.h"
#include "zuzax/equil/vcs_VolPhase.h"
#include "zuzax/equil/vcs_internal.h"
#include "zuzax/numerics/ResidEval.h"
#include "zuzax/numerics/NonlinearSolver_JAC.h"

#include "InterfacialMassTransfer_input.h"
#include "InterfacialMassTransfer.h"
#include "InterfacialMassTransfer_1to1Distrib.h"
#include "imtPSS_NoSurf_DiffBL.h"

#include <sstream>
#include <iomanip>

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


//=====================================================================================================


int main(int argc, char **argv)
{
  int retn = 0;
  //bool doCathode = false;
  string commandFileNet="";
  // printed usage

  NonlinearSolver_JAC::s_TurnOffTiming = true;

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

  
    //  retn = cell_input(commandFileNet);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }
  
  
    Zuzax::imtPSS_NoSurf_DiffBL* iface  = new Zuzax::imtPSS_NoSurf_DiffBL();

    
    Zuzax::IMT_KEY_INPUT *face_input = new IMT_KEY_INPUT();
    
    std::string commandFile = "interface.inp";
  
    /**	
     * Initialize a block input structure for the command file
     */
    BEInput::BlockEntry *cf = new BEInput::BlockEntry("command_file");

    /*	
     * Go get the problem description from the input file
     */
    retn = imt_input(face_input, commandFile, cf);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

 
    retn = iface->model_create(face_input);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }


    retn = iface->setInitialConditions(face_input);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }


    std::vector<double> molNum(30);
    int nT = 5;
    double deltaT = 1.0E-5;
    double Tinitial = 0.0;
    double Tfinal = 0.0;
    iface->relativeLocalToGlobalTimeStepMinimum_ = 0.05;
    for (int itimes = 0; itimes < nT; itimes++) {
      Tinitial = Tfinal;
      iface->resetStartingCondition(Tinitial);
      Tfinal = Tinitial + deltaT;
      iface->integrate(deltaT, 1.0E-3);
      iface->getMoleNumSpecies(DATA_PTR(molNum));
      double net[30];
      iface->getIntegratedProductionRates(net);
      cout << setw(15) << Tfinal << setw(15) << 0.0 << endl;
      iface->printInterfacialMassTransfer(1, false);
    }

    delete cf;
    delete face_input;
    delete iface;
    Zuzax::appdelete();

    return retn;

  } catch (ZuzaxError) {

    showErrors();
    return -1;
  }
} 
