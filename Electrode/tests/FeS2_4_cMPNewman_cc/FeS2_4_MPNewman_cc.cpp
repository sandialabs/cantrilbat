/*
 * $Id: FeS2_4_MPNewman_cc.cpp 496 2013-01-07 21:15:37Z hkmoffa $
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
#include "Electrode_MP_RxnExtent_FeS2.h"

#include <sstream>
#include <iomanip>

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

using namespace VCSnonideal;
using namespace std;

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
  
    ZZCantera::Electrode_MP_RxnExtent_FeS2 *electrodeC  = new ZZCantera::Electrode_MP_RxnExtent_FeS2();

    ELECTRODE_KEY_INPUT *electrodeC_input = new ELECTRODE_KEY_INPUT();
    
    std::string commandFileC = "cathode.inp";


	
    // Initialize a block input structure for the command file

    BEInput::BlockEntry *cfC = new BEInput::BlockEntry("command_file");

 
    // Go get the problem description from the input file
    electrodeC_input->printLvl_ = 3;
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
   
    double deltaTgoal = 0.1;

    double Tinitial = 0.0;
    double Tfinal = 0.0;

    electrodeC->setVoltages(1.95, 0.0);

 
    double molNum[10];



    electrodeC->setPhaseExistenceForReactingSurfaces(true);
 
    
    electrodeC->setVoltages(1.60, 0.0);

    double oc = electrodeC->openCircuitVoltageSSRxn(1);
    oc = electrodeC->openCircuitVoltage(1);
    printf("oc[0] = %g\n", oc);
    FILE * fp = fopen("outputTable.csv","w");


    int nT = 500;
    deltaTgoal = 2.0E-1;

    electrodeC->printElectrode();
    electrodeC->setPrintLevel(1);
    electrodeC->setDeltaTSubcycle(0.01);
    //electrodeC->detailedResidPrintFlag_ = 2;
    electrodeC->detailedResidPrintFlag_ = 0;
    electrodeC->enableExtraPrinting_ = 0;
    electrodeC->printCSVLvl_ = 1;
    //electrodeC->printCSVLvl_ = 9;

    deltaTgoal = 1.0;
    double coul = 0.0;

    double amps = -50.0;
    fprintf(fp, "Constant Current Curves \n");
    fprintf(fp, "amps = %g  \n", amps);
    fprintf(fp, "deltaT = %g  \n", deltaTgoal);
    fprintf(fp, "\n");
    fprintf(fp, "  Tfinal      ,      Coul     ,     Ah   , Volts \n");

    for (int itimes = 0; itimes < nT; itimes++) {
      Tinitial = Tfinal;

      electrodeC->resetStartingCondition(Tinitial);
     
      electrodeC->setPrintLevel(0);
      electrodeC->printCSVLvl_ = 0;
      double deltaT = deltaTgoal; 
      double volts =  electrodeC->integrateConstantCurrent(amps, deltaT, 2.2, 1.3);
      Tfinal = Tinitial + deltaT;

      electrodeC->printCSVLvl_ = 1;
      electrodeC->setPrintLevel(1);
      electrodeC->integrate(deltaT);
      electrodeC->setPrintLevel(0);
      electrodeC->printCSVLvl_ = 0;
    
      electrodeC->getMoleNumSpecies(molNum);
      doublereal net[12];
      double amps = electrodeC->getIntegratedProductionRatesCurrent(net);
      coul  += amps * deltaT;
      fprintf(fp, " %12.6E ,  %12.6E , %12.6E , %12.6E\n", Tfinal, coul, coul/3600. , volts);

    
      electrodeC->printElectrode();
      electrodeC->writeSolutionTimeIncrement();


      SubIntegrationHistory sih1 = electrodeC->timeHistory();
      //sih1.print(3);

      // test new capability
      {
	  electrodeC->setTimeHistoryBaseFromCurrent();

	  Electrode_Exterior_Field_Interpolation_Scheme_Enum ef =  T_FINAL_CONST_FIS;
	  Subgrid_Integration_RunType_Enum subIntegrationType =  FVDELTA_TIMEINTEGRATION_SIR;
	  electrodeC->integrate(deltaT, 1.0E-3, ef, subIntegrationType);

	  electrodeC->getMoleNumSpecies(molNum);
	  doublereal net[12];
	  double amps = electrodeC->getIntegratedProductionRatesCurrent(net);
	  coul  += amps * deltaT;
	  fprintf(fp, " %12.6E ,  %12.6E , %12.6E , %12.6E\n", Tfinal, coul, coul/3600. , volts);

    
	  //electrodeC->printElectrode();

	  SubIntegrationHistory&sih2 = electrodeC->timeHistory();
	  //sih2.print(3);
	  if (sih1 != sih2) {
	      printf("we have a prob\n");
	  }
      } 

    }

    fclose(fp);

      delete cfC;
      delete electrodeC_input;
      delete electrodeC;
      appdelete();

    return 0;

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
} 
