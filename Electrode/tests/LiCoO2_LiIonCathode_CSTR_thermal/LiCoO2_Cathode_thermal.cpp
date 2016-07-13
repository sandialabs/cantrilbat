/*
 * $Id: LiCoO2_Cathode_3_cc.cpp 496 2013-01-07 21:15:37Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "cantera/numerics/NonlinearSolver.h"
#include "cantera/equilibrium.h"
#include "Electrode_Factory.h"
#include "importPL.h"
#include "BE_BlockEntry.h"
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
  
    //ZZCantera::Electrode_CSTR_LiCoO2Cathode *electrodeC  = new ZZCantera::Electrode_CSTR_LiCoO2Cathode();
    ZZCantera::Electrode *electrodeC  = newElectrodeObject("CSTR_LiCoO2Cathode");

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

    electrodeC->doThermalPropertyCalculations_ = true;
  
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

    electrodeC->doThermalPropertyCalculations_ = true;

    double deltaT = 0.1;
    double Tinitial = 0.0;
    double Tfinal = 0.0;

    double molNum[10];

    electrodeC->setPhaseExistenceForReactingSurfaces(true);
 
    
    electrodeC->setVoltages(3.78, 0.0);

    double oc = electrodeC->openCircuitVoltageSSRxn(0);
    oc = electrodeC->openCircuitVoltage(0);
    printf("oc[0] = %g\n", oc);

    double ROP[10];
    int nT = 50;
    deltaT = 2.0E-4;
    electrodeC->printCSVLvl_ = 4;

    electrodeC->printElectrode();
    electrodeC->setPrintLevel(2);
    electrodeC->setPrintLevel(10);
    electrodeC->setDeltaTSubcycle(0.01);
    //electrodeC->detailedResidPrintFlag_ = 10;
    electrodeC->enableExtraPrinting_ = 10;
    double sa[10];

    //electrodeC->DO_NEW_METHOD_ = 1;
     nT = 2;
    for (int itimes = 0; itimes < nT; itimes++) {
      Tinitial = Tfinal;

      electrodeC->resetStartingCondition(Tinitial);
      
      Tfinal = Tinitial + deltaT;

      //electrodeC->integrateConstantCurrent(amps, deltaT, 4.0, 2.5);
      electrodeC->integrate(deltaT);

      electrodeC->getMoleNumSpecies(molNum);
      doublereal net[12];
      double amps = electrodeC->getIntegratedProductionRatesCurrent(net);
 
      cout << setw(15) << Tfinal << setw(15) << amps << endl;
      electrodeC->printElectrode();
      electrodeC->writeSolutionTimeIncrement();

 
      electrodeC->getSurfaceAreas(sa);
 
      ReactingSurDomain* rsd = electrodeC->reactingSurface(0);
      double cd0 = rsd->getCurrentDensityRxn();
      double current0 = cd0 * sa[0];

      printf("currentDensity0 = %g\n", cd0);
      printf("current0 = %g amps\n", current0);

      double cd1 = electrodeC->getNetSurfaceProductionRatesCurrent(0, ROP);
      printf("cd1 = %g\n", cd1);
      double current1 = cd1 * sa[0];
      printf("current1 = %g amps\n", current1);

      double nStoich, OCV, io, nu, beta, resist, cd2;
      double TT = electrodeC->temperature();
      
#ifdef DONOTREMOVE
      cd2 = rsd->getExchangeCurrentDensityFormulation(0, &nStoich, &OCV, &io, &nu, &beta, &resist);
#else
      bool okk = rsd->getExchangeCurrentDensityFormulation(0, nStoich, OCV, io, nu, beta, resist);
      if (okk) {
          cd2 = rsd->calcCurrentDensity(nu, nStoich, io, beta, TT, resist);
      } else {
          cd2 = 0.0;
      }
#endif

      printf("cd2 = %g\n", cd2);

      double current2 = cd2 * sa[0];
      printf("current2 = %g amps\n", current2);

      double q_over   = electrodeC->thermalEnergySourceTerm_overpotential(0);
      double q_entrop = electrodeC->thermalEnergySourceTerm_reversibleEntropy(0);
      double q_enth   = electrodeC->thermalEnergySourceTerm_EnthalpyFormulation(0);

      double q_enth2  = electrodeC->thermalEnergySourceTerm_EnthalpyFormulation_SingleStep();
      q_enth2 /= deltaT;

      double q_over2  =  electrodeC->thermalEnergySourceTerm_Overpotential_SingleStep();
      q_over2 /= deltaT;

      double q_entrop2 = electrodeC->thermalEnergySourceTerm_ReversibleEntropy_SingleStep();
      q_entrop2 /= deltaT;

      printf("    Check of Thermal Source terms: These should add up to each other:\n");
      printf(" q_over    = %g\n", q_over);
      printf(" q_over2   = %g\n", q_over2);
      printf(" q_entrop  = %g\n", q_entrop);
      printf(" q_entrop2 = %g\n", q_entrop2);
      printf(" q_enth    = %g\n", q_enth);
      printf(" q_enth2   = %g\n", q_enth2);
      double diff = q_enth - q_over - q_entrop;
      if (fabs(diff) < 1.0E-16) {
          diff = 0.0;
      }
      double diff2 = q_enth2 - q_over - q_entrop;
      if (fabs(diff2) < 1.0E-11) {
          diff2 = 0.0;
      }
      printf(" diff  =  q_enth -  q_over  - q_entrop = %g\n", diff);
      printf(" diff2 =  q_enth2 - q_over  - q_entrop = %g\n", diff2);
    }

    delete cfC;
    delete electrodeC_input;
    delete electrodeC;
    ZZCantera::appdelete();

    return retn;

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
} 
