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
#include "zuzax/thermo/IonsFromNeutralVPSSTP.h"
#include "zuzax/numerics/ResidEval.h"
#include "zuzax/numerics/RootFind.h"
#include "zuzax/numerics/NonlinearSolver_JAC.h"

#include "Electrode_input.h"
#include "Electrode_MP_RxnExtent.h"

#include <sstream>
#include <iomanip>
#include <cstdio>

using namespace std;
using namespace Zuzax;

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



int main(int argc, char **argv)
{
  int retn = 0;
  //bool doCathode = false;
  string commandFileNet = "cell.inp";
  string commandFileA = "electrodeAnode.inp";
  //  string commandFileC = "electrodeCathode.inp";
  string commandFileC = "cathode.inp";
  // printed usage

  vcs_nonideal::vcs_timing_print_lvl = 0;
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
  
    Electrode_MP_RxnExtent *electrodeC  = new Electrode_MP_RxnExtent();

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
    retn = 0;

    std::vector<double> valAmps; 
    valAmps.push_back(50.0);

    std::vector<double> rval;
    rval.push_back(1.0);


    std::vector<double> val;
    val.push_back(1.0E-7);
    val.push_back(1.0E-8);
    val.push_back(1.0E-9);
    val.push_back(1.0E-10);

   for (int aTimes = 0; aTimes < (int) valAmps.size(); aTimes++) {

    for (int iTimes = 0; iTimes < (int) val.size() ; iTimes++) {

    for (int rTimes = 0; rTimes < (int) rval.size(); rTimes++) {

    bool startNewRecord = true;
    electrodeC->diffusionCoeffRegions_[0] = val[iTimes];
    electrodeC->diffusionCoeffRegions_[1] = val[iTimes];
    electrodeC->diffusionCoeffRegions_[2] = val[iTimes];
    electrodeC->diffusionCoeffRegions_[3] = val[iTimes];
    electrodeC->diffusionCoeffRegions_[4] = val[iTimes];

    electrodeC->rxnPerturbRegions_[0] = rval[rTimes];
    electrodeC->rxnPerturbRegions_[1] = rval[rTimes];
    electrodeC->rxnPerturbRegions_[2] = rval[rTimes];
    electrodeC->rxnPerturbRegions_[3] = rval[rTimes];

    double deltaT = 0.1;
    double Tinitial = 0.0;
    double Tfinal = 0.0;

    electrodeC->setVoltages(1.95, 0.0);
    double molNum[10];

    electrodeC->setPhaseExistenceForReactingSurfaces(true);
    electrodeC->setVoltages(1.90, 0.0);
    electrodeC->setState_relativeExtentRxn(0.0);
    electrodeC->setTime(0.0);

    double oc = electrodeC->openCircuitVoltageSSRxn(1);
    oc = electrodeC->openCircuitVoltage(1);
    printf("oc[0] = %g\n", oc);

    double coul = 0.0;

    double ampsIn = -valAmps[aTimes];  
    int nT = 400;
    deltaT = 5.0E-1;
    deltaT = 0.5 * 50.0 / (-ampsIn);
 
    std::string sname = "outputTable_" +  int2str(aTimes) + "_" + int2str(iTimes) +  "_" + int2str(rTimes) + ".csv";
    FILE * fp = fopen(sname.c_str(),"w");

    electrodeC->updateState();
    electrodeC->extractInfo();
    electrodeC->printElectrode();
    electrodeC->setPrintLevel(0);
    electrodeC->setDeltaTSubcycle(0.01);
    electrodeC->detailedResidPrintFlag_ = 0;
    electrodeC->printCSVLvl_ = 2;
    electrodeC->enableExtraPrinting_ = 0;
    fprintf(fp, "Constant Current Curves \n");
    fprintf(fp, "amps = %g  \n", ampsIn);
    fprintf(fp, "DiffC = %g \n", val[iTimes]);
    fprintf(fp, "deltaT = %g  \n", deltaT);
    fprintf(fp, "Rperturb = %g \n", rval[rTimes]);
    fprintf(fp, "\n");
    fprintf(fp, "  Tfinal      ,      Coul     ,     Ah   , Volts , ElectronsPerMoleFeS2 \n");


    ampsIn =  -valAmps[aTimes];
    double deltatsmall = 1.0E-5 * deltaT;
    double relE = electrodeC->relativeExtentRxn(Tinitial);
    double volts =  electrodeC->integrateConstantCurrent(ampsIn, deltatsmall, 2.2, 0.3);
    fprintf(fp, " %13.6E ,  %13.6E , %13.6E , %13.6E, %13.6E\n", Tfinal, -coul, -coul/3600. , volts, relE);
    
    double deltaTgoal = deltaT;
    double TFINAL = nT * deltaTgoal;
    Tfinal = 0.0;
    int itimes = 0;
    while (Tfinal < TFINAL) {
      bool iprint = false;
      if (itimes < 20) {
        iprint = true;
      }
      if (itimes % 10 == 0) {
        iprint = true;
      }
      Tinitial = Tfinal;
      deltaT = deltaTgoal;
      ampsIn =  -valAmps[aTimes];
      volts =  electrodeC->integrateConstantCurrent(ampsIn, deltaT, 2.2, 0.0);
      Tfinal = Tinitial + deltaT;
      if (Tfinal >= 0.97 * TFINAL) {
        iprint = true;
      }

      electrodeC->getMoleNumSpecies(molNum);
      double net[12];
      double amps = electrodeC->getIntegratedProductionRatesCurrent(net);
 
      coul  += amps * deltaT;
      relE = electrodeC->relativeExtentRxn(Tfinal);
      fprintf(fp, " %13.6E ,  %13.6E , %13.6E , %13.6E, %13.6E\n", Tfinal, -coul, -coul/3600. , volts, relE);

      cout << setw(15) << Tfinal << setw(15) << amps << endl;
      if (iprint) {
        electrodeC->printElectrode();
      } else {
        printf("main: Skipping Tinitial %g to Tfinal %g printout\n", Tfinal, Tinitial); 
      }
      electrodeC->writeSolutionTimeIncrement(startNewRecord);
      startNewRecord = false;
      electrodeC->resetStartingCondition(Tfinal);
       
      itimes++; 
      if (itimes >= nT) {
        printf("extra its\n");
      }
    }

    fclose(fp);
    
    }
    }
    }

    delete cfC;
    delete electrodeC_input;
    delete electrodeC;
    appdelete();

    return retn;

  } catch (ZuzaxError) {

    showErrors();
    return -1;
  }
} 
