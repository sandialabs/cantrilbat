/*
 * $Id: LiSi_5_MultiPlat.cpp 274 2013-03-26 17:04:11Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "cantera/equilibrium.h"
#include "cantera/thermo/MolalityVPSSTP.h"

#include "cantera/thermo/IonsFromNeutralVPSSTP.h"
#include "cantera/numerics/ResidEval.h"
#include "cantera/numerics/RootFind.h"
#include "cantera/numerics/NonlinearSolver.h"

#include "Electrode_input.h"
#include "Electrode.h"
#include "Electrode_MP_RxnExtent.h"

#include <stdio.h>

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif
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


int main(int argc, char **argv)
{
  int retn = 0;
  //bool doCathode = false;
  string commandFileNet = "cell.inp";
  string commandFileA = "electrodeAnode.inp";
  string commandFileC = "electrodeCathode.inp";
  // printed usage

  ZZCantera::NonlinearSolver::s_TurnOffTiming = true;


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

    ZZCantera::Electrode_MP_RxnExtent *electrodeA  = new ZZCantera::Electrode_MP_RxnExtent();
    
    ELECTRODE_KEY_INPUT *electrodeA_input = new ELECTRODE_KEY_INPUT();
    
    std::string commandFileA = "anode.inp";
   
  
    /**	
     * Initialize a block input structure for the command file
     */
    BEInput::BlockEntry *cfA = new BEInput::BlockEntry("command_file");

    /*	
     * Go get the problem description from the input file
     */
    electrodeA_input->printLvl_ = 5;
    retn = electrodeA_input->electrode_input(commandFileA, cfA);

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

    double deltaT = 100;

    double Tinitial = 0.0;
    double Tfinal = 0.0;
    int inLi13Si4 = electrodeA->globalSpeciesIndex("Li13Si4(S)");
    int inLi7Si3 = electrodeA->globalSpeciesIndex("Li7Si3(S)");
    int inLi12Si7 = electrodeA->globalSpeciesIndex("Li12Si7(S)");
    int inSi = electrodeA->globalSpeciesIndex("Si(S)");
    if (inLi13Si4 < 0) {
      throw CanteraError("main", "species not found");
    }

    int nspGlobal = electrodeA->nSpecies();
    vector<double> molNum(nspGlobal, 0.0);
    vector<double> molFrac(nspGlobal, 0.0);
    electrodeA->setVoltages(0.0, -0.20);

    double x[30];
    double mu[30];
    double emu[30];
    for (int k = 0; k < 30; k++) {
      x[k] = 0.0;
      mu[k] = 0.0;
      emu[k] = 0.0;
    }
    double sa[10];

    // Set the Electrolyte mole fractions
    int ipLyte = electrodeA->globalPhaseIndex("LiKCl_electrolyte");
    ThermoPhase *tpLyte = &(electrodeA->thermo(ipLyte));
    int iLip = tpLyte->speciesIndex("Li+");
    int iKp = tpLyte->speciesIndex("K+");
    int iClm = tpLyte->speciesIndex("Cl-");
    x[iLip] = 0.3E-5;
    x[iKp] = 0.2E-5;
    x[iClm] = 0.5E-5;
    electrodeA->setPhaseMoleNumbers(ipLyte, x);

    electrodeA->getMoleNumSpecies(DATA_PTR(molNum));
    electrodeA->getMoleFractions(DATA_PTR(molFrac));
    electrodeA->getSurfaceAreas(sa);
    sa[1] = sa[0];
    printf("moleNum: Li13Si4(S) = %g,  Li7Si3() = %g\n",
	   molNum[inLi13Si4],  molNum[inLi7Si3]); 
    printf("molFrac: Li13Si4(S) = %g,  Li7Si3() = %g\n",
	   molFrac[inLi13Si4],  molFrac[inLi7Si3]); 
 


    for (size_t ip = 0; ip < electrodeA->nPhases(); ip++) {
      ThermoPhase *tp_ptr = & electrodeA->thermo(ip);
      string pname = tp_ptr->name();
      printf(" phase %d = %s\n", (int) ip, pname.c_str());
      tp_ptr->getChemPotentials(DATA_PTR(mu));
      tp_ptr->getElectrochemPotentials(DATA_PTR(emu));
      for (size_t k = 0; k < tp_ptr->nSpecies(); k++) {
	string sname = tp_ptr->speciesName(k);
	int gsi = electrodeA->globalSpeciesIndex(sname);
	printf("   %-16s %4d %-12.3E %-12.3E\n", sname.c_str(), gsi,  mu[k], emu[k]);

      }
    }
    
    //   electrodeA->followElectrolyteMoles = 0;
  

    FILE *fp;
    char fileName[80] = "voltHist.csv";
    fp = fopen(fileName, "w");

    electrodeA->setVoltages(0.10, 0.0);

    double oc[3];
    oc[0] = electrodeA->openCircuitVoltage(0);
    printf("oc[0] = %g\n", oc[0]);
    oc[1] = electrodeA->openCircuitVoltage(1);
    printf("oc[1] = %g\n", oc[1]);
    oc[2] = electrodeA->openCircuitVoltage(2);
    printf("oc[2] = %g\n", oc[2]);
    double coul = 0.0;
    double amps = 0.2455;
    deltaT = 100;

    fprintf(fp, "Constant Current Curves \n");
    fprintf(fp, "amps = %g  \n", amps);
    fprintf(fp, "deltaT = %g  \n", deltaT);
    fprintf(fp, "\n");
    fprintf(fp, "  Tfinal    , amps  ,      Coul     ,     Ah   , Volts   n_Li13Si4,  n_Li7Si3, n_Li12Si7, n_Si , OCV_0, OCV_1, OCV_2, sa_0, sa_1, sa_2\n");
    double volts;

    electrodeA->printElectrode();
    for (int itimes = 0; itimes < 130; itimes++) {

      Tinitial = Tfinal;

      electrodeA->resetStartingCondition(Tinitial);
      
      Tfinal = Tinitial + deltaT;

      volts =  electrodeA->integrateConstantCurrent(amps, deltaT);
      coul  += amps * deltaT;

      electrodeA->getMoleNumSpecies(x);
      electrodeA->getSurfaceAreas(sa);

      if (isnan(volts)) break;

      fprintf(fp, " %12.5E ,  %12.5E , %12.5E , %12.5E , %12.5E , %12.5E , %12.5E , %12.5E , %12.5E , %12.5E , %12.5E , %12.5E , %12.5E , %12.5E , %12.5E\n", Tfinal, amps, coul, coul/3600. , volts, x[inLi13Si4],  x[inLi7Si3], x[inLi12Si7], x[inSi] , oc[0], oc[1], oc[2] , sa[0], sa[1], sa[2] );
      
      //   electrodeA->getMoleNumSpecies(DATA_PTR(molNum));
      electrodeA->printElectrode();

      // accept step
     
    }
    int nSteps = 0;
    volts += 0.01;
    electrodeA->printCSVLvl_ = 1;
    electrodeA->setVoltages(volts, 0.0);
    for (int itimes = 0; itimes < 50; itimes++) {
      
      Tinitial = Tfinal;

      electrodeA->resetStartingCondition(Tinitial);
      
      Tfinal = Tinitial + deltaT;
    
      int iSteps = electrodeA->integrate(deltaT);
      nSteps += iSteps;
      
      amps = electrodeA->getIntegratedProductionRatesCurrent( x );
      coul  += amps * deltaT;

      electrodeA->getMoleNumSpecies( x ); 
      electrodeA->getSurfaceAreas(sa);

      if (isnan(volts)) break;
    
      fprintf(fp, " %12.5E ,  %12.5E , %12.5E , %12.5E , %12.5E , %12.5E , %12.5E , %12.5E , %12.5E ,"
                  " %12.5E , %12.5E , %12.5E , %12.5E , %12.5E , %12.5E\n",
                   Tfinal, amps, coul, coul/3600. , volts, x[inLi13Si4],  x[inLi7Si3], x[inLi12Si7], 
                   x[inSi] ,  oc[0], oc[1], oc[2], sa[0], sa[1], sa[2] );
      
      electrodeA->printElectrode();
      // accept step

    }


 
    fclose(fp);
    delete cfA;

    printf("      ===============================================================\n");
    printf("          RUN FINISHED Number of subcycle steps = %d\n", nSteps);
    printf("      ===============================================================\n");



    delete electrodeA_input;
    delete electrodeA;


    ZZCantera::appdelete();

    return retn;

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
} 
