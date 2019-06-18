/*
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "zuzax/equilibrium.h"
#include "zuzax/thermo/MolalityVPSSTP.h"

#include "zuzax/equil/vcs_MultiPhaseEquil.h"
//#include "zuzax/equil/vcs_prob.h"
//#include "zuzax/equil/vcs_solve.h"
#include "zuzax/equil/vcs_VolPhase.h"
#include "zuzax/thermo/IonsFromNeutralVPSSTP.h"
#include "zuzax/numerics/ResidEval.h"
#include "zuzax/numerics/NonlinearSolver_JAC.h"

#include "Electrode_input.h"
#include "Electrode_SimpleDiff.h"
#include "Electrode_RadialDiffRegions.h"

#include "Electrode_FD_Jacobian.h"

#include <sstream>
#include <iomanip>

using namespace std;
using namespace Zuzax;

// a lvl of one prints out the .csv file
int mpequil_debug_print_lvl = 1;
int VCS_Debug_Print_Lvl = 3;
double tfinalEnd = 20.;

void printUsage() {
    cout << "usage: LiCoO2_Cathode_Radial_cV [-h] [-help_cmdfile] [-d #] [anode.inp]"
         <<  endl;
    cout << "    -h               : Prints this help" << endl;
    cout << "    -help_cmdfile    : Prints a list of block commands understood by this parser "
      "- add cathode.inp for more information" << endl;
    cout << "   -d #              : Level of debug printing" << endl;
    cout << "   cathode.inp         : Command file (if missing, assume cathode.inp)" << endl;
    cout << endl;
}
//===================================================================================================================================
double getVoltage(double t)
{
    double Vstart = 3.8;
    double Vend   = 2.5;
    double v = Vstart + (Vend - Vstart) * t / (tfinalEnd);
    return v;
}
//===================================================================================================================================

int main(int argc, char **argv)
{
    int retn = 0;
    //bool doCathode = false;
    std::string commandFileNet = "cathode.inp";
    int numSubIntegrations;
    bool printInputFormat = false; // print cmdfile.txt format
    // printed usage

    //VCSnonideal::vcs_timing_print_lvl = 0;
    NonlinearSolver_JAC::s_TurnOffTiming = true;
    JacobianManager::s_print_Jac = true;

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
	    } else if (commandFileNet == "cathode.inp") {
		commandFileNet = tok;
	    } else {
		printUsage();
		exit(1);
	    }
	}
    }

    try {
	double volts;

	/*
	 *  Use this opportunity to set environment to get consistent printing for test cases
	 *  Note, this has moved to electrode_model_create routine, generally.
	 */ 
        Electrode::s_printLvl_DEBUG_SPECIAL = 1;
	Electrode::readEnvironmentalVariables();
  
	Zuzax::Electrode_SimpleDiff *electrodeC  = new Zuzax::Electrode_SimpleDiff();

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
	double molNum[10];
	double vvolts = electrodeC->openCircuitVoltageRxn(0, 0);
	printf(" OC volts = %g\n", vvolts);

	electrodeC->setPhaseExistenceForReactingSurfaces(true);
	volts = getVoltage(0.0);
	electrodeC->setVoltages(3.7, 0.0);
	double oc = electrodeC->openCircuitVoltageSSRxn(0);
	oc = electrodeC->openCircuitVoltage(0);
	printf("oc[0] = %g\n", oc);

	electrodeC->printCSVLvl_ = 3;

	double pmv[10];

	ThermoPhase *th = electrodeC->phasePtr("LiCoO2_Interstitials_cathode");
	th->getPartialMolarVolumes(pmv);
    
	electrodeC->setDeltaTSubcycle(2.5E-4);
	electrodeC->printElectrode();

	remove("soln.xml");

	//
	// We seek 20 amps / m2 for a total cross sectional area of 1 cm2
	//
	double amps;

        double rdel = 1.0E-4;
	Electrode_FD_Jacobian jac(electrodeC, rdel);
        bool baseAlreadyCalc = true;
        bool useDefaultDeltas = true;
        std::vector<double> centerpoint;
        std::vector<double> dof_Deltas(100, 0.0);
	jac.default_setup(centerpoint);
	electrodeC->setPrintLevel(1);
	//electrodeC->enableExtraPrinting_ = 10;
	//electrodeC->detailedResidPrintFlag_ = 10;
	//electrodeC->setMaxNumberSubCycles(40);

	double coul = 0.0;
 
	double net[12];
	int nT = 50;
	double deltaT = tfinalEnd / nT;

	FILE * fp = fopen("outputTable.csv","w");
	fprintf(fp, "Constant Current Curves \n");
	fprintf(fp, "deltaT = %g  \n", deltaT);
	fprintf(fp, "\n");
	fprintf(fp, "  Tfinal      ,       Coul     ,        Ah    ,        amps     ,  Volts , numSubs\n");

	double Tfinal = 0.0;
	for (int itimes = 0; itimes < nT; itimes++) {
	    Tinitial = Tfinal;
	    electrodeC->resetStartingCondition(Tinitial);
	    Tfinal = Tinitial + deltaT;
	    volts = getVoltage(Tfinal);
	    electrodeC->setVoltages(volts, 0.0);
	    volts = electrodeC->voltage();
	    //
	    // Integrate
	    numSubIntegrations = electrodeC->integrate(deltaT);

	    double Tfinalc = electrodeC->timeFinalFinal();
            if (Tfinalc != Tfinal) {
                printf("error in final time: %g %g\n", Tfinalc, Tfinal);
                exit(-1);
            }

	    amps = electrodeC->getIntegratedProductionRatesCurrent(net);
	    coul  += amps * deltaT;

	    fprintf(fp, " %12.6E ,  %12.6E , %12.6E , %12.6E , %12.6E ,  %12d\n", 
		    Tfinal, coul, coul/3600. , amps, volts, numSubIntegrations);

	    electrodeC->getMoleNumSpecies(molNum);
	    amps = electrodeC->getIntegratedProductionRatesCurrent(net);

	    cout << setw(15) << Tfinal << setw(15) << amps << numSubIntegrations << endl;

	    electrodeC->printElectrode();
	    printf(" NumSubs = %d\n", numSubIntegrations);
	    electrodeC->writeSolutionTimeIncrement();

	    //
	    // Calculate the jacobian
	    //  Fill in centerpoint vector, which are external DOFS. Then, calculate the perturbation delta for the 
	    //  numerical delta
	    //
	    jac.default_dofs_fill(centerpoint);
            jac.calc_Perturbations(centerpoint, dof_Deltas);
	    //
	    // Calculate the jacobian
	    //    Do a one-sided jacobian calculation
	    //
            jac.compute_oneSided_jacobian(centerpoint, deltaT, &dof_Deltas[0], useDefaultDeltas, baseAlreadyCalc);

	    //
	    // Print out the results
	    //
            jac.print_jacobian();

	}

	delete cfC;
	delete electrodeC_input;
	delete electrodeC;
	Zuzax::appdelete();

	return 0;

    } catch (ZuzaxError) {

	showErrors();
	return -1;
    }
} 
