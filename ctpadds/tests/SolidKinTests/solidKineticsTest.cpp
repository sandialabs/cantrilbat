/*
 * $Id: LiCoO2_Cathode_3_cc.cpp 496 2013-01-07 21:15:37Z hkmoffa $
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
#include "cantera/thermo.h"
#include "cantera/thermo/StoichSubstance.h"
#include "cantera/thermo/MetalSHEelectrons.h"

#include "cantera/kinetics.h"

#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/thermo/IonsFromNeutralVPSSTP.h"
#include "cantera/numerics/RootFind.h"
#include "cantera/numerics/NonlinearSolver.h"

#include "cantera/thermo/IdealSolidSolnPhase.h"

#include "SolidKinetics.h"

#include "BE_BlockEntry.h"

#include <cstdio>
#include <sstream>
#include <iomanip>

using namespace std;
using namespace Cantera;
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
    Cantera::NonlinearSolver::s_TurnOffTiming = true;

    int retn = 0;
    try {
	//Cantera::Electrode *electrodeC  = Cantera::Electrode_Factory_TB::factory()->newElectrodeObject("CSTR_ZnAnode");

 
	doublereal netROP[10], fwdROP[10], revROP[10];
 
	std::string commandFileC = "anode.inp";
	
	// Initialize a block input structure for the command file
//	BEInput::BlockEntry *cfC = new BEInput::BlockEntry("command_file");



        XML_Node* xx = Cantera::get_XML_File("IdealSolidSolnPhaseExample.xml");
        IdealSolidSolnPhase *ss = new IdealSolidSolnPhase(*xx);


        XML_Node& xmlPhase = ss->xml();


	SolidKinetics* eK = new SolidKinetics(ss);

        vector<ThermoPhase*> tpList;
        tpList.push_back(ss);
    
	bool ok = importKinetics(xmlPhase, tpList, eK);
	if (!ok) {
	    exit(-1);
	}

	double TT = 298.15;
	double pres = OneAtm;
	ss->setState_TP(TT, pres);


	double molNum[10];
	std::fill(molNum, molNum+10, 0.0);

	const int C2H2_index       = eK->kineticsSpeciesIndex("C2H2-graph");
	const int C_index       = eK->kineticsSpeciesIndex("C-graph");
	const int H2_index       = eK->kineticsSpeciesIndex("H2-solute");



	//
	// Integrate at constant voltage as a function of increasing KOH concentration
	//
	// Jump to failed step
	//molNum[oh_index] =0.23;
	double ocv = 0.0;


	int numTimes = 0; 
	//while( molNum[oh_index] < 0.4 )
	{

	    molNum[C2H2_index] = 0.9;
	    molNum[C_index] = 0.01;
	    molNum[H2_index] = 1.0 - molNum[C2H2_index] -  molNum[C_index] ;

            ss->setState_TPX(TT, pres, molNum);

            double xmole[20];
	    ss->getMoleFractions(xmole);
	    printf("moleF C2H2 = %g\n", xmole[C2H2_index]);
	    printf("moleF C  = %g\n", xmole[C_index]);
	    printf("moleF H2  = %g\n", xmole[H2_index]);

	    printf("============================================================================================\n");
	    printf("============================================================================================\n");

	    std::cout << " c_H2 = " << molNum[H2_index];
	    std::cout << " OCV = " << ocv << endl;

	    double deltaG[20];
	    eK->getDeltaGibbs(deltaG);





            eK->getNetRatesOfProgress(netROP); 

	    eK->getFwdRatesOfProgress(fwdROP);
	    eK->getRevRatesOfProgress(revROP);


            printf("  %g     %g    %g   \n", netROP[0], fwdROP[0] , revROP[0]);

            double Kappa[10];
            eK->getEquilibriumConstants(Kappa);

            printf(" %g \n", Kappa[0]);

            
            double ca[10];
            ss->getActivityConcentrations(ca);
            printf(" %g   %g   %g \n", ca[0], ca[1], ca[2]);

            double aa[10];
            ss->getActivities(aa);
            printf(" %g   %g   %g \n", aa[0], aa[1], aa[2]);

            double cs[10];
            cs[0] = ss->standardConcentration(0);
            cs[1] = ss->standardConcentration(1);
            cs[2] = ss->standardConcentration(2);
            printf(" %g   %g   %g \n", cs[0], cs[1], cs[2]);


            double kf = fwdROP[0] / ca[0];

            double revcalc = kf / Kappa[0] * ca[1] * ca[1] * ca[2];

            printf(" %g \n", revcalc);
	}

   
	Cantera::appdelete();

	return retn;

    } catch (CanteraError) {

	showErrors();
	return -1;
    }
} 
