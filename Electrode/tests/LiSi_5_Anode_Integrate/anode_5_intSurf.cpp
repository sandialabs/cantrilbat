/*
 * $Id: anode_5_intSurf.cpp 496 2013-01-07 21:15:37Z hkmoffa $
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

#include "Electrode_input.h"
#include "Electrode.h"
#include "Electrode_SimplePhaseChangeDiffusion.h"
#include "ExtraGlobalRxn.h"
#include "RxnMolChange.h"
#include "BlockEntry.h"

#include <stdio.h>

using namespace std;
using namespace Cantera;
using namespace VCSnonideal;
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


class ECurr : public Cantera::ResidEval
{
public:

  ECurr(Electrode *ee, double deltaT) :
    m_ee(ee),
    m_deltaT(deltaT)
  {

  }

  int nEquations() const {
    return 1;
  }
 int evalSS(const doublereal t, const doublereal * const x,
                       doublereal * const r) {
   m_ee->setVoltages(0.0, x[0]);
   m_ee->integrate(m_deltaT);
   double amps = m_ee->getIntegratedProductionRatesCurrent(srcNet);
   r[0] = amps;
   return 0;
  }


  Electrode *m_ee;
  double m_deltaT;
  double srcNet[50];
};


//======================================================================================================================


int main(int argc, char **argv)
{
  int i;
  int ip1 = 0;
 


  int retn = 0;
  //bool doCathode = false;
  string commandFileNet = "cell.inp";
  string commandFileA = "electrodeAnode.inp";
  string commandFileC = "electrodeCathode.inp";
  bool printInputFormat = false; // print cmdfile.txt format
  bool printedUsage = false; // bool indicated that we have already
  // printed usage

  VCSnonideal::vcs_timing_print_lvl = 0;

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
	    printedUsage = true;
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
		  if (lvl == 0) ip1 = 0;
		  else          ip1 = lvl; 
		  mpequil_debug_print_lvl = lvl;
		}
	      }  
	    }
	  } else {
	    printUsage();
	    printedUsage = true;
	    exit(1);
	  }
	}
      } else if (commandFileNet == "") {
	commandFileNet = tok;
      } else {
	printUsage();
	printedUsage = true;
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
    //IonsFromNeutralVPSSTP * ionicLiquid = new Cantera::IonsFromNeutralVPSSTP("LiKCl_recipmoltenSalt.xml");
    //  trans_ = Cantera::newTransportMgr("Liquid", ionicLiquid_, 1);

    Cantera::Electrode_SimplePhaseChangeDiffusion *electrodeA  = new Cantera::Electrode_SimplePhaseChangeDiffusion();
    
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

    double deltaT = 2.0;

    double Tinitial = 0.0;
    double Tfinal = 0.0;
    int inLi13Si4 = electrodeA->globalSpeciesIndex("Li13Si4(S)");
    int inLi7Si3 = electrodeA->globalSpeciesIndex("Li7Si3(S)");
    int inLi_i_ = electrodeA->globalSpeciesIndex("Li(i)");
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
    ThermoPhase *tpLyte = &(electrodeA->thermo(1));
    int iLip = tpLyte->speciesIndex("Li+");
    int iKp = tpLyte->speciesIndex("K+");
    int iClm = tpLyte->speciesIndex("Cl-");
    x[iLip] = 0.3;
    x[iKp] = 0.2;
    x[iClm] = 0.5;
    electrodeA->setPhaseMoleNumbers(1, x);

    electrodeA->getMoleNumSpecies(DATA_PTR(molNum));
    electrodeA->getMoleFractions(DATA_PTR(molFrac));
    electrodeA->getSurfaceAreas(sa);
    sa[1] = sa[0];
    printf("moleNum: Li13Si4(S) = %g,  Li7Si3() = %g\n",
	   molNum[inLi13Si4],  molNum[inLi7Si3]); 
    printf("molNum: Li(i) = %g\n", molNum[inLi_i_]);
    printf("molFrac: Li13Si4(S) = %g,  Li7Si3() = %g\n",
	   molFrac[inLi13Si4],  molFrac[inLi7Si3]); 
    printf("molFrac: Li(i) = %g\n", molFrac[inLi_i_]);

    //int inph_Interst = electrodeA->globalPhaseIndex("Li7Si3_Interstitial");
    //ThermoPhase &tpInt = electrodeA->thermo(inph_Interst);
    double emu_Lip = -1.0;
    double emu_Li13S = -1.0;
    double emu_Li7S = -1.0;
    double emu_LiInt = -1.0;
    double emu_electron = -1.0;

    for (int ip = 0; ip < electrodeA->nPhases(); ip++) {
      ThermoPhase *tp_ptr = & electrodeA->thermo(ip);
      string pname = tp_ptr->name();
      printf(" phase %d = %s\n", ip, pname.c_str());
      tp_ptr->getChemPotentials(DATA_PTR(mu));
      tp_ptr->getElectrochemPotentials(DATA_PTR(emu));
      for (int k = 0; k < tp_ptr->nSpecies(); k++) {
	string sname = tp_ptr->speciesName(k);
	int gsi = electrodeA->globalSpeciesIndex(sname);
	printf("   %-16s %4d %-12.3E %-12.3E\n", sname.c_str(), gsi,  mu[k], emu[k]);
	if (sname == "Li+") {
	  emu_Lip = emu[k];
	}
	if (sname == "Li13Si4(S)") {
	  emu_Li13S = emu[k];
	}
	if (sname == "Li7Si3(S)") {
	  emu_Li7S = emu[k];
	}
	if (sname == "Li(i)") {
	  emu_LiInt = emu[k];
	}
	if (sname == "electron_Li_LiCl") {
	  emu_electron = emu[k];
	}
      }
    }
    
    //   electrodeA->followElectrolyteMoles = 0;
    double deltag_ext = emu_electron +  emu_Lip -  emu_LiInt;
    printf("deltag_ext (Li(i) -> E + Li+) = %12.3E\n", deltag_ext); 

    double deltag_int = emu_LiInt + 4./11.*emu_Li7S - 3./11.*emu_Li13S;
    printf("deltag_int (3./11. Li13Si4(S) -> 4/11 Li7Si3(S) + Li(i)) = %12.3E\n", deltag_int); 

    electrodeA->printElectrode();
    for (int itimes = 0; itimes < 30; itimes++) {

      Tinitial = Tfinal;

      electrodeA->resetStartingCondition(Tinitial);
      
      Tfinal = Tinitial + deltaT;

      electrodeA->integrate(deltaT);
      electrodeA->getMoleNumSpecies(DATA_PTR(molNum));
      printf("Li13Si4(S) = %g,  Li7Si3() = %g\n", molNum[inLi13Si4],  molNum[inLi7Si3]);
      electrodeA->printElectrode();

      // accept step
     
    }
    Tinitial = Tfinal; 
    electrodeA->resetStartingCondition(Tinitial);
 
    delete cfA;



    printf("\n Results from mpequil:\n");
    printf("\n");
    printf("Temperature = %g Kelvin\n",  electrodeA->temperature());
    printf("Pressure    = %g ", electrodeA->pressure());
    printf("Pa\n");
 
    printf("\n");
    printf("----------------------------------------"
	   "---------------------\n");
    printf(" Name             Mole_Number     Mole_Fraction     Chem_Potential");
    printf(" (J/kmol)\n");
    printf("--------------------------------------------------"
	   "-----------\n");
    string spname;
   
    for (i = 0; i < electrodeA->nSpecies(); i++) {
      spname = electrodeA->speciesName(i);
      printf("%-12s", spname.c_str());
      printf("  %15.6e %15.6e  %15.6e\n", electrodeA->moleNumSpecies(i),
	     electrodeA->moleFraction(i), electrodeA->speciesElectrochemPotential(i));
    }

 
    printf("------------------------------------------"
	   "-------------------\n");   
   
    for (int iph = 0; iph < (int) electrodeA->nPhases(); iph++) {
      ThermoPhase *tp = &(electrodeA->thermo(iph));
      string phname = tp->id();
      printf("%-12s", phname.c_str());
      printf("  %15.6e %15.6e\n", electrodeA->phaseMoles(iph), electrodeA->phaseVoltage(iph));
    }


    delete electrodeA_input;
    delete electrodeA;


    Cantera::appdelete();

    return retn;

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
} 
