/* ======================================================================= */
/* $RCSfile: Margules_test_2.cpp,v $ */
/* $Author: hkmoffa $ */
/* $Date: 2013-01-07 15:24:38 -0700 (Mon, 07 Jan 2013) $ */
/* $Revision: 501 $ */
/* ======================================================================= */



#include "cantera/thermo.h"
#include "cantera/thermo/MargulesVPSSTP.h"

#include "cantera/equilibrium.h"

#include "cantera/thermo.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/base/logger.h"
#include <cmath>


#include <cstdio>



using namespace std;
using namespace Cantera;

int CHECK_DEBUG_MODE = 0;

void printUsage() {
    cout << "usage: Margules_test_2 [-h] [-d #] " <<  endl;
    cout <<"                -> Everything is hardwired" << endl;
}


int main(int argc, char **argv)
{
  VCSnonideal::vcs_timing_print_lvl = 0;
  int printLvl = 2;
  int retn = 0;
  bool printedUsage = false;
  bool printInputFormat = false;
  
  try {

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
		    printLvl = lvl;
		  }
		}
	      }
	    } else {
	      printUsage();
	      printedUsage = true;
	      exit(1);
	    }
	  }
	} else {
	  printUsage();
	  printedUsage = true;
	  exit(1);
	}
      }
    }


    double x[20];
    for (int k = 0; k < 20; k++) {
      x[k] = 0.0;
    }
    double um[20];
    double pres = OneAtm;

    MargulesVPSSTP *salt = new MargulesVPSSTP(1);

    //int iKCl_l = salt->speciesIndex("KCl(L)");
    //int iLiCl_l = salt->speciesIndex("LiCl(L)");

    string f_licl = "LiCl_solid.xml";
    string id = "LiCl(S)";
    Cantera::ThermoPhase *LiCl_solid = Cantera::newPhase(f_licl, id);
    LiCl_solid->setState_TP(298.15, pres);

    string f_kcl = "KCl_solid.xml";
    id = "KCl(S)";
    Cantera::ThermoPhase *KCl_solid = Cantera::newPhase(f_kcl, id);
    KCl_solid->setState_TP(298.15, pres);

    /*
     * set states
     */
    x[0] = 0.7;
    x[1] = 1.0 - x[0];
    double T = 273.15 + 353.;
    salt->setState_TPX(T, OneAtm, x);
    LiCl_solid->setState_TP(T, OneAtm);
    KCl_solid->setState_TP(T, OneAtm);

    salt->getChemPotentials(um);
    double um_B = um[0];
    double um_A = um[1];
    printf("um[0] = %g\n", um[0]);
    printf("um[1] = %g\n", um[1]);

    LiCl_solid->getChemPotentials(um);
    double um_licls = um[0];
    KCl_solid->getChemPotentials(um);
    double um_kcls = um[0];


    double um_solid = x[0] * um_licls + x[1] * um_kcls;       

    double um_liquid = x[0] * um_B + x[1] * um_A;

    bool solidStable = true;
    if (um_liquid < um_solid) {
      solidStable = false;
    }

    printf(" %g      %g    %g   %g", T-273.15, um_solid, um_liquid, um_liquid-um_solid);
    if (solidStable) {
      printf(" Solid Stable\n");
    } else {
      printf(" Liquid Stable\n");
    } 
    printf("         LiCl: solid mu = %g liquid mu  %g\n", um_licls , um_B );
     

    Cantera::MultiPhase mmm;

    mmm.addPhase(salt, 10.);
    mmm.addPhase(LiCl_solid, 0.);
    mmm.addPhase(KCl_solid, 0.);


    int   solver = 2;
    int estimateEquil = 0;

    vcs_equilibrate(mmm, "TP", estimateEquil, printLvl, solver);
    cout << mmm << endl;

    delete salt;
    salt = 0;
    Cantera::appdelete();

    return retn;
  

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
    
}
