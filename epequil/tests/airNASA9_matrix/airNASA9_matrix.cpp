/*
 *  $Author: hkmoffa $
 *  $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
 *  $Revision: 508 $
 *
 *  Copyright 2002 California Institute of Technology
 *
 */
#include "cantera/IdealGasMix.h"
#include "cantera/equilibrium.h"

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

int main(int argc, char **argv) {

#ifdef DEBUG_CHEMEQUIL
  ZZCantera::ChemEquil_print_lvl = 0;
#endif

  int numSucc = 0;
  int numFail = 0;
  try {
    int retnSub;
    double T = 500.;

    IdealGasMix g("airNASA9.xml", "airNASA9");
    double pres = OneAtm;

    int kk = g.nSpecies();
    vector_fp Xmol(kk, 0.0);
    int iO2 = g.speciesIndex("O2");
    int iN2 = g.speciesIndex("N2");

    /*
     * Do an initial calculation that can be debugged
     * easily
     */
    Xmol[iO2] = 0.6;
    Xmol[iN2] = 0.4;
    T = 6100.;
    pres = 0.01;
    g.setState_TPX(T, pres, DATA_PTR(Xmol));
    try {
      retnSub = equilibrate(g, "TP", -1);
      if (retnSub != 1) {
	cerr << "ERROR: ChemEquil equilibration step failed at " 
	     << " T    = " << T 
	     << " Pres = " << pres 
             << " X_O2 = " << Xmol[iO2]
	     << endl;
	cout << "ERROR: ChemEquil equilibration step failed at " 
	     << " T    = " << T 
	     << " Pres = " << pres
             << " X_O2 = " << Xmol[iO2]
	     << endl;
        cout << g;
	exit(-1);
      } else {
        cout << g;
      }
    } catch (CanteraError) {
      showErrors(cerr);
      cerr << "ERROR: equilibration step failed at " 
	   << " T    = " << T 
	   << " Pres = " << pres 
             << " X_O2 = " << Xmol[iO2]
	   << endl;
      cout << "ERROR: equilibration step failed at " 
	   << " T    = " << T 
	   << " Pres = " << pres
             << " X_O2 = " << Xmol[iO2]
	   << endl;
      cout << g;
      exit(-1);
    }

    /*
     * Do an initial calculation that can be debugged
     * easily
     */
    Xmol[iO2] = 0.0;
    Xmol[iN2] = 1.0;
    g.setState_TPX(T, pres, DATA_PTR(Xmol));
    try {
      retnSub = equilibrate(g, "TP", -1);
      if (retnSub != 1) {
	cerr << "ERROR: ChemEquil equilibration step failed at " 
	     << " T    = " << T 
	     << " Pres = " << pres 
             << " X_O2 = " << Xmol[iO2]
	     << endl;
	cout << "ERROR: ChemEquil equilibration step failed at " 
	     << " T    = " << T 
	     << " Pres = " << pres
             << " X_O2 = " << Xmol[iO2]
	     << endl;
        cout << g;
	exit(-1);
      } else {
        cout << g;
      }
    } catch (CanteraError) {
      showErrors(cerr);
      cerr << "ERROR: equilibration step failed at " 
	   << " T    = " << T 
	   << " Pres = " << pres 
           << " X_O2 = " << Xmol[iO2]
	   << endl;
      cout << "ERROR: equilibration step failed at " 
	   << " T    = " << T 
	   << " Pres = " << pres
           << " X_O2 = " << Xmol[iO2]
	   << endl;
      cout << g;
      exit(-1);
    }

    /*
     * OK do the matrix.
     */
    bool showResults = false;
    for (int iS = 0; iS < 6; iS++) {
      Xmol[iO2] = 0.0 + 0.1 * iS;
      Xmol[iN2] = 1.0 - Xmol[iO2];
      for (int iP = -4; iP < 10; iP++) {
	pres = pow(10.0, iP) *1.0E-2;
	for (int iT = 0; iT < 45; iT++) {
	  double T = 200. * iT + 300.; 

	  /*
	   * Reset mole fractions to a base state
	   */
	  g.setState_TPX(T, pres, DATA_PTR(Xmol));

	  retnSub = 0;
	  try {
	    retnSub = equilibrate(g, "TP", -1);
	    if (retnSub != 1) {
	      cerr << "ERROR: ChemEquil equilibration step failed at " 
		   << " T    = " << T 
		   << " Pres = " << pres 
                   << " x_O2 = " << Xmol[iO2]
		   << endl;
	      cout << "ERROR: ChemEquil equilibration step failed at " 
		   << " T    = " << T 
		   << " Pres = " << pres
                   << " X_O2 = " << Xmol[iO2]
		   << endl;
	    } else {
              cout << "PASS: ChemEquil equilibration passed at "
                   << " T    = " << T
                   << " Pres = " << pres
                   << " x_O2 = " << Xmol[iO2]
                   << endl;
            }
	  } catch (CanteraError) {
	    showErrors(cerr);
	    cerr << "ERROR: equilibration step failed at " 
		 << " T    = " << T 
		 << " Pres = " << pres 
                 << " x_O2 = " << Xmol[iO2]
		 << endl;
	    cout << "ERROR: equilibration step failed at " 
		 << " T    = " << T 
		 << " Pres = " << pres
                 << " x_O2 = " << Xmol[iO2]
		 << endl;
	  }
	  if (showResults) {
	    cout << g;
	  }
	  if (retnSub ==1) {
	    numSucc++;
	  } else {
	    numFail++;
	  }

	}
      }
    }
    
    cout << "NUMBER OF SUCCESSES =  " << numSucc << endl;
    cout << "NUMBER OF FAILURES  =  " << numFail << endl;

    return numFail;
  }
  catch (CanteraError) {
    showErrors(cerr);
    cerr << "ERROR: program terminating due to unforeseen circumstances." << endl;
    return -1;
  }
}
