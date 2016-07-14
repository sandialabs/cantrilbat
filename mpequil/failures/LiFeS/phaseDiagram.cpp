///////////////////////////////////////////////////////////////////////
//
//     This program was designed to examine the phase diagram of 
//       Li[x]FeS2 as a funcion of temperature.  However, multiple 
//       errors occur for different cases as described on L74 below.
//                                            - CAL  6/29/2010
//
///////////////////////////////////////////////////////////////////////


// Include cantera header files. They should be included in the form
// These headers are designed for use in C++ programs and provide a
// simplified interface to the Cantera kernel header files. If you
// need to include kernel headers directly, use the format
// <cantera/kernel/*.h>.  /* n */

#include <cantera/thermo.h>
#include <cantera/thermo/MargulesVPSSTP.h>

#include "cantera/equilibrium.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/base/logger.h"

#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;
// All Cantera kernel names are in namespace Cantera. You can either
// reference everything as ZZCantera::<name>, or include the following
// 'using namespace' line.
#ifdef useZuzaxNamespace
using namespace Zuzax;
#define ZZCantera Zuzax
#else
using namespace Cantera;
#define ZZCantera Cantera
#endif


// The program is put into a function so that error handling code can
// be conveniently put around the whole thing. See main() below.

void createPhaseDiagram() {

  // name output files
  ofstream phaseFile;
  ofstream gasPhaseFile;
  ofstream Fe1_xSPhaseFile;
  ofstream Li2_xFe1_xS2PhaseFile;
  int tab = 20;
  int nsp = 3;
  int printLvl = 0;     
  int outputLvl = 1;

  if (outputLvl) {
    phaseFile.open("phaseFile.dat");
    gasPhaseFile.open("gasPhaseFile.dat");
    Fe1_xSPhaseFile.open("Fe1_xSPhaseFile.dat");
    Li2_xFe1_xS2PhaseFile.open("Li2_xFe1_xS2PhaseFile.dat");
    phaseFile.precision(5);
    gasPhaseFile.precision(5);
    Fe1_xSPhaseFile.precision(5);
    Li2_xFe1_xS2PhaseFile.precision(5);
  }

  // create phses to be added to multiphase
  ThermoPhase *S_solid = newPhase("S_Solid.xml","S(S)");
  ThermoPhase *S_liquid = newPhase("S_Liquid.xml","S(L)");
  ThermoPhase *gas = newPhase("Cathode_Gas.xml","Cathode(G)");
  ThermoPhase *FeS2 = newPhase("FeS2.xml","FeS2(S)");
  MargulesVPSSTP *Fe1_xS = new MargulesVPSSTP( "Fe1_xS.xml" );
  ThermoPhase *Li3Fe2S4 = newPhase("Li3Fe2S4.xml","Li3Fe2S4(S)");
  MargulesVPSSTP *Li2_xFe1_xS2;
  //////////////////////////////////////////////////////////////
  /*
    Two different failure cases
    Argument of 1 will create a "vcs_RxnStepSizes:: we shouldn't 
      be here!" error 
    Argument of 0 will converge to the wrong solution for 
      xN[iLi] = 0.3012.  By examining the ouput in phaseFile the 
      Gibbs energy at this Li concentration the Gibbs energy 
      0.44898*(-1034.6)+0.10204*(-86.298) = -473.  However, if 
      one was to take 0.05102 moles of Li3Fe2S4 and create a mixture 
      of Fe[1-x]S and Li[2+x]Fe[1-y]S2, the Gibbs energy would be -528.
  */    
  //////////////////////////////////////////////////////////////
  int argument = 0;

  if (argument){
    Li2_xFe1_xS2 = new MargulesVPSSTP("Li2_xFe1_xS2_FAILURE.xml");
  }
  else {
    Li2_xFe1_xS2 = new MargulesVPSSTP("Li2_xFe1_xS2_CONVERGENCE.xml");
  }

  int iS = 1;
  int iFe = 0;
  int iLi = 2;

  double pres = OneAtm;
  double temp;
  double xN[nsp];

  double xmin = 0;
  double xmax = 0.4;
  double tmin = 523.15;
  double tmax = 1073.15;
  int Nsmax = 50;
  int Ntmax = 11;

  // populate multiphase
  MultiPhase mmm;
  mmm.addPhase(S_solid,0.);
  mmm.addPhase(S_liquid,0.);
  mmm.addPhase(gas,0.);
  mmm.addPhase(FeS2,1.);
  mmm.addPhase(Fe1_xS,0);
  mmm.addPhase(Li3Fe2S4,0);
  mmm.addPhase(Li2_xFe1_xS2,0);
  mmm.init();
  int nPhases = mmm.nPhases();

  int iS_solid = 0;
  int iS_liquid = 1;
  int igas = 2;
  int iFeS2 = 3;
  int iFe1_xS = 4;
  int iLi3Fe2S4 = 5;
  int iLi2_xFe1_xS2 = 6;

  for ( int t = 9; t < Ntmax; t++ ){
    for ( int s = 1; s < Nsmax; s++ ){

      if (Ntmax > 1){
	temp = (tmin + (tmax-tmin)*t/(Ntmax-1));
      }
      else{
	temp = tmin;
      }
      if (Nsmax > 1){
	xN[iLi] = (xmin + (xmax-xmin)*s/(Nsmax-1));
      }
      else{
	xN[iLi] = xmin;
      }
      xN[iFe] = (1 - xN[iLi])/3;
      xN[iS] = 2*xN[iFe];

      // set phases to proper mole numbers
      double a = 2.5*xN[iLi];///na;
      double b = 1-a;//(na+1)*a)/nb;
      mmm.setPhaseMoles(iS_solid,0);
      mmm.setPhaseMoles(iS_liquid,0);
      mmm.setPhaseMoles(igas,0);
      mmm.setPhaseMoles(iFeS2,b);
      mmm.setPhaseMoles(iFe1_xS,0);
      mmm.setPhaseMoles(iLi3Fe2S4,0);
      mmm.setPhaseMoles(iLi2_xFe1_xS2,a);
      // set species mole fractions for margules phases
      double junk[2] = {1,0};
      mmm.setPhaseMoleFractions(iFe1_xS,junk);
      double crap[2] = {1,0};
      mmm.setPhaseMoleFractions(iLi2_xFe1_xS2,crap);

      mmm.setState_TP(temp,pres);

      int   solver = 2;
      int estimateEquil = 0;

      // equilibrate
      if (t == 9 && s == 1) {
	printLvl = 10;
      } else {
	printLvl = 1;
      }
      vcs_equilibrate(mmm, "TP", estimateEquil, printLvl, solver);
    
      gas = &mmm.phase(igas);

      double mu[10];
      phaseFile << endl;
      phaseFile << "Mole Fraction Li: " << xN[iLi] << endl;
      for (int j = 0; j < nPhases; j++ ) {
	phaseFile << mmm.phase(j).name() << ": ";
	mmm.phase(j).getChemPotentials(mu);
	for (int i = 0; i < mmm.phase(j).nSpecies(); i++ ){
	  phaseFile << mu[i]/1e6 << ", ";
	}
	phaseFile << endl;
      }
      phaseFile << setw(2*tab) << " ";
      for (int j = 0; j < nPhases; j++ ) {
	phaseFile << setw(tab) << mmm.phase(j).name();
      }
      phaseFile << endl;

      // output to files
      if (outputLvl) {

	phaseFile << setw(tab) << temp << setw(tab) << pres ;
	for (int ph = 0; ph < nPhases; ph++ ){
	  phaseFile << setw(tab) << mmm.phaseMoles(ph);
	}
	phaseFile << endl;

	gasPhaseFile << setw(tab) << temp << setw(tab) << pres ;
	for (int sp = 0; sp < gas->nSpecies(); sp++ ){
	  if (mmm.phaseMoles(igas) > 1e-8){
	    gasPhaseFile << setw(tab) << gas->moleFraction(sp);
	  }
	  else{
	    gasPhaseFile << setw(tab) << 0;
	  }
	}
	gasPhaseFile << endl;
	
	Fe1_xSPhaseFile << setw(tab) << temp << setw(tab) << pres ;
	for (int sp = 0; sp < Fe1_xS->nSpecies(); sp++ ){
	  if (mmm.phaseMoles(iFe1_xS) > 1e-8){
	    Fe1_xSPhaseFile << setw(tab) << Fe1_xS->moleFraction(sp);
	  }
	  else{
	    Fe1_xSPhaseFile << setw(tab) << 0;
	  }
	}
	Fe1_xSPhaseFile << endl;

	Li2_xFe1_xS2PhaseFile << setw(tab) << temp << setw(tab) << pres ;
	for (int sp = 0; sp < Li2_xFe1_xS2->nSpecies(); sp++ ){
	  if (mmm.phaseMoles(iLi2_xFe1_xS2) > 1e-8){
	    Li2_xFe1_xS2PhaseFile << setw(tab) << Li2_xFe1_xS2->moleFraction(sp);
	  }
	  else{
	    Li2_xFe1_xS2PhaseFile << setw(tab) << 0;
	  }
	}
	Li2_xFe1_xS2PhaseFile << endl;
	
      }
      
    }
  }

}

int main() {

  try {
    createPhaseDiagram();
    ZZCantera::appdelete();
  }
  catch (CanteraError) {
    showErrors(std::cout);
  }

}
