/**
 *  @file example2.cpp
 *
 * $Id: electrodeCell_kin.cpp 504 2013-01-07 22:32:48Z hkmoffa $
 * 
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

//  Example 2
//
//  Read a mechanism, and print to the standard output stream a
//  well-formatted Chemkin ELEMENT section.
//

//#include "IdealReactingGas.h"

//#include "TemperatureTable.h"

#include "electrodeCell_kin.h"

#include "LE_PickList.h"
#include "BE_MoleComp.h"

#include "mdp_allo.h"

#include "PhaseList.h"

//#include "IdealSolidSolnPhase.h"


#include "cantera/numerics/DenseMatrix.h"

// Kinetics includes
#include "cantera/kinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/thermo/SurfPhase.h"
//#include "SolidKinetics.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "ApplBase_print.h"

#include "Electrode.h"
#include "Electrode_input.h"

#include <iostream>
#include <new>
#include <string>
#include <vector>
#include <typeinfo>



using namespace Cantera;
using namespace std;
using namespace ca_ab;

#  define MIN(x,y)     (( (x) < (y) ) ? (x) : (y))


void
setAllBathSpeciesConditions(Electrode *electrode) {

  double *electrodeMF = new double[electrode->nSpecies()];
  electrode->getMoleFractions(electrodeMF);

  for (size_t iph = 0; iph < electrode->nPhases(); iph++) {
 
    ThermoPhase *tphase = &(electrode->thermo(iph));

    size_t kstart =  electrode->getGlobalSpeciesIndex(iph,0);
    double *xmol = &(electrodeMF[kstart]);
  /*
   * We need to set the state here. Note, mass fractions
   * don't matter, since we are only querying the
   * standard state of the species.
   */
    tphase->setState_TPX(electrode->temperature(), electrode->pressure(), xmol);

  //double pot = gThermoMainPhase->potential();
  /*
   * Set the electric potential
   */
  //gThermoMainPhase->setElectricPotential(electrode->);
  }
  delete [] electrodeMF;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/*
 * Calculate the arrhenius parameters from a fit from two points.
 * Ea/R is returned. This has units of Kelvin.
 *
 */
void calcArrhParam(double k1, double T1,
		   double k2, double T2,
		   double &Afac, double &EadivR) {
    if (T1 <= 0.0 || T2 <= 0.0) {
      Afac = k1;
      EadivR = 0.0;
      return;
    }
    if (T1 == T2) {
      Afac = k1;
      EadivR = 0.0;
      return;
    }
    if (k1 <= 0.0 || k2 <= 0.0) {
      Afac = k1;
      EadivR = 0.0;
      return;
    }
    EadivR = log(k1/k2) / (1.0/T2 - 1.0/T1);
    Afac = exp(log(k1) + EadivR/T1);
    return;
}

//!  Fit two current vs voltage points to an exponential relation
//!  based on the overpotential
/*!
 *   if iform = 1, we solve for the anodic branch
 *
 *     i = +- i0 exp (alpha_a F nu / RT)
 *
 *         where nu = V - Erxn
 *               F = faraday's constant
 *               R = Gas Constant
 *               T = temperature
 *     -> solve for i0 and alpha_a
 *
 *  if iform = 0, we solve for the cathodic branch
 *
 *     i = +- i0 exp (- alpha_c F nu / RT)
 *
 *         where nu = V - Erxn
 *               F = faraday's constant
 *               R = Gas Constant
 *               T = temperature
 *     -> solve for i0 and alpha_c
 *
 * Currently, the signs are passed through unchanged
 */
void calcBF(double i1, double V1, double i2, 
	    double V2, double Erxn, int iform, double FdRT, 
	    double &i0, double &alpha) {

  double n1 = V1 - Erxn;
  double n2 = V2 - Erxn;
  double sgn = 1.0;
  AssertTrace(iform == 1 || iform == 0);

  if (iform == 1) {
    // anodic branch
    if (i1 < 0.0) {
      if (i2 < 0.0) {
        sgn = -1.0;
        i1 = -i1;
        i2 = -i2;
      } else {
        throw CanteraError("calcBF"," two currents have different signs");
      }
    }
    double alpha_prime = (log(i1) - log(i2)) / (n1 - n2);
    alpha = alpha_prime / FdRT;
    i0 = exp( log(i1) - n1 * alpha_prime);
  } else {
    // Cathodic branch
    if (i1 < 0.0) {
      if (i2 < 0.0) {
        sgn = -1.0;
        i1 = -i1;
        i2 = -i2;
      } else {
        throw CanteraError("calcBF"," two currents have different signs");
      }
    }
    double alpha_prime = - (log(i1) - log(i2)) / (n1 - n2);
    alpha = alpha_prime / FdRT;
    i0 = exp( log(i1)  + n1 * alpha_prime);
  }
  i0 = i0 * sgn;

  return;
}
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/**
 * Setup the tables for the reactions.
 */
void
getKineticsTables(TemperatureTable& TT, Electrode *electrode,
		  Cantera::Kinetics &kin,
		  DenseMatrix& kfwd_Table, 
		  DenseMatrix& krev_Table,
		  DenseMatrix& deltaG_Table,
		  DenseMatrix& deltaH_Table,
		  DenseMatrix& deltaS_Table,
		  DenseMatrix& Afwd_Table,
		  DenseMatrix& EafwddivR_Table,
		  DenseMatrix& Arev_Table,
		  DenseMatrix& EarevdivR_Table,
		  DenseMatrix& kfwdPrime_Table, 
		  DenseMatrix& krevPrime_Table) {
  int i, j, iphKin;
  int nReactions = kin.nReactions();
  int nPhases = kin.nPhases();
  double *Rarray = new double [nReactions];
  double deltaT = 0.001;
  double Afac, EadivR, csc;
  int kindex;
  //Cantera::MultiPhase *mp = electrode->m_mp;
  Cantera::PhaseList *pl = electrode;

  double *electrodeMF = new double[electrode->nSpecies()];
  electrode->getMoleFractions(electrodeMF);

  double T, Tdelta;
  for (i = 0; i < TT.size(); i++) {
    T = TT[i];
    for (iphKin = 0; iphKin < nPhases; iphKin++) {
      ThermoPhase &gRef = kin.thermo(iphKin);
      //int kstart = mp->getGlobalSpeciesIndex(&gRef);
      ThermoPhase &gRef2 = pl->thermo(iphKin);
      if (&gRef != &gRef2) {
	printf("unexpected\n");
	exit(-1);
      }

      int kstart =  electrode->getGlobalSpeciesIndex(iphKin,0);
      gRef.setState_TPX(T, electrode->pressure(), &electrodeMF[kstart]);
    }

   

    kin.getFwdRateConstants(Rarray);
    for (j = 0; j < nReactions; j++) {
      kfwd_Table(i, j) = Rarray[j];	
    }

    kin.getRevRateConstants(Rarray);
    for (j = 0; j < nReactions; j++) {
      krev_Table(i, j) = Rarray[j];	
    }

    kin.getDeltaSSGibbs(Rarray);
    for (j = 0; j < nReactions; j++) {
      deltaG_Table(i, j) = Rarray[j];	
    }

    kin.getDeltaSSEnthalpy(Rarray);
    for (j = 0; j < nReactions; j++) {
      deltaH_Table(i, j) = Rarray[j];	
    }

    kin.getDeltaSSEntropy(Rarray);
    for (j = 0; j < nReactions; j++) {
      deltaS_Table(i, j) = Rarray[j];	
    }

    Tdelta = T + deltaT;
    for (iphKin = 0; iphKin < nPhases; iphKin++) {
      ThermoPhase &gRef = kin.thermo(iphKin);
      int kstart =  electrode->getGlobalSpeciesIndex(iphKin,0);
      gRef.setState_TPX(Tdelta, electrode->pressure(), &electrodeMF[kstart]);
    }


    kin.getFwdRateConstants(Rarray);
    for (j = 0; j < nReactions; j++) {
      calcArrhParam(kfwd_Table(i, j), T,
		    Rarray[j], Tdelta,
		    Afac, EadivR);
      Afwd_Table(i,j) = Afac;
      EafwddivR_Table(i,j) = EadivR;
    }

    kin.getRevRateConstants(Rarray);
    for (j = 0; j < nReactions; j++) {
      calcArrhParam(krev_Table(i, j), T,
		    Rarray[j], Tdelta,
		    Afac, EadivR);
      Arev_Table(i,j) = Afac;
      EarevdivR_Table(i,j) = EadivR;
    }

    for (j = 0; j < nReactions; j++) {
      /*
       * reactants
       */
      const vector<size_t>& reactants = kin.reactants(j);
      int nReac = reactants.size();
      /*
       * products
       */
      const vector<size_t>& products = kin.products(j);
      int nProd = products.size();
	
      kfwdPrime_Table(i, j) = kfwd_Table(i,j);
      for (int k = 0; k < nReac; k++) {
	kindex = reactants[k];
	int iphase = kin.speciesPhaseIndex(kindex);
	ThermoPhase &tpRef = kin.thermo(iphase);

	int kLocal = kindex - kin.kineticsSpeciesIndex(0, iphase);
	csc = tpRef.standardConcentration(kLocal);
	kfwdPrime_Table(i, j) *= csc;
      }

      krevPrime_Table(i, j) = krev_Table(i,j);
      for (int k = 0; k < nProd; k++) {
	kindex = products[k];
	int iphase = kin.speciesPhaseIndex(kindex);
	ThermoPhase &tpRef = kin.thermo(iphase);
	int kLocal = kindex - kin.kineticsSpeciesIndex(0, iphase);
	csc = tpRef.standardConcentration(kLocal);
	krevPrime_Table(i, j) *= csc;
      }

    }

  }
     
  for (iphKin = 0; iphKin < nPhases; iphKin++) {
    ThermoPhase &gRef = kin.thermo(iphKin);
    int kstart =  electrode->getGlobalSpeciesIndex(iphKin,0);
    gRef.setState_TPX(Tdelta, electrode->pressure(), &electrodeMF[kstart]);
    
  
  }


   delete [] electrodeMF; 
   delete [] Rarray;
}
  /*****************************************************************/
  /*****************************************************************/
  /*****************************************************************/

  string iunit_string(int i) {
    string val;
    switch (i) {
    case 0:
      val = "kmol";
      break;
    case 1:
      val = "m";
      break;
    case 2:
      val = "Kg";
      break;
    case 3:
      val = "Pa";
      break;
    case 4:
      val = "K";
      break;
    case 5:
      val = "s";
      break;
    default:
      val = "";
      break;
    }
    return val;
  }

  bool isInteger(double num) {
    int inum = (int) num;
    double rnum = inum;
    if (rnum == num) return true;
    return false;
  }
  /*****************************************************************/
  /*****************************************************************/
  /*****************************************************************/

  string formUnitsString(double unitsA[6]) {
    string num;
    string denom;
    string val;
    int numNum = 0;
    int numDenom = 0;
    for (int i = 0; i < 6; i++) {
      if (unitsA[i] > 0.0) numNum++;
      if (unitsA[i] < 0.0) numDenom++;
    }
    if ((numNum + numDenom) == 0) return "";
    if (numNum > 0) {
      int itmp = 0;
      for (int i = 0; i < 6; i++) {
	if (unitsA[i] > 0.0) {
	  val = iunit_string(i);
	  if (itmp > 0) {
	    num = num + " ";
	  }
	  num += val;
          if (isInteger(unitsA[i])) {
	    int np = (int) unitsA[i];
	    if (np != 1) {
	      num += "^" + int2str(np);
	    }
	  } else {
	    if (unitsA[i] == 0.5) {
	      num += "^0.5";
	    } else {
	      num += "^" + fp2str(unitsA[i]);
	    }
	  }
	  itmp++;
	}
      }
    } else {
      num = "1";
    }
    if (numDenom > 0) {
      int idtmp = 0;
      for (int i = 0; i < 6; i++) {
	if (unitsA[i] < 0.0) {
	  val = iunit_string(i);
	  double dp = - unitsA[i];
	  if (idtmp > 0) {
	    denom = denom + " ";
	  }
	  denom += val;
          if (isInteger(dp)) {
	    int np = (int) dp;
	    if (np != 1) {
	      denom += "^" + int2str(np);
	    }
	  } else {
	    if (dp == 0.5) {
	      denom += "^0.5";
	    } else {
	      denom += "^" + fp2str(dp);
	    }
	  }
	  idtmp++;
	}
      }
    }
    if (numNum > 0 && numDenom > 0) {
      if (numNum > 1) {
	string tmp = "(" + num + ")";
	num = tmp;
      }
      if (numDenom > 1) {
	string tmp = "(" + denom + ")";
	denom = tmp;
      }
    }
    if (numDenom > 0) {
      val = num + "/" + denom;
    } else {
      val = num;
    }
    return val;
  }

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/**
 *  This routine will print out a table of information a single
 *  reaction, j.
 *
 *
 */
void printKineticsTable(Electrode *electrode, int j,
			TemperatureTable &TT,
			Cantera::Kinetics &kin, 
			DenseMatrix& kfwd_Table, 
			DenseMatrix& krev_Table,
			DenseMatrix& deltaG_Table,
			DenseMatrix& deltaH_Table,
			DenseMatrix& deltaS_Table,
			DenseMatrix& Afwd_Table,
			DenseMatrix& EafwddivR_Table,
			DenseMatrix& Arev_Table,
			DenseMatrix& EarevdivR_Table,
			DenseMatrix& kfwdPrime_Table, 
			DenseMatrix& krevPrime_Table,
			Cantera::RxnMolChange *rmc) {

    /*
     *  Get the species data object from the Mixture object
     *  this is defined in the constituents.h file, and is 
     *  inherited by Mixture through BaseMix
     */

    int tableWidth = 128;

    /*
     * Conversion factor between Joules/kmol, cantera's units,
     * and kcal / gmol.
     */
    double conv = 1.0 / (4.184E6);
    double convS = conv * 1.0E3;
    double Rkcal = 1.98721E-3;

    ThermoPhase *gThermo = &kin.thermo(0);
    /*
     *  Dump out all of the information about the reaction
     */
    cout << endl;
    print_char('=', tableWidth);
    cout << endl;
    cout << endl;
    cout << "INFORMATION TABLE FOR Reaction \""  << j;
    cout << "\" IN PHASE \"";
    cout << gThermo->id();
    cout << "\"" << endl;
    cout << endl;

    /*
     * Reaction String
     */
    string rs = kin.reactionString(j);
    dnt(1); 
    cout << rs << endl;
    cout << endl;

    /*
     * Print out reactants
     */
    const vector<size_t>& reactants = kin.reactants(j);
    int nReac = reactants.size();
    dnt(1); 
    if (nReac == 1) {
      cout << "There is one reactant: ";
    } else {
      cout << "There are " << nReac << " reactants: ";
    }
    for (int k = 0; k < nReac; k++) {
      int kindex = reactants[k];
      cout << kin.kineticsSpeciesName(kindex) << " ";
    }
    cout << endl;


    /*
     * Print out products
     */
    const vector<size_t>& products = kin.products(j);
    int nProd = products.size();
    dnt(1); 
    if (nProd == 1) {
      cout << "There is one product: ";
    } else {
      cout << "There are " << nProd << " products: ";
    }
    for (int k = 0; k < nProd; k++) {
      int kindex = products[k];
      cout << kin.kineticsSpeciesName(kindex) << " ";
    }
    cout << endl;

    /*
     *  reaction type
     */
    double pCurrent = gThermo->pressure();
    int rType = kin.reactionType(j);
    dnt(1); 
    cout << "Reaction Type = ";
    switch (rType) {
    case ELEMENTARY_RXN:
	cout << "Elementary Rxn" << endl;
	break;
    case THREE_BODY_RXN:

	cout << "Three Body Rxn: "<< endl;
	dnt(1); 
	cout <<"        Rate Constants below include third body collision partner "
	     << endl;
	dnt(1);
	cout <<"        and are effected by current pressure, " 
	     << pCurrent << " Pa " << endl;
	break;
    case FALLOFF_RXN:
	cout << "Falloff Rxn" << endl;
	dnt(1); 
	cout <<"        Rate Constants below include third body collision partner effects"
	     << endl;
	dnt(1);
	cout <<"        and are effected by current pressure, " 
	     << pCurrent << " Pa " << endl;
	break;
    case CHEMACT_RXN:
	cout << "Chemically activated rxn" << endl;
	dnt(1); 
	cout <<"        Rate Constants below include third body collision partner effects"
	     << endl;
	dnt(1);
	cout <<"        and are effected by current pressure, " 
	     << pCurrent << " Pa " << endl;
	break;
    case SURFACE_RXN:
	cout << "Surface rxn" << endl;
	break;
    case GLOBAL_RXN:
	cout << "Global Rxn" << endl;
	break;
    default:
	cout << "Unknown_Type - (" << rType << ")" << endl;
	break;
    }
  

    /*
     * Reaction Reversible?
     */
    bool rReversible = kin.isReversible(j);
    dnt(1); 
    if (rReversible) {
      cout << "Reaction is reversible" << endl;
    } else {
      cout << "Reaction is not reversible" << endl;
    }

    cout << endl;

    /*
     * Figure out the units for the forward reaction rate constant
     * First start out with the units for the rop
     */
    int nphase = kin.nPhases();
    int ndim = 3;
    for (int iph = 0; iph < nphase; iph++) {
      ThermoPhase &tpRef = kin.thermo(iph);
      int idim = tpRef.nDim();
      if (idim < ndim) ndim = idim;
    }
    doublereal unitsROP[6] = { 1.0,(double) (-ndim), 0.0, 0.0, 0.0, -1.0 };
    doublereal unitskfwd[6];
    doublereal unitskrev[6] ;
    doublereal unitsSpecies[6];
    for (int i = 0; i < 6; i++) {
      unitskfwd[i] = unitsROP[i];
      unitskrev[i] = unitsROP[i];
    }
    for (int k = 0; k < nReac; k++) {
      int kindex = reactants[k];
      int iphase = kin.speciesPhaseIndex(kindex);
      ThermoPhase &tpRef = kin.thermo(iphase);
      tpRef.getUnitsStandardConc(unitsSpecies, kindex);
      for (int i = 0; i < 6; i++) {
	unitskfwd[i] -= unitsSpecies[i];
      }
    }
    for (int k = 0; k < nProd; k++) {
      int kindex = products[k];
      int iphase = kin.speciesPhaseIndex(kindex);
      ThermoPhase &tpRef = kin.thermo(iphase);
      tpRef.getUnitsStandardConc(unitsSpecies, kindex);
      for (int i = 0; i < 6; i++) {
	unitskrev[i] -= unitsSpecies[i];
      }
    }
 
    /****************** phases change table ****************************/
    int pwidth = 52;
    if (rmc->m_ChargeTransferInRxn != 0.0) {
      dnt(1);
      printf("This reaction is charge transfer reaction, n = %g\n",
	     rmc->m_ChargeTransferInRxn);
      dnt(1); printf("Electrochemical Beta = %g\n", rmc->m_beta); 
      pwidth = 86;
    }
 
    dnt(2); print_char('-', pwidth); printf("\n");
 
    dnt(2); printf("|              Phase | Dim Delta_Moles  Delta_Mass |");
    if (rmc->m_ChargeTransferInRxn != 0.0) {
      printf("Electron_Transfer Base_Potential |");
    }
    printf("\n");
    for (int iph = 0; iph < rmc->m_nPhases; iph++) {
      ThermoPhase &tpRef = kin.thermo(iph);
      string sname = tpRef.name();
      dnt(2);
      printf("|%20s|", sname.c_str()); 
      printf("%4d",    rmc->m_phaseDims[iph]);
      printf("%12.5g", rmc->m_phaseMoleChange[iph]);
      printf("%12.5g", rmc->m_phaseMassChange[iph]);
      if (rmc->m_ChargeTransferInRxn != 0.0) {
	printf(" |");
	printf("  %12.5g", rmc->m_phaseChargeChange[iph]);
        double phi = tpRef.electricPotential();
	printf("      %12.5g", phi);
      }
      printf(" |\n");
    }
    
    dnt(2);  print_char('-', pwidth); printf("\n\n");
 

    /*******************************************************************/


    string units_kfwd = formUnitsString(unitskfwd);
    string units_krev = formUnitsString(unitskrev);
    dnt(1);
    cout << "units_kfwd = Units for forward reaction rate constant = " 
	 << units_kfwd << endl;
    dnt(1);
    cout << "units_krev = Units for reverse reaction rate constant = " 
	 << units_krev << endl;

    /*
     * Dump out the kinetics table as a function of the
     * temperature
     */
    cout << endl;
    cout << "|";
    print_char('-', tableWidth); cout << "|" << endl;
    cout << "|        |                              "
	 << " |                                 ";
    if (rReversible) {
      cout << "|                               ";
    } else {
      cout << "| -- REVERSE RXN NOT IN MECH -- ";
    }
    cout << "|                     ";
    cout << "|" << endl;
    cout << "|  Temp  |   kfwd      Afwd     Ea_fwd  " 
	 << " |  DeltaG0    DeltaH0    DeltaS0  "
	 << "|  k_rev     A_rev     Ea_rev   "
	 << "|  kstar   kstar_rev  ";
 
    cout << "|" << endl;
    cout << "|   (K)  |   (units_kfwd)     (kcal/gmol)"
	 << "|(kcal/gmol)(kcal/gmol)(cal/gmolK)"
	 << "|  (units_krev)      (kcal/gmol)|"
	 << "    (kmol / m^" << ndim << " s)   |";
    cout << endl;
    cout << "|--------|----------"
	 << "---------------------|-----------------"
	 << "----------------|-----------------------------"
	 << "--|---------------------|" << endl;
  
    for (int i = 0; i < TT.size(); i++) {
      cout << "|";
      pr_df(TT[i], 7, 1);
      cout << " |";
      pr_de(kfwd_Table(i,j), 10, 2);
      pr_de(Afwd_Table(i,j), 10, 2);
      pr_dg(EafwddivR_Table(i,j) * Rkcal, 10, 2);
      cout << " | ";
      pr_dg(deltaG_Table(i,j) * conv, 10, 2);
      pr_dg(deltaH_Table(i,j) * conv, 10, 2);
      pr_dg(deltaS_Table(i,j) * convS, 10, 2);
      cout << "  |";
      pr_de(krev_Table(i,j), 10, 2);
      pr_de(Arev_Table(i,j), 10, 2);
      pr_dg(EarevdivR_Table(i,j) * Rkcal, 10, 2);
      cout << " |";
      pr_de(kfwdPrime_Table(i,j), 10, 2);
      pr_de(krevPrime_Table(i,j), 10, 2);
      cout << " |" << endl;
    }
    cout << "|" ;
    print_char('-', tableWidth);
    cout << "|" << endl;
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

void doKineticsTablesHetero(Electrode *electrode, 
			    InterfaceKinetics *gKinetics,   
			    TemperatureTable &TT) {
    
    int nReactions = gKinetics->nReactions();
    // int nPhases = gKinetics->nPhases();
    DenseMatrix kfwd_Table(TT.size(), nReactions);
    DenseMatrix krev_Table(TT.size(), nReactions);
    DenseMatrix deltaG_Table(TT.size(), nReactions);
    DenseMatrix deltaH_Table(TT.size(), nReactions);
    DenseMatrix deltaS_Table(TT.size(), nReactions);
    DenseMatrix Afwd_Table(TT.size(), nReactions);
    DenseMatrix EafwddivR_Table(TT.size(), nReactions);
    DenseMatrix Arev_Table(TT.size(), nReactions);
    DenseMatrix EarevdivR_Table(TT.size(), nReactions);
    DenseMatrix kfwdPrime_Table(TT.size(), nReactions);
    DenseMatrix krevPrime_Table(TT.size(), nReactions);


    getKineticsTables(TT, electrode,
		      *gKinetics, kfwd_Table, krev_Table,
		      deltaG_Table, deltaH_Table, deltaS_Table,
		      Afwd_Table, EafwddivR_Table, 
		      Arev_Table, EarevdivR_Table,
		      kfwdPrime_Table, krevPrime_Table);

    std::vector<Cantera::RxnMolChange *> rmcVector;
    rmcVector.resize(nReactions,0);

    for (int i = 0; i < nReactions; i++) {
      rmcVector[i] = new Cantera::RxnMolChange(gKinetics, i);

      printKineticsTable(electrode, i, TT,
			 *gKinetics, kfwd_Table, krev_Table,
			 deltaG_Table, deltaH_Table, deltaS_Table,
			 Afwd_Table, EafwddivR_Table, 
			 Arev_Table, EarevdivR_Table,
			 kfwdPrime_Table, krevPrime_Table,
			 rmcVector[i]);

      if ((rmcVector[i])->m_ChargeTransferInRxn != 0.0) {
	processCurrentVsPotTable(rmcVector[i],
				 electrode, i, TT,
				 *gKinetics, kfwd_Table, krev_Table,
				 deltaG_Table, deltaH_Table, deltaS_Table,
				 Afwd_Table, EafwddivR_Table, 
				 Arev_Table, EarevdivR_Table,
				 kfwdPrime_Table, krevPrime_Table);
      }

    }
    setAllBathSpeciesConditions(electrode);

    int  numExtraGlobalRxns = electrode->numExtraGlobalRxnPathways();
    for (int iextra = 0; iextra < numExtraGlobalRxns; iextra++) {
      /*
      struct EGRInput * egr_ptr = electrode->m_EGRList[iextra];
      ExtraGlobalRxn *egr = new ExtraGlobalRxn(gKinetics);
      double *RxnVector  = new double[nReactions];
      for (int i = 0; i < nReactions; i++) {
	RxnVector[i] = 0.0;
      }
      for (int iErxn = 0; iErxn < egr_ptr->m_numElemReactions; iErxn++) {
	struct ERSSpec * ers_ptr = egr_ptr->m_ERSList[iErxn];
	RxnVector[ers_ptr->m_reactionIndex] = ers_ptr->m_reactionMultiplier;
      }
      egr->setupElemRxnVector(RxnVector, 
			      egr_ptr->m_SS_KinSpeciesKindex);
      RxnMolChange *rmcEGR = new RxnMolChange(gKinetics, egr);
      */

      ExtraGlobalRxn *egr = electrode->extraGlobalRxnPathway(iextra);
      Cantera::RxnMolChange *rmcEGR =  electrode->rxnMolChangesEGR(iextra);
      RxnTempTableStuff * rts = new RxnTempTableStuff(-1,0);
      getGERKineticsTables(TT, electrode, *gKinetics, *egr, *rts); 
      setAllBathSpeciesConditions(electrode);

      printGERKineticsTable(electrode, -1, TT, *gKinetics, *egr, rmcEGR, *rts);

      if ((rmcEGR)->m_ChargeTransferInRxn != 0.0) {
	processGERCurrentVsPotTable(rmcEGR, electrode, 0, TT,
				    *gKinetics, *egr, *rts);
      }

      //delete egr; egr=0;
      //delete [] RxnVector;
      //delete rmcEGR; rmcEGR = 0;
      delete rts; rts=0;
    }
  

    for (int i = 0; i < nReactions; i++) {
      delete  rmcVector[i];
    }
}


void doKineticsTablesHomog(Electrode *electrode, Cantera::Kinetics *gKinetics,   
			   TemperatureTable &TT) {
    
  int nReactions = gKinetics->nReactions();
  DenseMatrix kfwd_Table(TT.size(), nReactions);
  DenseMatrix krev_Table(TT.size(), nReactions);
  DenseMatrix deltaG_Table(TT.size(), nReactions);
  DenseMatrix deltaH_Table(TT.size(), nReactions);
  DenseMatrix deltaS_Table(TT.size(), nReactions);
  DenseMatrix Afwd_Table(TT.size(), nReactions);
  DenseMatrix EafwddivR_Table(TT.size(), nReactions);
  DenseMatrix Arev_Table(TT.size(), nReactions);
  DenseMatrix EarevdivR_Table(TT.size(), nReactions);
  DenseMatrix kfwdPrime_Table(TT.size(), nReactions);
  DenseMatrix krevPrime_Table(TT.size(), nReactions);

  getKineticsTables(TT, electrode,
		    *gKinetics, kfwd_Table, krev_Table,
		    deltaG_Table, deltaH_Table, deltaS_Table,
		    Afwd_Table, EafwddivR_Table, 
		    Arev_Table, EarevdivR_Table,
		    kfwdPrime_Table, krevPrime_Table);
  vector<Cantera::RxnMolChange *> rmcVector;
  rmcVector.resize(nReactions,0);
  for (int i = 0; i < nReactions; i++) {
    rmcVector[i] = new Cantera::RxnMolChange(gKinetics, i);
    printKineticsTable(electrode, i, TT,
		       *gKinetics, kfwd_Table, krev_Table,
		       deltaG_Table, deltaH_Table, deltaS_Table,
		       Afwd_Table, EafwddivR_Table, 
		       Arev_Table, EarevdivR_Table,
		       kfwdPrime_Table, krevPrime_Table,	 rmcVector[i]);
  }

  for (int i = 0; i < nReactions; i++) {
    delete  rmcVector[i];
  }
}


/*************************************************************************/
void processCurrentVsPotTable(Cantera::RxnMolChange *rmc,
			      Electrode *electrode, int irxn,
			      TemperatureTable &TT,
			      Kinetics &kin, 
			      DenseMatrix& kfwd_Table, 
			      DenseMatrix& krev_Table,
			      DenseMatrix& deltaG_Table,
			      DenseMatrix& deltaH_Table,
			      DenseMatrix& deltaS_Table,
			      DenseMatrix& Afwd_Table,
			      DenseMatrix& EafwddivR_Table,
			      DenseMatrix& Arev_Table,
			      DenseMatrix& EarevdivR_Table,
			      DenseMatrix& kfwdPrime_Table, 
			      DenseMatrix& krevPrime_Table)
{
  int iph, k, e;
  ThermoPhase *tp = 0;
  ThermoPhase *tpMetal = 0;
  int nTotalSpecies = kin.nTotalSpecies();
  int nRxns         = kin.nReactions();
  vector<double> Rfwd(nTotalSpecies, 0.0);
  vector<double> Rrev(nTotalSpecies, 0.0);
  vector<double> Rnet(nTotalSpecies, 0.0);
  vector<double> deltaG(nRxns, 0.0);
  vector<double> deltaSSG(nRxns, 0.0);

  InterfaceKinetics *iK = dynamic_cast<InterfaceKinetics *>(&kin);
  if (!iK) {
    throw CanteraError("RxnMolChange", "unknown condition on charge");
  }
  /*
   * First thing is to figure out how to calculate the voltage
   *    V = phi_metal - phi_soln
   * Look for the electron species. The phase where this occurs will be
   * called the metal, iMetal;
   */
  int iMetal = 0;
  int iSoln = 0;
  int kElectron = -1;
  int nPhases = iK->nPhases();
  for (iph = 0; iph < nPhases; iph++) {
    tp = &(iK->thermo(iph));
    int nSpecies = tp->nSpecies();
    int nElements = tp->nElements();
    int eElectron = tp->elementIndex("E");
    if (eElectron >= 0) {
      for (k = 0; k < nSpecies; k++) {
	if (tp->nAtoms(k,eElectron) == 1) {
	  int ifound = 1;
	  for (e = 0; e < nElements; e++) {
	    if (tp->nAtoms(k,e) != 0.0) {
	      if (e != eElectron) {
		ifound = 0;
	      }
	    }
	  }
	  if (ifound == 1) {
	    iMetal = iph;
	    kElectron = iK->kineticsSpeciesIndex(k, iph);
	  }
	}
      }
    }
    if (iph != iMetal) {
      if (rmc->m_phaseChargeChange[iph] != 0) {
	iSoln = iph;
      }
    }
  }
  int iElectrodeMetal = electrode->metalPhaseIndex();
  int iElectrodeLyte  = electrode->solnPhaseIndex();

  double phi0Metal = electrode->phaseVoltage(iElectrodeMetal);
  //double phi0Metal = rmc->m_phasePotentials[iMetal];
  double phi0Soln = electrode->phaseVoltage(iElectrodeLyte);
  double nStoichElectrons = - rmc->m_phaseChargeChange[iMetal];
  double V0 = phi0Metal  -phi0Soln;
  tpMetal = &(iK->thermo(iMetal));
  ThermoPhase *tpSoln =  &(iK->thermo(iSoln));
  int iSurf = iK->reactionPhaseIndex();
  ThermoPhase *tps = &(iK->thermo(iSurf));
  SurfPhase *tpSurface = dynamic_cast<SurfPhase *>(tps);

  int nSurfSpecies = tpSurface->nSpecies();
  vector<double> thetaSurf(nSurfSpecies, 0.0);
  tpSurface->getCoverages(DATA_PTR(thetaSurf));

  int pwidth = 100;
  int twidth = 93;
  dnt(1); print_char('-', pwidth); printf("\n");
  string rs = kin.reactionString(irxn);
  dnt(1); printf("Voltage Dependence for Rxn: %s", rs.c_str());
  printf("\n\n");
  dnt(2); printf("ThermoPhase index pertaining to the solution = %d\n", iSoln);
  dnt(2); printf("ThermoPhase index pertaining to the metal = %d\n", iMetal); 
  dnt(2); printf("Kinetic Species Index pertaining to the electrons = %d\n", kElectron);
  dnt(2); printf("Input Electric Potential of the Metal = %g Volts\n", phi0Metal);
  dnt(2); printf("Input Electric Potential of the Soln  = %g Volts\n", phi0Soln);
  dnt(2); printf("Input Voltage = phi_Metal - phi_Soln  = %g Volts\n", V0);
  dnt(2); printf("nStoichElectrons = %g \n", nStoichElectrons); 
  /*
   * the 1.0E-4 is to change from m-2 to cm-2.
   */
  iK->getDeltaGibbs(DATA_PTR(deltaG));
  iK->getDeltaSSGibbs(DATA_PTR(deltaSSG));
  double deltaGrxn = deltaG[irxn];
  double deltaSSGrxn = deltaSSG[irxn];
  double Erxn =  deltaGrxn/ Faraday/ nStoichElectrons;
  double ESSrxn =  deltaSSGrxn/ Faraday/  nStoichElectrons;
  double nElectronProducts = iK->productStoichCoeff(kElectron, irxn);
  double nElectronReactants = iK->reactantStoichCoeff(kElectron, irxn);
  bool AnodicDir = false;
  if (nElectronProducts > nElectronReactants) {
    AnodicDir = true;
    Erxn = - Erxn;
    ESSrxn = - ESSrxn;
  }
  dnt(2); printf("Delta G (rxn)   = %13.5g J/kmol   --->   ", deltaGrxn);
  dnt(2); printf("E    = %g Volts", Erxn);
  if (AnodicDir) {
    printf(" (Switching to electrons as reactants)");
  }
  printf("\n");
  dnt(2); printf("Delta G_SS (rxn)= %13.5g J/kmol   --->   ", deltaSSGrxn);
  dnt(2); printf("E_SS = %g Volts", ESSrxn);
  if (AnodicDir) {
    printf(" (Switching to electrons as reactants)");
  }
  printf("\n");
  dnt(2); printf("Bath Gas conditions:\n");
  dnt(3); printf("T = %g K\n", tpSoln->temperature());
  dnt(3); printf("P = %g Pascal\n", tpSoln->pressure());
  dnt(3); printf("     SurfSpecName         Theta_k_Bath\n");
  for (k = 0; k < nSurfSpecies; k++) {
    string sName = tpSurface->speciesName(k);
    dnt(3); printf("%16s %12.5g\n", sName.c_str(), thetaSurf[k]);
  }
  printf("\n");

  double startingVoltage = -0.6;
  if (Erxn < startingVoltage) {
    startingVoltage = Erxn - 0.2;
    startingVoltage = (int(startingVoltage*10.0))/10.0;
  }

  double nF = Faraday * rmc->m_phaseChargeChange[iMetal] * 1.0E-4; 
 
  double Voltage0 = V0 + startingVoltage;
  double phiMetal0 = phi0Metal + startingVoltage;
  dnt(2); print_char('-', twidth); printf("\n");
  dnt(2);
  cout << "|  Voltage  |    Rfwd          Rrev         Rnet   " 
       << "|        Ifwd         Irev         Inet  ";
  cout << "|" << endl;

  dnt(2);
  cout << "|   Volts   |             (kmol/m**2 s)            "
       << "|               (amps/cm2)               "
       << "|" << endl;
  dnt(2); print_char('-', twidth); printf("\n");
  bool iEdone = false;
  double Voltage, phiMetal;
  for (int iV = 0; iV < 12; iV++) {
    Voltage0 += 0.1;
    phiMetal0 += 0.1;
    
    if (!iEdone && (Voltage0 > Erxn)) {
      iEdone = true;
      Voltage = Erxn;
      ThermoPhase &thRef = electrode->thermo(iElectrodeMetal);
      phiMetal =  Erxn + thRef.electricPotential();
     // phiMetal = Erxn + rmc->m_phasePotentials[iSoln];
    } else {
      Voltage = Voltage0;
      phiMetal = phiMetal0;
    }
    if (fabs(Voltage) < 1.0E-10) Voltage = 0.0;
    tpMetal->setElectricPotential(phiMetal);

    iK->getFwdRatesOfProgress(DATA_PTR(Rfwd));
    iK->getRevRatesOfProgress(DATA_PTR(Rrev));
    iK->getNetRatesOfProgress(DATA_PTR(Rnet));
    dnt(2);
    printf("!%10.4g |%11.3g %12.3g %12.3g |%12.3g %12.3g  %12.3g |\n", 
	   Voltage, Rfwd[irxn], Rrev[irxn], Rnet[irxn], nF*Rfwd[irxn], 
	   nF*Rrev[irxn], nF*Rnet[irxn] );

  }
  dnt(2); print_char('-', twidth); printf("\n");
  
  dnt(2); printf(" TABLE that uses the Surface Solver to "
		 "Solve for Surface Site Concentrations\n");
 
  Voltage0 = V0 - 0.6;
  phiMetal0 = phi0Metal - 0.6;
  int nSS = MIN(nSurfSpecies, 6);
  twidth = 54 + nSS * 11;
  dnt(2); print_char('-', twidth); printf("\n");
  dnt(2);
  cout << "|  Voltage  " 
       << "|       Ifwd         Irev        Inet  ";
  cout << "| ";

  for (k = 0; k < nSS; k++) {
    string sname =  tpSurface->speciesName(k);
    printf("%10s ", sname.c_str());
  }
  cout <<"|" << endl;
  dnt(2);
  cout << "|   Volts   "
       << "|               (amps/cm2)             "
       << "| " ;
  for (k = 0; k < nSS; k++) {
    printf("           ");
  }
  cout <<"|" << endl;
  dnt(2); print_char('-', twidth); printf("\n");
  iEdone = false;

  for (int iV = 0; iV < 12; iV++) {
    Voltage0 += 0.1;
    phiMetal0 += 0.1;
    
    if (!iEdone && (Voltage0 > Erxn)) {
      iEdone = true;
      Voltage = Erxn;
     // phiMetal = Erxn + rmc->m_phasePotentials[iSoln];
      ThermoPhase &thRef = electrode->thermo(iElectrodeMetal);
      phiMetal =  Erxn + thRef.electricPotential();

    } else {
      Voltage = Voltage0;
      phiMetal = phiMetal0;
    }
    if (fabs(Voltage) < 1.0E-10) Voltage = 0.0;
    tpMetal->setElectricPotential(phiMetal);
    iK->advanceCoverages(100.0);
    iK->getFwdRatesOfProgress(DATA_PTR(Rfwd));
    iK->getRevRatesOfProgress(DATA_PTR(Rrev));
    iK->getNetRatesOfProgress(DATA_PTR(Rnet));
    
    dnt(2);
    printf("!%10.4g |%11.3g %12.3g %12.3g |", 
	   Voltage,
	   nF*Rfwd[irxn], nF*Rrev[irxn], nF*Rnet[irxn] );
    tpSurface->getCoverages(DATA_PTR(thetaSurf));
    for (k = 0; k < nSS; k++) {
      printf("%10.4g ", thetaSurf[k]);
    }
    printf(" |\n");
  }
  dnt(2); print_char('-', twidth); printf("\n");
  


  dnt(1); print_char('-', pwidth); printf("\n");
}


RxnTempTableStuff::RxnTempTableStuff(int irxn, int irxnGE) :
  m_irxn(irxn),
  m_irxnGlobalExtra(irxnGE)
 {
 }

RxnTempTableStuff::~RxnTempTableStuff() {
}

/*
 *
 */
void
getGERKineticsTables(TemperatureTable& TT, Electrode *electrode,
		     Cantera::Kinetics &kin,
		     ExtraGlobalRxn &egr,
		     RxnTempTableStuff &rts) {
  

  int i, iphKin;
  int nReactions = kin.nReactions();
  int nPhases = kin.nPhases();
  int nTotSpecies = kin.nTotalSpecies();
  double *Rarray = new double [nReactions];
  double *R2array = new double [nReactions];
  double *Sarray = new double [nTotSpecies];
  double deltaT = 0.001;
  double Afac, EadivR;

  vector<double> & kfwd_Table      = rts.kfwd_Table;
  vector<double> & krev_Table      = rts.krev_Table;
  vector<double> & deltaG_Table    = rts.deltaG_Table;
  vector<double> & deltaH_Table    = rts.deltaH_Table;
  vector<double> & deltaS_Table    = rts.deltaS_Table;
  vector<double> & deltaGss_Table    = rts.deltaGss_Table;
  vector<double> & deltaHss_Table    = rts.deltaHss_Table;
  vector<double> & deltaSss_Table    = rts.deltaSss_Table;
  vector<double> & Afwd_Table      = rts.Afwd_Table;
  vector<double> & EafwddivR_Table = rts.EafwddivR_Table;
  vector<double> & Arev_Table      = rts.Arev_Table;
  vector<double> & EarevdivR_Table = rts.EarevdivR_Table;
  vector<double> & kfwdPrime_Table = rts.kfwdPrime_Table;
  vector<double> & krevPrime_Table = rts.krevPrime_Table;
  vector<double> & NetROP_Table    = rts.NetROP_Table;
  vector<double> & FwdROP_Table    = rts.FwdROP_Table;
  vector<double> & RevROP_Table    = rts.RevROP_Table;
  vector<double> & Anet_Table      = rts.Anet_Table;
  vector<double> & EanetdivR_Table = rts.EanetdivR_Table;

  int m = TT.size();
  kfwd_Table.resize(m,0.0);
  krev_Table.resize(m,0.0);
  deltaG_Table.resize(m,0.0);
  deltaH_Table.resize(m,0.0);
  deltaS_Table.resize(m,0.0);
  deltaGss_Table.resize(m,0.0);
  deltaHss_Table.resize(m,0.0);
  deltaSss_Table.resize(m,0.0);
  Afwd_Table.resize(m,0.0);
  EafwddivR_Table.resize(m,0.0);
  Arev_Table.resize(m,0.0);
  EarevdivR_Table.resize(m,0.0);
  kfwdPrime_Table.resize(m,0.0);
  krevPrime_Table.resize(m,0.0);
  NetROP_Table.resize(m,0.0);
  FwdROP_Table.resize(m,0.0);
  RevROP_Table.resize(m,0.0);
  Anet_Table.resize(m,0.0);
  EanetdivR_Table.resize(m,0.0);

  InterfaceKinetics *iKin = dynamic_cast<InterfaceKinetics *>(&kin);

  setAllBathSpeciesConditions(electrode);
  double T, Tdelta;

  double *electrodeMF = new double[electrode->nSpecies()];
  electrode->getMoleFractions(electrodeMF);
 
  for (i = 0; i < TT.size(); i++) {
    T = TT[i];
   
    electrode->getMoleFractions(electrodeMF);
    for (iphKin = 0; iphKin < nPhases; iphKin++) {
      ThermoPhase &gRef = kin.thermo(iphKin);
      int kstart =  electrode->getGlobalSpeciesIndex(iphKin,0);
      gRef.setState_TPX(T, electrode->pressure(), &(electrodeMF[kstart]));
    }


    if (iKin) {
      iKin->advanceCoverages(100.0);
    }
    
    kin.getFwdRatesOfProgress(Rarray);
    kin.getRevRatesOfProgress(R2array);
    FwdROP_Table[i] = egr.FwdROPValue(Rarray, R2array);
    RevROP_Table[i] = egr.RevROPValue(Rarray, R2array);

    kin.getNetRatesOfProgress(Rarray);
    NetROP_Table[i] = egr.ROPValue(Rarray);

    kin.getDeltaSSGibbs(Rarray);
    deltaGss_Table[i] = egr.deltaRxnVecValue(Rarray);
 
    kin.getDeltaSSEnthalpy(Rarray);
    deltaHss_Table[i] = egr.deltaRxnVecValue(Rarray);

    kin.getDeltaSSEntropy(Rarray);
    deltaSss_Table[i] = egr.deltaRxnVecValue(Rarray);
    
    kin.getDeltaGibbs(Rarray);
    deltaG_Table[i] = egr.deltaRxnVecValue(Rarray);
 
    kin.getDeltaEnthalpy(Rarray);
    deltaH_Table[i] = egr.deltaRxnVecValue(Rarray);

    kin.getDeltaEntropy(Rarray);
    deltaS_Table[i] = egr.deltaRxnVecValue(Rarray);

    Tdelta = T + deltaT;
 
    electrode->getMoleFractions(electrodeMF);
    for (iphKin = 0; iphKin < nPhases; iphKin++) {
      ThermoPhase &gRef = kin.thermo(iphKin);
      int kstart =  electrode->getGlobalSpeciesIndex(iphKin,0);
      gRef.setState_TPX(T, electrode->pressure(), &electrodeMF[kstart]);
    }

    if (iKin) {
      iKin->advanceCoverages(100.0);
    }
    kin.getFwdRatesOfProgress(Rarray);
    kin.getRevRatesOfProgress(R2array);
 
    double FwdROP_delta = egr.FwdROPValue(Rarray, R2array);
    double RevROP_delta = egr.RevROPValue(Rarray, R2array);

    kin.getNetRatesOfProgress(Rarray);
    double NetROP_delta = egr.ROPValue(Rarray);

    calcArrhParam(FwdROP_Table[i], T, FwdROP_delta, Tdelta, Afac, EadivR);
    Afwd_Table[i] = Afac;
    EafwddivR_Table[i] = EadivR;
   
    calcArrhParam(RevROP_Table[i], T, RevROP_delta,  Tdelta, Afac, EadivR);
    Arev_Table[i] = Afac;
    EarevdivR_Table[i] = EadivR;

    if (NetROP_Table[i] >= 0.0) {
      calcArrhParam(NetROP_Table[i], T, NetROP_delta, Tdelta, Afac, EadivR);
      Anet_Table[i] = Afac;
    } else {
      calcArrhParam(-NetROP_Table[i], T, -NetROP_delta, Tdelta, Afac, EadivR);
      Anet_Table[i] = -Afac;
    }
    EanetdivR_Table[i] = EadivR;


  }
    

  setAllBathSpeciesConditions(electrode); 



  delete [] electrodeMF;
  delete [] Sarray;
  delete [] Rarray;
  delete [] R2array;
}
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/**
 *  This routine will print out a table of information a single
 *  reaction, j.
 *
 *
 */
void printGERKineticsTable(Electrode *electrode, int iGER,
			   TemperatureTable& TT,
			   Cantera::Kinetics &kin,
			   ExtraGlobalRxn &egr,
			   Cantera::RxnMolChange *rmc,
			   RxnTempTableStuff &rts) {
  
  // vector<double> & kfwd_Table      = rts.kfwd_Table;
  // vector<double> & krev_Table      = rts.krev_Table;
  // vector<double> & deltaG_Table    = rts.deltaG_Table;
  // vector<double> & deltaH_Table    = rts.deltaH_Table;
  // vector<double> & deltaS_Table    = rts.deltaS_Table;
  vector<double> & deltaGss_Table  = rts.deltaGss_Table;
  vector<double> & deltaHss_Table  = rts.deltaHss_Table;
  vector<double> & deltaSss_Table  = rts.deltaSss_Table;
  vector<double> & Afwd_Table      = rts.Afwd_Table;
  vector<double> & EafwddivR_Table = rts.EafwddivR_Table;
  vector<double> & Arev_Table      = rts.Arev_Table;
  vector<double> & EarevdivR_Table = rts.EarevdivR_Table;
  // vector<double> & kfwdPrime_Table = rts.kfwdPrime_Table;
  // vector<double> & krevPrime_Table = rts.krevPrime_Table;
  vector<double> & NetROP_Table    = rts.NetROP_Table;
  vector<double> & FwdROP_Table    = rts.FwdROP_Table;
  vector<double> & RevROP_Table    = rts.RevROP_Table;

  vector<double> & Anet_Table      = rts.Anet_Table;
  vector<double> & EanetdivR_Table = rts.EanetdivR_Table;
 
  /*
   *  Get the species data object from the Mixture object
   *  this is defined in the constituents.h file, and is 
   *  inherited by Mixture through BaseMix
   */

  int tableWidth = 138;

  /*
   * Conversion factor between Joules/kmol, cantera's units,
   * and kcal / gmol.
   */
  double conv = 1.0 / (4.184E6);
  double convS = conv * 1.0E3;
  double Rkcal = 1.98721E-3;

  ThermoPhase *gThermo = &kin.thermo(0);
  /*
   *  Dump out all of the information about the reaction
   */
  cout << endl;
  print_char('=', tableWidth);
  cout << endl;
  cout << endl;
  cout << "INFORMATION TABLE FOR Global Extra Reaction Reaction \""  << iGER;
  cout << "\" IN PHASE \"";
  cout << gThermo->id();
  cout << "\"" << endl;
  cout << endl;

  /*
   * Reaction String
   */
  string rs = egr.reactionString();
  dnt(1); 
  cout << rs << endl;
  cout << endl;

  /*
   * Print out reactants
   */
  const vector_int& reactants = egr.reactants();
  int nReac = reactants.size();
  dnt(1); 
  if (nReac == 1) {
    cout << "There is one reactant: ";
  } else {
    cout << "There are " << nReac << " reactants: ";
  }
  for (int k = 0; k < nReac; k++) {
    int kindex = reactants[k];
    cout << kin.kineticsSpeciesName(kindex) << " (" << egr.m_ReactantStoich[k] << ") ";
  }
  cout << endl;


  /*
   * Print out products
   */
  const vector_int& products = egr.products();
  int nProd = products.size();
  dnt(1); 
  if (nProd == 1) {
    cout << "There is one product: ";
  } else {
    cout << "There are " << nProd << " products: ";
  }
  for (int k = 0; k < nProd; k++) {
    int kindex = products[k];
    cout << kin.kineticsSpeciesName(kindex) << " (" << egr.m_ProductStoich[k] << ") ";
  }
  cout << endl;

  /*
   *  reaction type
   */
  // double pCurrent = gThermo->pressure();

  dnt(1); 
  cout << "Reaction Type = ";
  cout << "Global Extra Rxn - not in mechanism" << endl;

  /*
   * Reaction Reversible?
   */
  bool rReversible = egr.isReversible();
  dnt(1); 
  if (rReversible) {
    cout << "Reaction is reversible" << endl;
  } else {
    cout << "Reaction is not reversible" << endl;
  }
  cout << endl;

  /*
   * Figure out the units for the forward reaction rate constant
   * First start out with the units for the ROP
   */
  int nphase = kin.nPhases();
  int ndim = 3;
  for (int iph = 0; iph < nphase; iph++) {
    ThermoPhase &tpRef = kin.thermo(iph);
    int idim = tpRef.nDim();
    if (idim < ndim) ndim = idim;
  }
  doublereal unitsROP[6] = { 1.0, (double)(-ndim), 0.0, 0.0, 0.0, -1.0 };
  doublereal unitskfwd[6];
  doublereal unitskrev[6] ;
  doublereal unitsSpecies[6];
  for (int i = 0; i < 6; i++) {
    unitskfwd[i] = unitsROP[i];
    unitskrev[i] = unitsROP[i];
  }
  for (int k = 0; k < nReac; k++) {
    int kindex = reactants[k];
    int iphase = kin.speciesPhaseIndex(kindex);
    ThermoPhase &tpRef = kin.thermo(iphase);
    tpRef.getUnitsStandardConc(unitsSpecies, kindex);
    for (int i = 0; i < 6; i++) {
      unitskfwd[i] -= unitsSpecies[i];
    }
  }
  for (int k = 0; k < nProd; k++) {
    int kindex = products[k];
    int iphase = kin.speciesPhaseIndex(kindex);
    ThermoPhase &tpRef = kin.thermo(iphase);
    tpRef.getUnitsStandardConc(unitsSpecies, kindex);
    for (int i = 0; i < 6; i++) {
      unitskrev[i] -= unitsSpecies[i];
    }
  }
 
  /****************** phases change table ****************************/
  int pwidth = 52;
  if (rmc->m_ChargeTransferInRxn != 0.0) {
    dnt(1);
    printf("This reaction is a charge transfer reaction, n = %g\n",
	   rmc->m_ChargeTransferInRxn);
    pwidth = 86;
  }
 
  dnt(2); print_char('-', pwidth); printf("\n");
 
  dnt(2); printf("|              Phase | Dim Delta_Moles  Delta_Mass |");
  if (rmc->m_ChargeTransferInRxn != 0.0) {
    printf("Electron_Transfer Base_Potential |");
  }
  printf("\n");
  for (int iph = 0; iph < rmc->m_nPhases; iph++) {
    ThermoPhase &tpRef = kin.thermo(iph);
    string sname = tpRef.name();
    dnt(2);
    printf("|%20s|", sname.c_str()); 
    printf("%4d",    rmc->m_phaseDims[iph]);
    printf("%12.5g", rmc->m_phaseMoleChange[iph]);
    printf("%12.5g", rmc->m_phaseMassChange[iph]);
    if (rmc->m_ChargeTransferInRxn != 0.0) {
      printf(" |");
      printf("  %12.5g", rmc->m_phaseChargeChange[iph]);
      double phi = tpRef.electricPotential();
      printf("      %12.5g", phi);
    }
    printf(" |\n");
  }
    
  dnt(2);  print_char('-', pwidth); printf("\n\n");
 

  /*******************************************************************/


  string units_kfwd = formUnitsString(unitskfwd);
  string units_krev = formUnitsString(unitskrev);
  dnt(1);
  string units_ROP =  "kmol / m^" + int2str(ndim) + " s";
  cout << "units_ROP = Units for Global Reaction Rate of progress = " 
       << units_ROP << endl;
 

  /*
   * Dump out the kinetics table as a function of the
   * temperature
   */
  cout << endl;
  cout << "|";
  print_char('-', tableWidth); cout << "|" << endl;
  cout << "|        |                              "
       << " |                                 ";
  if (rReversible) {
    cout << "|                               ";
  } else {
    cout << "| -- REVERSE RXN NOT IN MECH -- ";
  }
  cout << "|                               ";
  cout << "|" << endl;
  cout << "|  Temp  |  FwdROP     Afwd     Ea_fwd  " 
       << " |  DeltaG0    DeltaH0    DeltaS0  "
       << "|  RevROP    A_rev     Ea_rev   "
       << "|  NetROP       Anet     Ea_net ";
  cout << "|" << endl;

  cout << "|   (K)  |    (units_ROP)     (kcal/gmol)"
       << "|(kcal/gmol)(kcal/gmol)(cal/gmolK)"
       << "|  (units_ROP)      (kcal/gmol) |"
       << "   (units_ROP)     (kcal/gmol) |";
  cout << endl;
  cout << "|--------|----------"
       << "---------------------|-----------------"
       << "----------------|-----------------------------"
       << "--|-------------------------------|" << endl;
  
  for (int i = 0; i < TT.size(); i++) {
    cout << "|";
    pr_df(TT[i], 7, 1);
    cout << " |";
    pr_de(FwdROP_Table[i], 10, 2);
    pr_de(Afwd_Table[i], 10, 2);
    pr_dg(EafwddivR_Table[i] * Rkcal, 10, 2);
    cout << " | ";
    pr_dg(deltaGss_Table[i] * conv, 10, 2);
    pr_dg(deltaHss_Table[i] * conv, 10, 2);
    pr_dg(deltaSss_Table[i] * convS, 10, 2);
    cout << "  |";
    pr_de(RevROP_Table[i], 10, 2);
    pr_de(Arev_Table[i], 10, 2);
    pr_dg(EarevdivR_Table[i] * Rkcal, 10, 2);
    cout << " |";
    pr_de(NetROP_Table[i], 10, 2);
    pr_de(Anet_Table[i], 10, 2);
    pr_dg(EanetdivR_Table[i] * Rkcal, 10, 2);
    cout << " |" << endl;
  }
  cout << "|" ;
  print_char('-', tableWidth);
  cout << "|" << endl;
}

// Calculate current at a single voltage
/*
 *
 */
double processGERCurrent(Cantera::RxnMolChange *rmc,
			 Electrode *electrode, int iGERrxn,
			 Kinetics &kin,
			 ExtraGlobalRxn &egr)
{
  int k;
  //InterfaceKinetics *iK = electrode->m_rSurDomain;
  ReactingSurDomain *iK = electrode->currOuterReactingSurface();
  int nr = iK->nReactions();
  std::vector<double> deltaG;
  deltaG.resize(nr, 0.0);
  std::vector<double> deltaSSG;
  deltaSSG.resize(nr, 0.0);
  int iphMetal = electrode->metalPhaseIndex();
  int RSphMetal = iK->PLtoKinPhaseIndex_[iphMetal];
  //int RSphMetal = electrode->CurrReactingSurfacePhaseIndex(iphMetal);
  double phi0Metal = electrode->phaseVoltage(iphMetal);
  double nStoichElectrons = - rmc->m_phaseChargeChange[RSphMetal];

  int iphSoln = electrode->solnPhaseIndex();
  int jphSoln =  iK->PLtoKinPhaseIndex_[iphSoln];
  //int jphSoln = electrode->CurrReactingSurfacePhaseIndex(iphSoln);
  ThermoPhase& tpJSoln = kin.thermo(jphSoln);
  double phi0Soln = tpJSoln.electricPotential();
  ThermoPhase& tpRSMetal = kin.thermo(RSphMetal);
  double phiRSMetal = tpRSMetal.electricPotential();

  double V0 = phiRSMetal - phi0Soln;
  ThermoPhase *tpMetal = &(iK->thermo(RSphMetal));
  ThermoPhase *tpSoln =  &(iK->thermo(jphSoln));
  int iSurf = iK->reactionPhaseIndex();
  ThermoPhase *tps = &(iK->thermo(iSurf));
  SurfPhase *tpSurface = dynamic_cast<SurfPhase *>(tps);

  int nSurfSpecies = tpSurface->nSpecies();
  vector<double> thetaSurf(nSurfSpecies, 0.0);
  tpSurface->getCoverages(DATA_PTR(thetaSurf));

  int pwidth = 100;
  int twidth = 93;
  dnt(1); print_char('-', pwidth); printf("\n");
  string rs = egr.reactionString();
  dnt(1); printf("Voltage Dependence for Rxn: %s", rs.c_str());
  printf("\n\n");
  dnt(2); printf("ThermoPhase index pertaining to the solution = %d\n", jphSoln);
  dnt(2); printf("ThermoPhase index pertaining to the metal = %d\n", RSphMetal); 
  dnt(2); printf("Kinetic Species Index pertaining to the electrons = %d\n", 
		 electrode->kKinSpecElectron(0));
  dnt(2); printf("Input Electric Potential of the Metal = %g Volts\n", phi0Metal);
  dnt(2); printf("Input Electric Potential of the Soln  = %g Volts\n", phi0Soln);
  dnt(2); printf("Input Voltage                         = %g Volts\n", V0);



  // dnt(2); printf("Number of reactant Electrons          = %g \n", nElectrons);
  /*
   * the 1.0E-4 is to change from m-2 to cm-2.
   */
  iK->getDeltaGibbs(DATA_PTR(deltaG));
  double deltaGrxn = egr.deltaRxnVecValue(DATA_PTR(deltaG));

  iK->getDeltaSSGibbs(DATA_PTR(deltaSSG));
  double deltaSSGrxn = egr.deltaRxnVecValue(DATA_PTR(deltaSSG));

 
  double Erxn   = deltaGrxn/Faraday/nStoichElectrons;
  double ESSrxn = deltaSSGrxn/Faraday/nStoichElectrons;
  dnt(2); printf("Delta G (rxn)   = %13.5g J/kmol   --->   ", deltaGrxn);
  dnt(2); printf("E    = %g Volts\n", Erxn);
  dnt(2); printf("Delta G_SS (rxn)= %13.5g J/kmol   --->   ", deltaSSGrxn);
  dnt(2); printf("E_SS = %g Volts\n", ESSrxn);

  dnt(2); printf("Bath Gas conditions:\n");
  dnt(3); printf("T = %g K\n", tpSoln->temperature());
  dnt(3); printf("P = %g Pascal\n", tpSoln->pressure());
  dnt(3); printf("     SurfSpecName         Theta_k_Bath\n");
  for (k = 0; k < nSurfSpecies; k++) {
    string sName = tpSurface->speciesName(k);
    dnt(3); printf("%16s %12.5g\n", sName.c_str(), thetaSurf[k]);
  }
  printf("\n");


  /*
   *  n here is defined as the stoichiometric coefficient for the electron in the global
   *  reaction. Therefore n = -1 means the reaction is written in the cathodic direction,
   *  where electrons are consumed.
   *  n = 1 means that the reaction is written in the anodic direction, where electrons
   *  are created.
   */
  double nF = - Faraday * rmc->m_phaseChargeChange[RSphMetal] * 1.0E-4;

  //bool iEdone = false;
 
  double Voltage, phiMetal;
 
#ifdef DEBUG_HKM_ETABLE
  FILE *EFP = fopen("polarizationData.csv", "w");
#endif
  dnt(2); printf(" TABLE that uses the Surface Solver to Solve for Surface Site Concentrations\n");
 

  int nSS = MIN(nSurfSpecies, 6);
  twidth = 104 + nSS * 11;
  dnt(2); print_char('-', twidth); printf("\n");
  dnt(2);
  cout << "|  Voltage  " 
       << "|       Ifwd   Iexchange_f   Beta_f   "
       << "|       Irev   Iexchange_r   Beta_r   "
       << "|    Inet     "
       << "| ";

  for (k = 0; k < nSS; k++) {
    string sname =  tpSurface->speciesName(k);
    printf("%10s ", sname.c_str());
  }
  cout <<"|" << endl;

#ifdef DEBUG_HKM_ETABLE
  fprintf(EFP,"  Voltage,    Ifwd,   Iexchange_f,   Beta_f,   "
	 "   Irev,   Iexchange_r,   Beta_r,      Inet,  ");
  for (k = 0; k < nSS; k++) {
    string sname =  tpSurface->speciesName(k);
    fprintf(EFP, "%10s ", sname.c_str());
    if (k < nSS - 1) {
      fprintf(EFP,", ");
    }
  }
  fprintf(EFP,"\n");
#endif

  dnt(2);
  cout << "|   Volts   "
       << "|         (amps/cm2)                  "
       << "|         (amps/cm2)                  "
       << "|             " 
       << "| ";
  for (k = 0; k < nSS; k++) {
    printf("           ");
  }
  cout <<"|" << endl;
  dnt(2); print_char('-', twidth); printf("\n");


  std::vector<double> Rfwd;
  Rfwd.resize(nr, 0.0);
  std::vector<double> Rrev;
  Rrev.resize(nr, 0.0);
  std::vector<double> Rnet;
  Rnet.resize(nr, 0.0);


  Voltage = electrode->voltage();
  phiMetal = electrode->phaseVoltage(iphMetal);



  iK->advanceCoverages(1000.0);

  iK->getFwdRatesOfProgress(DATA_PTR(Rfwd));
  iK->getRevRatesOfProgress(DATA_PTR(Rrev));
  double FwdROP = egr.FwdROPValue(DATA_PTR(Rfwd), DATA_PTR(Rrev));
  double RevROP = egr.RevROPValue(DATA_PTR(Rfwd), DATA_PTR(Rrev));

  iK->getNetRatesOfProgress(DATA_PTR(Rnet));
  double NetROP = egr.ROPValue(DATA_PTR(Rnet));

  double deltaPhi = 1.0E-4;
  double phiMetalDelta = phiMetal + deltaPhi;
  double  VoltageDelta = Voltage + deltaPhi;
  tpMetal->setElectricPotential(phiMetalDelta);
  iK->advanceCoverages(1000.0);
  iK->getFwdRatesOfProgress(DATA_PTR(Rfwd));
  iK->getRevRatesOfProgress(DATA_PTR(Rrev));
   
  double FwdROP_delta = egr.FwdROPValue(DATA_PTR(Rfwd), DATA_PTR(Rrev));
  double RevROP_delta = egr.RevROPValue(DATA_PTR(Rfwd), DATA_PTR(Rrev));
  double i0_f;
  double alpha_f;
  double i0_r;
  double alpha_r;
  double FdRT = Faraday/(GasConstant * tpSoln->temperature());

  if (nF > 0.0) {
    calcBF(nF*FwdROP, Voltage, nF*FwdROP_delta, VoltageDelta,  Erxn, 1, FdRT,
	   i0_f, alpha_f);
  } else {
      calcBF(nF*FwdROP, Voltage, nF*FwdROP_delta, VoltageDelta,  Erxn, 0, FdRT,
             i0_f, alpha_f);
  }
  if (nF > 0.0) {
    calcBF(nF*RevROP, Voltage, nF*RevROP_delta, VoltageDelta,  Erxn, 0, FdRT,
	   i0_r, alpha_r);
  } else {
    calcBF(nF*RevROP, Voltage, nF*RevROP_delta, VoltageDelta,  Erxn, 1, FdRT,
	   i0_r, alpha_r);
  }
  

  //  calcBF(-nF*FwdROP, Voltage, -nF*FwdROP_delta, VoltageDelta,  Erxn, 0, FdRT,
  //	 i0_f, alpha_f);
  //calcBF(nF*RevROP, Voltage, nF*RevROP_delta, VoltageDelta,  Erxn, 1, FdRT,
  //	 i0_r, alpha_r);

  dnt(2);
  printf("!%10.6g | %11.5g %11.5g %11.5g | %11.5g %11.5g %11.5g | %11.5g |",
	 Voltage, nF*FwdROP, i0_f, alpha_f, -nF*RevROP, i0_r, alpha_r, nF*NetROP );

  // printf("!%10.4g | %11.5g %11.5g %11.5g | %11.5g %11.5g %11.5g | %11.5g |", 
  //	 Voltage, -nF*FwdROP,i0_f, alpha_f,  nF*RevROP, i0_r, alpha_r, -nF*NetROP );

  double Icurr =  nF*NetROP;

#ifdef DEBUG_HKM_ETABLE
  fprintf(EFP,"!%10.6g | %11.5g %11.5g %11.5g | %11.5g %11.5g %11.5g | %11.5g |",
	Voltage, nF*FwdROP, i0_f, alpha_f, -nF*RevROP, i0_r, alpha_r, nF*NetROP );
  // fprintf(EFP, "%10.4g , %11.5g, %11.5g, %11.5g, %11.5g, %11.5g, %11.5g, %11.5g,", 
  //	  Voltage, -nF*FwdROP,i0_f, alpha_f,  nF*RevROP, i0_r, alpha_r, -nF*NetROP );
#endif

  tpSurface->getCoverages(DATA_PTR(thetaSurf));
  for (k = 0; k < nSS; k++) {
    printf("%10.5g ", thetaSurf[k]);
  }
  printf(" |\n");

#ifdef DEBUG_HKM_ETABLE
  for (k = 0; k < nSS; k++) {
    fprintf(EFP,"%10.5g ", thetaSurf[k]);
    if (k < nSS - 1) {
      fprintf(EFP,", ");
    }
  }
  fprintf(EFP,"\n");
#endif
  //  }
  dnt(2); print_char('-', twidth); printf("\n");
  
  tpMetal->setElectricPotential(phiMetal);
  iK->advanceCoverages(1000.0);

  dnt(1); print_char('-', pwidth); printf("\n");
#ifdef DEBUG_HKM_ETABLE
  fclose(EFP);
#endif
  return Icurr;
}

/*************************************************************************/
void processGERCurrentVsPotTable(Cantera::RxnMolChange *rmc,
				 Electrode *electrode, int iGERrxn,
				 TemperatureTable &TT,
				 Kinetics &kin,  
				 ExtraGlobalRxn &egr,
				 RxnTempTableStuff &rts)
{
  int iph, k, e;
  ThermoPhase *tp = 0;
  ThermoPhase *tpMetal = 0;
  int nTotalSpecies = kin.nTotalSpecies();
  int nRxns         = kin.nReactions();
  vector<double> Rfwd(nTotalSpecies, 0.0);
  vector<double> Rrev(nTotalSpecies, 0.0);
  vector<double> Rnet(nTotalSpecies, 0.0);
  vector<double> deltaG(nRxns, 0.0);
  vector<double> deltaSSG(nRxns, 0.0);

  InterfaceKinetics *iK = dynamic_cast<InterfaceKinetics *>(&kin);
  if (!iK) {
    throw CanteraError("RxnMolChange", "unknown condition on charge");
  }
  /*
   * First thing is to figure out how to calculate the voltage
   *    V = phi_metal - phi_soln
   * Look for the electron species. The phase where this occurs will be
   * called the metal, iMetal;
   */
  int iMetal = 0;
  int iSoln = 0;
  int kElectron = -1;
  int nPhases = iK->nPhases();
  for (iph = 0; iph < nPhases; iph++) {
    tp = &(iK->thermo(iph));
    int nSpecies = tp->nSpecies();
    int nElements = tp->nElements();
    int eElectron = tp->elementIndex("E");
    if (eElectron >= 0) {
      for (k = 0; k < nSpecies; k++) {
	if (tp->nAtoms(k,eElectron) == 1) {
	  int ifound = 1;
	  for (e = 0; e < nElements; e++) {
	    if (tp->nAtoms(k,e) != 0.0) {
	      if (e != eElectron) {
		ifound = 0;
	      }
	    }
	  }
	  if (ifound == 1) {
	    iMetal = iph;
	    kElectron = iK->kineticsSpeciesIndex(k, iph);
	  }
	}
      }
    }
    if (iph != iMetal) {
      if (rmc->m_phaseChargeChange[iph] != 0) {
	iSoln = iph;
      }
    }
  }
 
 // int iElectrodeMetal = electrode->metalPhaseIndex();
  //int iElectrodeLyte  = electrode->solnPhaseIndex();

  //double phi0Metal = electrode->phaseVoltage(iElectrodeMetal);

  ThermoPhase &tpRefMetal = kin.thermo(iMetal);
  double phi0Metal = tpRefMetal.electricPotential();
  ThermoPhase &tpRefSoln = kin.thermo(iSoln);
  double phi0Soln = tpRefSoln.electricPotential();


  double V0 = phi0Metal = phi0Soln;
  tpMetal = &(iK->thermo(iMetal));
  ThermoPhase *tpSoln =  &(iK->thermo(iSoln));
  int iSurf = iK->reactionPhaseIndex();
  ThermoPhase *tps = &(iK->thermo(iSurf));
  SurfPhase *tpSurface = dynamic_cast<SurfPhase *>(tps);

  int nSurfSpecies = tpSurface->nSpecies();
  vector<double> thetaSurf(nSurfSpecies, 0.0);
  tpSurface->getCoverages(DATA_PTR(thetaSurf));

  int pwidth = 100;
  int twidth = 93;
  dnt(1); print_char('-', pwidth); printf("\n");
  string rs = egr.reactionString();
  dnt(1); printf("Voltage Dependence for Rxn: %s", rs.c_str());
  printf("\n\n");
  dnt(2); printf("ThermoPhase index pertaining to the solution = %d\n", iSoln);
  dnt(2); printf("ThermoPhase index pertaining to the metal = %d\n", iMetal); 
  dnt(2); printf("Kinetic Species Index pertaining to the electrons = %d\n", kElectron);
  dnt(2); printf("Input Electric Potential of the Metal = %g Volts\n", phi0Metal);
  dnt(2); printf("Input Electric Potential of the Soln  = %g Volts\n", phi0Soln);
  dnt(2); printf("Input Voltage                         = %g Volts\n", V0);



  // dnt(2); printf("Number of reactant Electrons          = %g \n", nElectrons);
  /*
   * the 1.0E-4 is to change from m-2 to cm-2.
   */
  iK->getDeltaGibbs(DATA_PTR(deltaG));
  double deltaGrxn = egr.deltaRxnVecValue(DATA_PTR(deltaG));

  iK->getDeltaSSGibbs(DATA_PTR(deltaSSG));
  double deltaSSGrxn = egr.deltaRxnVecValue(DATA_PTR(deltaSSG));

 
  double Erxn = - deltaGrxn/Faraday/rmc->m_ChargeTransferInRxn;
  double ESSrxn = - deltaSSGrxn/Faraday/rmc->m_ChargeTransferInRxn;
  dnt(2); printf("Delta G (rxn)   = %13.5g J/kmol   --->   ", deltaGrxn);
  dnt(2); printf("E    = %g Volts\n", Erxn);
  dnt(2); printf("Delta G_SS (rxn)= %13.5g J/kmol   --->   ", deltaSSGrxn);
  dnt(2); printf("E_SS = %g Volts\n", ESSrxn);

  dnt(2); printf("Bath Gas conditions:\n");
  dnt(3); printf("T = %g K\n", tpSoln->temperature());
  dnt(3); printf("P = %g Pascal\n", tpSoln->pressure());
  dnt(3); printf("     SurfSpecName         Theta_k_Bath\n");
  for (k = 0; k < nSurfSpecies; k++) {
    string sName = tpSurface->speciesName(k);
    dnt(3); printf("%16s %12.5g\n", sName.c_str(), thetaSurf[k]);
  }
  printf("\n");



  double nF = Faraday * rmc->m_phaseChargeChange[iMetal] * 1.0E-4; 
  bool iEdone = false;
 
  double Voltage, phiMetal;
 
#ifdef DEBUG_HKM_ETABLE
  FILE *EFP = fopen("polarizationData.csv", "w");
#endif
  dnt(2); printf(" TABLE that uses the Surface Solver to Solve for Surface Site Concentrations\n");
 

  int nSS = MIN(nSurfSpecies, 6);
  twidth = 104 + nSS * 11;
  dnt(2); print_char('-', twidth); printf("\n");
  dnt(2);
  cout << "|  Voltage  " 
       << "|       Ifwd   Iexchange_f   Beta_f   "
       << "|       Irev   Iexchange_r   Beta_r   "
       << "|    Inet     "
       << "| ";

  for (k = 0; k < nSS; k++) {
    string sname =  tpSurface->speciesName(k);
    printf("%10s ", sname.c_str());
  }
  cout <<"|" << endl;

#ifdef DEBUG_HKM_ETABLE
  fprintf(EFP,"  Voltage,    Ifwd,   Iexchange_f,   Beta_f,   "
	  "   Irev,   Iexchange_r,   Beta_r,      Inet,  ");
  for (k = 0; k < nSS; k++) {
    string sname =  tpSurface->speciesName(k);
    fprintf(EFP, "%10s ", sname.c_str());
    if (k < nSS - 1) {
      fprintf(EFP,", ");
    }
  }
  fprintf(EFP,"\n");
#endif

  dnt(2);
  cout << "|   Volts   "
       << "|         (amps/cm2)                  "
       << "|         (amps/cm2)                  "
       << "|             " 
       << "| ";
  for (k = 0; k < nSS; k++) {
    printf("           ");
  }
  cout <<"|" << endl;
  dnt(2); print_char('-', twidth); printf("\n");
  iEdone = false;

  double Voltage0 = V0 - 0.60;
  double phiMetal0 = phi0Metal - 0.60;
  for (int iV = 0; iV < 120; iV++) {
    Voltage0 += 0.02;
    phiMetal0 += 0.02;
    
    if (!iEdone && (Voltage0 > Erxn)) {
      iEdone = true;
      Voltage = Erxn;
      ThermoPhase &tpRefSoln = kin.thermo(iSoln);
      double phiTmp = tpRefSoln.electricPotential();
      phiMetal = Erxn + phiTmp;
    } else {
      Voltage = Voltage0;
      phiMetal = phiMetal0;
    }
    if (fabs(Voltage) < 1.0E-10) Voltage = 0.0;
    tpMetal->setElectricPotential(phiMetal);
    iK->advanceCoverages(1000.0);

    iK->getFwdRatesOfProgress(DATA_PTR(Rfwd));
    iK->getRevRatesOfProgress(DATA_PTR(Rrev));
    double FwdROP = egr.FwdROPValue(DATA_PTR(Rfwd), DATA_PTR(Rrev));
    double RevROP = egr.RevROPValue(DATA_PTR(Rfwd), DATA_PTR(Rrev));

    iK->getNetRatesOfProgress(DATA_PTR(Rnet));
    double NetROP = egr.ROPValue(DATA_PTR(Rnet));

    double deltaPhi = 1.0E-4;
    double phiMetalDelta = phiMetal + deltaPhi;
    double  VoltageDelta = Voltage + deltaPhi;
    tpMetal->setElectricPotential(phiMetalDelta);
    iK->advanceCoverages(1000.0);
    iK->getFwdRatesOfProgress(DATA_PTR(Rfwd));
    iK->getRevRatesOfProgress(DATA_PTR(Rrev));
   
    double FwdROP_delta = egr.FwdROPValue(DATA_PTR(Rfwd), DATA_PTR(Rrev));
    double RevROP_delta = egr.RevROPValue(DATA_PTR(Rfwd), DATA_PTR(Rrev));
    double i0_f;
    double alpha_f;
    double i0_r;
    double alpha_r;
    double FdRT = Faraday/(GasConstant * tpSoln->temperature());
    calcBF(-nF*FwdROP, Voltage, -nF*FwdROP_delta, VoltageDelta,  Erxn, 0, FdRT,
	   i0_f, alpha_f);
    calcBF(nF*RevROP, Voltage, nF*RevROP_delta, VoltageDelta,  Erxn, 1, FdRT,
	   i0_r, alpha_r);

    dnt(2);
    printf("!%10.4g | %11.5g %11.5g %11.5g | %11.5g %11.5g %11.5g | %11.5g |", 
	   Voltage, -nF*FwdROP,i0_f, alpha_f,  nF*RevROP, i0_r, alpha_r, -nF*NetROP );

#ifdef DEBUG_HKM_ETABLE
    fprintf(EFP, "%10.4g , %11.5g, %11.5g, %11.5g, %11.5g, %11.5g, %11.5g, %11.5g,", 
	    Voltage, -nF*FwdROP,i0_f, alpha_f,  nF*RevROP, i0_r, alpha_r, -nF*NetROP );
#endif

    tpSurface->getCoverages(DATA_PTR(thetaSurf));
    for (k = 0; k < nSS; k++) {
      printf("%10.5g ", thetaSurf[k]);
    }
    printf(" |\n");

#ifdef DEBUG_HKM_ETABLE
    for (k = 0; k < nSS; k++) {
      fprintf(EFP,"%10.5g ", thetaSurf[k]);
      if (k < nSS - 1) {
	fprintf(EFP,", ");
      }
    }
    fprintf(EFP,"\n");
#endif
  }
  dnt(2); print_char('-', twidth); printf("\n");
  


  dnt(1); print_char('-', pwidth); printf("\n");
#ifdef DEBUG_HKM_ETABLE
  fclose(EFP);
#endif
}

