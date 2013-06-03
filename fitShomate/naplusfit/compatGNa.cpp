/**
 *
 *  @file HMW_graph_1.cpp
 */

/*
 *  $Author: hkmoffa $
 *  $Date: 2013-01-07 15:32:48 -0700 (Mon, 07 Jan 2013) $
 *  $Revision: 504 $
 */
#include <stdio.h>

#include "vcsc.h"
#include "ThermoPhase.h"


#include "vcs_Cantera_input.h"
#include "vcs_Cantera_convert.h"
#include "cantera/base/logger.h"
#include "cantera/thermo.h"


#include "TemperatureTable.h"
#include "HMWSoln.h"

using namespace Cantera;

class fileLog: public Logger {
public:
    fileLog(string fName) {
	m_fName = fName;
	m_fs.open(fName.c_str());
    }
    
    virtual void write(const string& msg) {
          m_fs << msg;
	  m_fs.flush();
    }

    virtual ~fileLog() {
	m_fs.close();
    }
    
    string m_fName;
    ofstream m_fs;

};

void printUsage() {
    cout << "usage: HMW_test " <<  endl;
    cout <<"                -> Everything is hardwired" << endl;
}

void pAtable(HMWSoln *HMW) {
    int nsp = HMW->nSpecies();
    double *acMol = mdp_alloc_dbl_1(nsp, 1.0);
    double *mf = mdp_alloc_dbl_1(nsp, 0.0);
    mf[1] = 0.0;
    double *activities =  mdp_alloc_dbl_1(100, 1.0);
    double *moll =  mdp_alloc_dbl_1(100, 0.0);
    
    HMW->getMolalityActivityCoefficients(acMol);
    HMW->getMoleFractions(mf);
    HMW->getActivities(activities);
    HMW->getMolalities(moll);
    string sName;
    printf("            Name      Activity  ActCoeffMolal "
	   "   MoleFract      Molality\n");
    for (int k = 0; k < nsp; k++) {
      sName = HMW->speciesName(k);
      printf("%16s %13g %13g %13g %13g\n", 
	     sName.c_str(), activities[k], acMol[k], mf[k], moll[k]);
    }
    mdp_safe_free((void **) &acMol);
    mdp_safe_free((void **) &mf);
    mdp_safe_free((void **) &activities);
    mdp_safe_free((void **) &moll);
}


double q16 = 41587.11;
double q17 = -315.90;
double q18 = 0.8514;
double q19 = -8.3637E-4;
double q20 = -1729.93;
double q20delta = -0.001669;

double delG(double T) {
    double tt = T * T;
    double tmp = q16 - q17 * T * log(T) - q18 * tt - 0.5 * q19 * tt * T
	+ q20 * T + q20delta * T;
    return tmp;
}

int main(int argc, char **argv)
{

   int retn;
   int i;
 
   try {
     Cantera::ThermoPhase *tp = 0;
     char iFile[80];
     //strcpy(iFile, "HMW_NaCl.xml");
     //strcpy(iFile, "HMW_NaCl_alt.xml");
     strcpy(iFile, "HMW_NaCl_sp1977.xml");
     if (argc > 1) {
       strcpy(iFile, argv[1]);
     }
 
     //fileLog *fl = new fileLog("HMW_graph_1.log");
     //setLogger(fl);

     HMWSoln *HMW = new HMWSoln(iFile, "NaCl_electrolyte");


     /*
      * Load in and initialize the 
      */
     Cantera::ThermoPhase *solid = Cantera::newPhase("NaCl_Solid.xml","NaCl(S)");
 
     
     int nsp = HMW->nSpecies();
     double *acMol = mdp_alloc_dbl_1(100, 1.0);
     double *act = mdp_alloc_dbl_1(100, 1.0);
     //pAtable(HMW);
     double *mf = mdp_alloc_dbl_1(100, 0.0);
     mf[1] = 0.0;
 
     double *moll =  mdp_alloc_dbl_1(100, 0.0);

     HMW->getMoleFractions(mf);
     string sName;

     TemperatureTable TTable(31, true, 293.15, 10., 0, 0);

     HMW->setState_TP(298.15, 1.01325E5);
  
     int i1 = HMW->speciesIndex("Na+");
     int i2 = HMW->speciesIndex("Cl-");
     //int i3 = HMW->speciesIndex("H2O(L)");
     for (i = 1; i < nsp; i++) {
       moll[i] = 0.0;
     }
     HMW->setMolalities(moll);

     double ISQRT;
     double Is = 0.0;

     /*
      * Set the Pressure
      */
     double pres = OneAtm;

     /*
      * Fix the molality
      */
     Is = 6.2;
     ISQRT = sqrt(Is);
     moll[i1] = Is;
     moll[i2] = Is;
     HMW->setState_TPM(298.15, pres, moll);
     double Xmol[30];
     HMW->getMoleFractions(Xmol);

     /*
      * ThermoUnknowns
      */
     double mu0_RT[20], mu[20];
     double mu0_NaCl, mu0_Naplus_base, mu0_Clminus, Delta_G0;
     double mu_NaCl, mu_Naplus, mu_Clminus, Delta_G;
     double Delta_G0_base, mu0_Naplus;
     double molarGibbs0, molarGibbs;

     /*
      * Create a Table of NaCl Enthalpy Properties as a Function
      * of the Temperature
      */

     printf("              T,        Aphi,      Delta_GO_base Delta_G0,"
	    "      mu0_Na+_base,"
	    "   mu0_Na+,   "
	    "   \n");
     printf("        Kelvin, sqrt(kg/gmol),      kJ/gmol       kJ/gmol,"
	    "       kJ/gmol,"
	    "       kJ/gmol"
	    "   \n");
     for (i = 0; i < TTable.NPoints; i++) {
       double T = TTable.T[i];
       double RT = GasConstant * T;
       HMW->setState_TPM(T, pres, moll);
       solid->setState_TP(T, pres);
       /*
	* Get the Standard State DeltaH
	*/
       solid->getGibbs_RT(mu0_RT);
       mu0_NaCl = mu0_RT[0] * RT * 1.0E-6;
       HMW->getGibbs_RT(mu0_RT);
       mu0_Naplus_base = mu0_RT[i1] * RT * 1.0E-6;
       mu0_Clminus = mu0_RT[i2] * RT * 1.0E-6;
       /*
	* Calculate the base delta G0, from the xml files.
	* This is before the transformation.
	*/
       Delta_G0_base =  (mu0_Naplus_base + mu0_Clminus) - mu0_NaCl;

       HMW->getMolalityActivityCoefficients(acMol);
       HMW->getActivities(act);
       double meanAC = sqrt(acMol[i1] * acMol[i2]);
       solid->getChemPotentials(mu);
       mu_NaCl = mu[0] * 1.0E-6;
       HMW->getChemPotentials(mu);
       mu_Naplus  = mu[i1] * 1.0E-6;
       mu_Clminus = mu[i2] * 1.0E-6;
       Delta_G = (mu_Naplus + mu_Clminus) - mu_NaCl;

       /*
	* Calculate the delta G value for the reaction from
	* Silvester and Pitzer's polynomial.
	*   units are in J/gmol.
	*/
       Delta_G0 = 4.184 * delG(T) / 1.0E3;
  
  
       mu0_Naplus = Delta_G0 + mu0_NaCl - mu0_Clminus;

       molarGibbs = HMW->gibbs_mole() * 1.0E-6;
       double Aphi = HMW->A_Debye_TP() / 3.0;
    
       for (int k = 0; k < nsp; k++) {
	 mu0_RT[k] *= RT * 1.0E-6;
       }

       molarGibbs0 = 0.0;
       for (int k = 0; k < nsp; k++) {
	 molarGibbs0 += Xmol[k] * mu0_RT[k];
       }
       double G_ex = molarGibbs - molarGibbs0;

       if (fabs (T-298.15) < 1.0) {
	 printf("mu0_Naplus  = %g\n", mu0_Naplus_base);
	 printf("mu0_Clminus = %g\n", mu0_Clminus);
	 printf("mu0_NaCl(s) = %g,   mu_NaCl(s) = %g\n",mu0_NaCl, mu_NaCl);
       }

       printf("%13g, %13g, %13g, %13g, %13g, %13g\n",
	      T, Aphi, Delta_G0_base, Delta_G0, mu0_Naplus_base, mu0_Naplus
	      );
     }

     mdp_safe_free((void **) &mf);
     mdp_safe_free((void **) &moll);
     delete HMW;
     HMW = 0;
     delete solid;
     solid = 0;
     Cantera::appdelete();

     return retn;

   } catch (CanteraError) {

     showErrors();
     Cantera::appdelete();
     return -1;
   }
} 
