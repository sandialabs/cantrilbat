/* ======================================================================= */
/* $RCSfile: DH_test.cpp,v $ */
/* $Author: hkmoffa $ */
/* $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $ */
/* $Revision: 508 $ */
/* ======================================================================= */

#include <stdio.h>

#include "mdp_allo.h"

#include "cantera/thermo/DebyeHuckel.h"

using namespace Cantera;
using namespace std;
using namespace mdpUtil;


void printUsage() {
    cout << "usage: DH_test " <<  endl;
    cout <<"                -> Everything is hardwired" << endl;
}

void pAtable(DebyeHuckel *DH) {
    int nsp = DH->nSpecies();
    double *acMol = mdp_alloc_dbl_1(nsp, 1.0);
    double *mf = mdp_alloc_dbl_1(nsp, 0.0);
    mf[1] = 0.0;
    double *activities =  mdp_alloc_dbl_1(100, 1.0);
    double *moll =  mdp_alloc_dbl_1(100, 0.0);
    
    DH->getMolalityActivityCoefficients(acMol);
    DH->getMoleFractions(mf);
    DH->getActivities(activities);
    DH->getMolalities(moll);
    string sName;
    printf("            Name      Activity  ActCoeffMolal "
	   "   MoleFract      Molality\n");
    for (int k = 0; k < nsp; k++) {
      sName = DH->speciesName(k);
      printf("%16s %13g %13g %13g %13g\n", 
	     sName.c_str(), activities[k], acMol[k], mf[k], moll[k]);
    }
    mdp_safe_free((void **) &acMol);
    mdp_safe_free((void **) &mf);
    mdp_safe_free((void **) &activities);
    mdp_safe_free((void **) &moll);
}

int main(int argc, char **argv)
{

  int retn = 0;
  string commandFile = "vcs_Cantera.inp";
  try {
  
    DebyeHuckel *DH = new DebyeHuckel("DH_NaCl.xml", "NaCl_electrolyte");
     
    int nsp = DH->nSpecies();
    
    /*
     *
     */
    double a1 = DH->AionicRadius(1);
    printf("a1 = %g\n", a1);
    double a2 = DH->AionicRadius(2);
    printf("a2 = %g\n", a2);
    double *acMol = mdp_alloc_dbl_1(100, 1.0);
    double *mu0 = mdp_alloc_dbl_1(100, 0.0);
    // pAtable(DH);
    double *mf = mdp_alloc_dbl_1(100, 0.0);
    mf[1] = 0.0;
 
    double *moll =  mdp_alloc_dbl_1(100, 0.0);

    DH->getMolalityActivityCoefficients(acMol);
    DH->getMoleFractions(mf);
    DH->getMolalities(moll);
    string sName;
  
    DH->getMolalities(moll);
    moll[1] = 9.3549;
    moll[2] = 9.3549;
    DH->setMolalities(moll);
    DH->m_useHelgesonFixedForm = true;
    pAtable(DH);

    DH->setState_TP(298.15, 1.01325E5);
    DH->getStandardChemPotentials(mu0);
    // translate from J/kmol to kJ/gmol
    int k;
    for (k = 0; k < nsp; k++) {
      mu0[k] *= 1.0E-6;                                                                                                                  
    }

    printf("           Species   Standard chemical potentials (kJ/gmol) \n");
    printf("------------------------------------------------------------\n");
    for (k = 0; k < nsp; k++) {
      sName = DH->speciesName(k);
      printf("%16s %16.9g\n", sName.c_str(), mu0[k]);
    }
    printf("------------------------------------------------------------\n");
    printf(" Some DeltaSS values:               Delta(mu_0)\n");
    double deltaG;
    int i1, i2, i3, j1, i4;
    double RT = 8.314472E-3 * 298.15;

    i1 = DH->speciesIndex("Na+");
    i2 = DH->speciesIndex("OH-");
    j1 = DH->speciesIndex("NaOH(aq)");
    if (i1 < 0 || i2 < 0 || j1 < 0) {
      printf("problems\n");
      exit(-1);
    }
    deltaG = mu0[j1] - mu0[i1] - mu0[i2];
    printf(" NaOH(aq): Na+ + OH- -> NaOH(aq): %14.7g kJ/gmol \n",
	   deltaG);
    printf("                                : %14.7g (dimensionless) \n",
	   deltaG/RT);
     
    i1 = DH->speciesIndex("Na+");
    i2 = DH->speciesIndex("H2O(L)");
    i3 = DH->speciesIndex("H+");
    j1 = DH->speciesIndex("NaOH(aq)");
    if (i1 < 0 || i2 < 0 || j1 < 0) {
      printf("problems\n");
      exit(-1);
    }
    deltaG = mu0[j1] - mu0[i1] - mu0[i2] + mu0[i3];
    printf(" NaOH(aq): Na+ + H2O(L) - H+ -> NaOH(aq): %14.7g kJ/gmol \n",
	   deltaG);
    printf("                                : %14.7g (dimensionless) \n",
	   deltaG/RT);

    i1 = DH->speciesIndex("Na+");
    i2 = DH->speciesIndex("Cl-");
    j1 = DH->speciesIndex("NaCl(aq)");
    if (i1 < 0 || i2 < 0 || j1 < 0) {
      printf("problems\n");
      exit(-1);
    }
    deltaG = mu0[j1] - mu0[i1] - mu0[i2];
    printf(" NaCl(aq): Na+ + Cl- -> NaCl(aq): %14.7g kJ/gmol \n",
	   deltaG);
    printf("                                : %14.7g (dimensionless) \n",
	   deltaG/RT);
    printf("                                : %14.7g (dimensionless/ln10) \n",
	   deltaG/(RT * log(10.0)));


    i1 = DH->speciesIndex("Na+");
    i2 = DH->speciesIndex("Cl-");
    deltaG = -432.6304 - mu0[i1] - mu0[i2];
    printf(" NaCl(S): Na+ + Cl- -> NaCl(S): %14.7g kJ/gmol \n",
	   deltaG);
    printf("                                : %14.7g (dimensionless) \n",
	   deltaG/RT);
    printf("                                : %14.7g (dimensionless/ln10) \n",
	   deltaG/(RT * log(10.0)));

    i1 = DH->speciesIndex("H+");
    i2 = DH->speciesIndex("H2O(L)");
    j1 = DH->speciesIndex("OH-");
    if (i1 < 0 || i2 < 0 || j1 < 0) {
      printf("problems\n");
      exit(-1);
    }
    deltaG = mu0[j1] + mu0[i1] - mu0[i2];
    printf(" OH-: H2O(L) - H+ -> OH-: %14.7g kJ/gmol \n",
	   deltaG);
    printf("                                : %14.7g (dimensionless) \n",
	   deltaG/RT);
    printf("                                : %14.7g (dimensionless/ln10) \n",
	   deltaG/(RT * log(10.0)));

    i1 = DH->speciesIndex("H+");
    i2 = DH->speciesIndex("H2O(L)");
    i3 = DH->speciesIndex("Na+");
    i4 = DH->speciesIndex("SiO2(aq)");
    j1 = DH->speciesIndex("NaH3SiO4(aq)");
    if (i1 < 0 || i2 < 0 || i3 < 0 || i4 < 0 || j1 < 0) {
      printf("problems\n");
      exit(-1);
    }
    deltaG = mu0[j1] -(- mu0[i1] + 2*mu0[i2] + mu0[i3] + mu0[i4]);
    printf(" NaH3SiO4(aq): -H+ +Na+ +2H2O + SiO2: %14.7g kJ/gmol \n",
	   deltaG);
    printf("                                : %14.7g (dimensionless) \n",
	   deltaG/RT);
    printf("                                : %14.7g (dimensionless/ln10) \n",
	   deltaG/(RT * log(10.0)));

    i1 = DH->speciesIndex("H+");
    i2 = DH->speciesIndex("H2O(L)");
    i3 = DH->speciesIndex("SiO2(aq)");
    j1 = DH->speciesIndex("H3SiO4-");
    if (i1 < 0 || i2 < 0 || i3 < 0  || j1 < 0) {
      printf("problems\n");
      exit(-1);
    }
    deltaG = mu0[j1] -(- mu0[i1] + 2*mu0[i2] + mu0[i3] );
    printf(" H3SiO4-: -H+  +2H2O + SiO2: %14.7g kJ/gmol \n",
	   deltaG);
    printf("                                : %14.7g (dimensionless) \n",
	   deltaG/RT);
    printf("                                : %14.7g (dimensionless/ln10) \n",
	   deltaG/(RT * log(10.0)));
    printf("------------------------------------------------------------\n");

    mdp_safe_free((void **) &mf);

    mdp_safe_free((void **) &moll);
    delete DH;
    DH = 0;
    Cantera::appdelete();

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
  return retn;
} 
