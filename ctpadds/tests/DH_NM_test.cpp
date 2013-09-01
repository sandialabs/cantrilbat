/* ======================================================================= */
/* $RCSfile: DH_NM_test.cpp,v $ */
/* $Author: hkmoffa $ */
/* $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $ */
/* $Revision: 508 $ */
/* ======================================================================= */
/*
 * Small test problem for the Newman formulation of the DH theory
 */
#include <stdio.h>

#include "../../util_src/src/mdp_allo.h"

#include "cantera/thermo/DebyeHuckel.h"
#include "cantera/base/Array.h"

using namespace std;
using namespace Cantera;
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
  
     DebyeHuckel *DH = new DebyeHuckel("DH_NaCl_NM.xml", "NaCl_electrolyte");
     
     int nsp = DH->nSpecies();
     const vector<string>& sn = DH->speciesNames();
    
     /*
      *
      */
     double a1 = DH->AionicRadius(1);
     printf("a1 = %g\n", a1);
     double a2 = DH->AionicRadius(2);
     printf("a2 = %g\n", a2);

     Array2D& beta = DH->get_Beta_ij();
     for (int i = 0; i < nsp; i++) {
       for (int j = 0; j < nsp; j++) {
	 printf("Beta(%s, %s) = %g\n", sn[i].c_str(), 
		sn[j].c_str(), beta(i,j));
       }
     }

     double *acMol = mdp_alloc_dbl_1(100, 1.0);

     pAtable(DH);
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

     mdp_safe_free((void **) &acMol);
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
