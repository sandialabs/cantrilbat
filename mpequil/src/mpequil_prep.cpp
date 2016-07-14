/*
 * $Id: mpequil_prep.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */


#include "mpequil_prep.h"
#include "mpequil_input.h"

#include <string>
using std::string;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif


                              
/*****************************************************************************
 *
 */
int mpequil_prep(MPEQUIL_INPUT *prob) 
{
    MultiPhase *mp = prob->m_mp;
    mp->setTemperature(prob->T);
    mp->setPressure(prob->Pres);
    if (prob->specifiedElementAbundances) {
      printf("specified Element Abundances not handled yet\n");
    }
    mp->setMoles(prob->spMoles);
    return 0;
}

void mpequil_query(MPEQUIL_INPUT *prob) 
{
    int k, p;
    ThermoPhase *tphase;
    MultiPhase *mp = prob->m_mp;
    for (k = 0; k < prob->nspecies; k++) {
      prob->spMoles[k] = mp->speciesMoles(k);
    }
  
    int loc = 0;
    for (p = 0; p < prob->nphase; p++) {
      tphase = &(mp->phase(p));
      tphase->getMoleFractions(prob->spMf + loc);
      loc += tphase->nSpecies();
    }

    loc = 0;
    for (p = 0; p < prob->nphase; p++) {
      tphase = &(mp->phase(p));
      tphase->getChemPotentials(prob->spChemPot + loc);
      loc += tphase->nSpecies();
    }


}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
   
   
   
