/*
 * $Id: epequil_prep.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "epequil_prep.h"
#include "epequil_input.h"

#include <string>
using std::string;
using namespace Cantera;
                              
/*****************************************************************************
 *
 */
int epequil_prep(EPEQUIL_INPUT *prob) 
{
    MultiPhase *mp = prob->m_mp;
    mp->setTemperature(prob->T);
    mp->setPressure(prob->Pres);
    if (prob->specifiedElementAbundances) {
      printf("specified Element Abundances not handled yet\n");
      exit(-1);
    } else {
      mp->setMoles(prob->spMoles);
    }     
    return 0;
}


/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
   
   
   
