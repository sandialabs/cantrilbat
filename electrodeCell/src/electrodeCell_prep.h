
/*
 * $Id: electrodeCell_prep.h 21 2012-03-01 00:06:01Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */


#ifndef _ELECTRODE_PREP_H
#define _ELECTRODE_PREP_H

#include "cantera/base/config.h"

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif 
{
  class Electrode;
}
#ifdef useZuzaxNamespace
extern int electrode_prep(Zuzax::Electrode *);
extern void electrode_query(Zuzax::Electrode *prob);  
#else
extern int electrode_prep(Cantera::Electrode *);
extern void electrode_query(Cantera::Electrode *prob);  

#endif


#endif
