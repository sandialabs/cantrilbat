/**
 *  @file example2.cpp
 *
 * $Id: cttables_thermo.cpp 497 2013-01-07 21:17:04Z hkmoffa $
 * 
 */


//  Example 2
//
//  Read a mechanism, and print to the standard output stream a
//  well-formatted Chemkin ELEMENT section.
//

#include "IdealReactingGas.h"

#include "TemperatureTable.h"

#include "LE_PickList.h"
#include "BE_MoleComp.h"

#include "mdp_allo.h"

#include "cantera/thermo/IdealSolidSolnPhase.h"

#include "cantera/thermo/StoichSubstance.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/VPStandardStateTP.h"
#include "cantera/thermo/PDSS.h"



#include <iostream>
#include <new>
#include <string>
#include <vector>
#include <typeinfo>

using namespace Cantera;
using namespace std;



#include "cttables.h"

/***********************************************************************
 * printVolSpecies():
 *
 *  This routine adds onto a species table the molar volume information
 *  supplied for the species. This is not part of the reference state.
 *  This is part of the description of the phase. 
 * 
 */
void printThermoCoeffSpecies(ThermoPhase *g_ptr, int k) {

  int j, nint, nzones = 1, i;
  string sName = g_ptr->speciesName(k);
  double *cptr;
  int index;
  size_t kindex;
  int type;
  PDSS * pdss_ptr = 0;
  SpeciesThermo& sp = g_ptr->speciesThermo();
  VPStandardStateTP *  vpss = 0;
  PDSS_enumType ptype;

  int rt = sp.reportType(k);
  double c[200];
  int rt2;
  double minTemp, maxTemp, refPressure;
  sp.reportParams(k, rt2, c, minTemp, maxTemp, refPressure);
  if (rt != rt2) {
    printf("report types don't match ERROR\n");
    exit(-1);
  }


  switch (rt) {
  case NASA2:
    dnt(1); printf("NASA Polynomial format: 2 zones\n");
    dnt(2); printf("Low  temperature polynomials: %g < T < %g: \n",
		   minTemp, c[0]);
    //dnt(2); printf("%12.5g %12.5g %12.5g %12.5g\n",
    //		   c[1], c[2],  c[3],  c[4]);
    dnt(2); printf("%12.5g %12.5g %12.5g %12.5g\n",
    		   c[3], c[4],  c[5],  c[6]);
    dnt(2); printf("%12.5g %12.5g %12.5g\n",
		   c[7], c[1], c[2] );
    dnt(2); printf("High temperature polynomials: %g < T < %g: \n",
		   c[0], maxTemp);
    dnt(2); printf("%16.10g %16.10g %16.10g %16.10g\n",
		   c[10], c[11], c[12],  c[13]);
    dnt(2); printf("%16.10g %16.10g %16.10g\n",
		   c[14], c[8], c[9] );
    break;

  case SHOMATE2:
    dnt(1); printf("SHOMATE Polynomial format: 2 zones\n");
    dnt(2); printf("Low  temperature polynomials: %g < T < %g: \n",
		   minTemp, c[0]);
    dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n",
		   c[1],  c[2], c[3], c[4]);
    dnt(2); printf("%17.11g %17.11g %17.11g\n",
		   c[5],  c[6], c[7] );
    dnt(2); printf("High temperature polynomials: %g < T < %g: \n",
		   c[0], maxTemp);
    dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n",
		   c[8],  c[9], c[10], c[11]);
    dnt(2); printf("%17.11g %17.11g %17.11g\n",
		   c[12], c[13], c[14] );
    break;

  case NASA1:
    dnt(1); printf("NASA Polynomial format: 1 zones\n");
    dnt(2); printf("temperature polynomials: %g < T < %g: \n",
		   minTemp, maxTemp);
    dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n",
		   c[1], c[2],  c[3],  c[4]);
    dnt(2); printf("%17.11g %17.11g %17.11g\n",
		   c[5], c[6], c[7] );
    break;

  case SHOMATE1:
    dnt(1); printf("SHOMATE Polynomial format: 1 zones\n");
    dnt(2); printf("temperature polynomials: %g < T < %g: \n",
		   minTemp, maxTemp);
    dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n",
		   c[1],  c[2], c[3], c[4]);
    dnt(2); printf("%17.11g %17.11g %17.11g\n",
		   c[5],  c[6], c[7] );
    break;

  case CONSTANT_CP:
  case SIMPLE:
    dnt(1); printf("CONSTANT_CP format:\n");
    dnt(2); printf(" Valid Range: %g < T < %g: \n",
		   minTemp, maxTemp);
    dnt(2); printf(" at T = %5g K, H  = %17.11g J kmol-1\n", c[0], c[1]);
    dnt(2); printf("                 S  = %17.11g J kmol-1 K-1\n", c[2]);
    dnt(2); printf("                 Cp = %17.11g J kmol-1 K-1\n", c[3]);
    break;

  case MU0_INTERP:
    dnt(1); printf("MU0_POLY format:\n");
    nint = (int) c[0];
    dnt(2); printf(" Valid Range: %g < T < %g: \n",
		   minTemp, maxTemp);
    dnt(2); printf(" T        mu_0\n");
    j = 2;
    for (i =0 ; i < nint; i++) {
      dnt(2); printf(" %g   %17.11g\n", c[j], c[j+1]);
      j += 2;
    }
    break;

  case NASA9:
    nzones = 1;
    dnt(1); printf("NASA9 Polynomial format: %d zones\n", nzones);
    cptr = c;
    dnt(2); printf("%2d:  Temperature polynomials: %g < T < %g: \n",
		   0, minTemp, maxTemp);
    dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n",
		   cptr[3], cptr[4],  cptr[5],  cptr[6]);
    dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n",
		   cptr[7], cptr[8], cptr[9], cptr[10] );
    dnt(2); printf("%17.11g\n", cptr[11] );     
    break;

  case NASA9MULTITEMP:
    nzones = (int) c[0];
    dnt(1); printf("NASA9 Polynomial format: %d zones\n", nzones);
    index = 1;
    for (i = 0; i < nzones; i++) {
      minTemp = c[index];
      maxTemp = c[index+1];
      cptr = c+index+2;
      dnt(2); printf("%2d:  Temperature polynomials: %g < T < %g: \n",
		     i, minTemp, maxTemp);
      dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n",
		     cptr[0], cptr[1],  cptr[2],  cptr[3]);
      dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n",
		     cptr[4], cptr[5], cptr[6], cptr[7] );
      dnt(2); printf("%17.11g\n", cptr[8] );
      index += 11;
    }
    break;
  case PDSS_TYPE:
    dnt(1); printf("Presure Dependent Standard State form\n");
    vpss = dynamic_cast<VPStandardStateTP *>(g_ptr);
    pdss_ptr = vpss->providePDSS(k);
    ptype = pdss_ptr->reportPDSSType();
 
    doublereal minTemp, maxTemp, refPressure;
    pdss_ptr->reportParams(kindex, type, &c[0], minTemp, maxTemp, refPressure);

    dnt(2); printf("PDSS Species index    = %d\n", (int) kindex);
    dnt(2); printf("PDSS min Temperature  = %g\n", minTemp);
    dnt(2); printf("PDSS max Temperature  = %g\n", maxTemp);
    dnt(2); printf("PDSS reference pres   = %g\n", refPressure);
    switch (ptype) {
    case cPDSS_IDEALGAS:
      dnt(2); printf("Ideal Gas PDSS thermo type\n");
      break;
    case  cPDSS_CONSTVOL:
      dnt(2); printf("ConstVol PDSS thermo type\n");
      break;
    case cPDSS_MOLAL_CONSTVOL:
      dnt(2); printf("Molal ConstVol PDSS thermo type\n");
      break;
    case cPDSS_WATER:
      dnt(2); printf("Water PDSS thermo type\n");
      break;
    case cPDSS_MOLAL_HKFT:
      dnt(2); printf("PDSS_MOLAL_HKFT PDSS thermo type\n");
      dnt(3); printf("deltaG_formation_tr_pr = %g\n", c[0]);
      dnt(3); printf("deltaH_formation_tr_pr = %g\n", c[1]);
      dnt(3); printf("Mu0_tr_pr              = %g\n", c[2]);	
      dnt(3); printf("Entrop_tr_pr           = %g\n", c[3]);
      dnt(3); printf("a1                     = %g\n", c[4]);
      dnt(3); printf("a2                     = %g\n", c[5]);
      dnt(3); printf("a3                     = %g\n", c[6]);
      dnt(3); printf("a4                     = %g\n", c[7]);
      dnt(3); printf("c1                     = %g\n", c[8]);
      dnt(3); printf("c2                     = %g\n", c[9]);
      dnt(3); printf("omega_pr_tr            = %g\n", c[10]);
	
      break;
    case cPDSS_UNDEF:
      dnt(2); printf("Undefined PDSS thermo type\n");
      break;
    default:
      printf("unknown species PDSS thermo type %d\n", (int) ptype);
    }
    break;
  default:
    printf("unknown species reference thermo type %d\n", rt);
  }

    
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

