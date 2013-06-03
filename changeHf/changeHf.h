       


/**
 *  @file changeHf.h
 *
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

/*
 * $Id: changeHf.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef CHANGEHF_H
#define CHANGEHF_H

#include <string>
using namespace std;

#include "chInput.h"

const double R_kcalmol = 8.314472E7 /  4.184E7 * 1.0E-3;
const double R_Jgmol   = 8.314472;
const double R_kJgmol  = 8.314472 * 1.0E-3;

struct UnitsIO {
    string sGibbs;
    double mGibbs;
    string sEntropy;
    double mEntropy;
                   
    void setup(int unit) {
        if (unit == UNITS_KCAL_CGS) {
          sGibbs = "kcal/gmol";
          mGibbs = R_kcalmol;
          sEntropy = "cal/mol*K";
          mEntropy = R_kcalmol * 1.0E3;
        } else if (unit == UNITS_KJOULE) {
	  sGibbs   = " kJ/gmol ";
          mGibbs   = R_kJgmol;                            
          sEntropy = "J/gmolK";
          mEntropy = R_kJgmol * 1.0E3;
	  
        } else {
          printf("ERROR unknown units");
          exit(-1);
        }
    }
};
extern UnitsIO UIO;

#endif
