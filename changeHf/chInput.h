
/**
 *  @file chInput.h
 *
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

/*
 * $Id: chInput.h 508 2013-01-07 22:54:04Z hkmoffa $
 */

#ifndef CHINPUT_H
#define CHINPUT_H

#include "BlockEntry.h"
#include "mdp_allo.h"
#include <string>

#define UNITS_KCAL_CGS 0
#define UNITS_KJOULE   1
#define UNITS_CGS      2

#define HFVALUE_UNSET -1000.

void setup_input(BEInput::BlockEntry *cf);
void process_input(BEInput::BlockEntry *cf, FILE *cmdptr);

bool doubleEqual(double d1, double d2, double atol = 1.0E-100);

class IOoptions {
public:
    /*
     * Constructor():
     */
    IOoptions() : 

	OutputUnits(UNITS_KCAL_CGS),
        DeltaValue(0.0),
        HfValue(HFVALUE_UNSET)
	{
	}
    /*
     * Destructor():
     */
    ~IOoptions() {
    }

    /*
     * File name
     */
    std::string FileName;
    std::string DestFileName;

    /*
     * OutputUnits
     */
    int OutputUnits;

    /*
     * SpeciesName
     */
    std::string SpeciesName;

    /*
     *  Delta value of the heat of formation
     */
    double DeltaValue;

    /*
     * Absolute value of the heat of formation to set it to
     */
    double HfValue;
};

extern IOoptions IOO;


extern int DebugPrinting;

#endif
