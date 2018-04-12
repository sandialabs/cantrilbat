/**
 * @file m1d_materials.h  Declaration file for functions which provide material properties
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _M1D_MATERIALS_H
#define _M1D_MATERIALS_H

#include  "m1d_ProblemStatementCell.h"

#include "m1d_CanteraElectrodeGlobals.h"



namespace m1d
{
    //! Returns the resistivity of aluminum
    /*!
     *  Aluminum is used as the cathode current collector, usually.
     *
     *    @param t temperature
     *
     *    @return  resistivity (Ohm m)
     */
    double resistivity_aluminum(double T = 293.15);

    //! Returns the resistivity of copper
    /*!
     *  Copper is used as the anode current collector, usually.
     *
     *    @param T temperature
     *
     *    @return  resistivity (Ohm m)
     */
    double resistivity_copper(double T = 293.15);

    //!  returns the artificial compressibility of a liquid electrolyte in Pascals-1.
    /*!
     *
     *   @return artificial compressibility of material in pa-1.
     */
    double artificialCompressibility();
    
}
#endif
