/**
 *  @file m1d_materials.cpp   Definitions for miscellaneous functions which provide material properties
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_materials.h"


namespace m1d
{
extern ProblemStatementCell *PSCinput_ptr;
//====================================================================================================
double resistivity_aluminum(double T)
{
    return 2.82E-8 * (1.0 + 0.0039 * (T - 293.15));
}
//====================================================================================================
double resistivity_copper(double T)
{
    return 1.68E-8 * (1.0 + 0.003862 * (T - 293.15));
}
//====================================================================================================
double artificialCompressibility()
{
#ifdef DEBUG_DARCY
   double beta = 0.69 / 1.01E5;
#else
   double beta = PSCinput_ptr->artificialCompressibilityInvAtm_ / OneAtm;
#endif
   return beta;

}
}
//====================================================================================================
