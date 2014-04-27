/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_materials.h"

namespace m1d
{
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





}
//====================================================================================================
