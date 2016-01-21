/*
 * @file         GFCEO_Electrode.cpp
 */
/*
 * Copywrite 2015 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "GFCEO_Electrode.h"


//-----------------------------------------------------------------------------------------------------------------------------------
namespace Cantera
{
//===================================================================================================================================
GFCEO_Electrode::GFCEO_Electrode(doublereal atol) :
   Cantera::ResidJacEval(atol)
{

}
//===================================================================================================================================
GFCEO_Electrode::~GFCEO_Electrode()
{

}
//===================================================================================================================================
GFCEO_Electrode::GFCEO_Electrode(const GFCEO_Electrode& right) :
    Cantera::ResidJacEval()
{
    operator=(right);
}
//===================================================================================================================================
GFCEO_Electrode&  GFCEO_Electrode::operator=(const GFCEO_Electrode& right)
{
    if (this == &right) {
        return *this;
    }
    Cantera::ResidJacEval::operator=(right);

    return *this;
}
//===================================================================================================================================
} // End of namespace Cantera
//-----------------------------------------------------------------------------------------------------------------------------------
