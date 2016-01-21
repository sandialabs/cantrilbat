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
GFCEO_Electrode::GFCEO_Electrode(Electrode* ee, doublereal atol, int iOwn) :
   Cantera::ResidJacEval(atol),
   ee_(ee),
   iOwnObject_(iOwn)
{

}
//===================================================================================================================================
GFCEO_Electrode::~GFCEO_Electrode()
{
    if (iOwnObject_) {
        delete ee_;
	ee_ = 0;
    }
}
//===================================================================================================================================
GFCEO_Electrode::GFCEO_Electrode(const GFCEO_Electrode& right) :
    Cantera::ResidJacEval()
{
    operator=(right);
}
//===================================================================================================================================
GFCEO_Electrode& GFCEO_Electrode::operator=(const GFCEO_Electrode& right)
{
    if (this == &right) {
        return *this;
    }
    Cantera::ResidJacEval::operator=(right);

    if (right.iOwnObject_) {
        delete ee_;
        ee_ = right.ee_->duplMyselfAsElectrode();
        iOwnObject_ = 1;
    } else {
        ee_ = right.ee_;
        iOwnObject_ = 0;
    }

    return *this;
}
//===================================================================================================================================
Electrode& GFCEO_Electrode::electrode()
{
    return *ee_;
}
//===================================================================================================================================
void GFCEO_Electrode::assertOwnership()
{
     if (! iOwnObject_) {
        Electrode* e = ee_->duplMyselfAsElectrode();
        ee_ = e;
        iOwnObject_ = 1;
    } 
}
//===================================================================================================================================
} // End of namespace Cantera
//-----------------------------------------------------------------------------------------------------------------------------------
