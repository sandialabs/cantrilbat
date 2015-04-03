/*
 * m1d_porousLiKCl_dom1D.h
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 */
#include "m1d_defs.h"
#include "m1d_cellTmps_PorousFlow.h"

namespace m1d {


//=====================================================================================================================
Value_CellRL::Value_CellRL(double val) :
    center(val), right(val), left(val)
{
}
//=====================================================================================================================
Value_CellRL::Value_CellRL(const Value_CellRL &r) :
    center(r.center), 
    right(r.right), 
    left(r.left)
{
}
//=====================================================================================================================
Value_CellRL& Value_CellRL::operator=(const Value_CellRL &r)
{
    if (this == &r) {
	return *this;
    }
    center = r.center; 
    right = r.right; 
    left = r.left;
    return *this;
}
//=====================================================================================================================
void Value_CellRL::fillNextRight(const Value_CellRL &valCellLeft, double newRightVal)
{
    center = valCellLeft.right;
    left   = valCellLeft.center;
    right  = newRightVal;
}
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
valCellTmps::valCellTmps() :
    AxialVeloc(0.0),
    Temperature(-1.0)
{
}
//=====================================================================================================================
valCellTmps::valCellTmps(const valCellTmps &r) :
    AxialVeloc(r.AxialVeloc),
    Temperature(r.Temperature)
{
}
//=====================================================================================================================
valCellTmps::~valCellTmps()
{
}
//=====================================================================================================================
valCellTmps& valCellTmps::operator=(const valCellTmps &r)
{
    if (this == &r) {
        return *this;
    }
    AxialVeloc = r.AxialVeloc;
    Temperature = r.Temperature;

    return *this;
}
//=====================================================================================================================
}
