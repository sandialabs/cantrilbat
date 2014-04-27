/**
 * file m1d_BoundaryCondition.h
 * 
 * Header for class BoundaryCondition and subclasses.
 */

/*  $Author: hkmoffa $
 *  $Revision: 540 $
 *  $Date: 2013-02-27 15:18:26 -0700 (Wed, 27 Feb 2013) $
 *
 */
// Copyright 2010 Sandia National Laboratories
#ifndef M1D_BC_BATTERY_H
#define M1D_BC_BATTERY_H

#include "m1d_BoundaryCondition.h"

namespace m1d {

//!  Boundary condition to apply to the current equation that takes into
//!  account of the resistance of the current collector
class BC_anodeCC: public BoundaryCondition
{

public:

    BC_anodeCC(double thickness, double anodeCC_volts);
    BC_anodeCC(const BC_anodeCC& right);
    virtual ~BC_anodeCC();
    BC_anodeCC& operator=(const BC_anodeCC& right);

    virtual double valueAtTime(double time, double voltsAnode, int interval);

protected:

    double anodeCC_volts_;

    double thickness_;
};

} //namespace m1d

#endif // M1D_BOUNDARYCONDITION
