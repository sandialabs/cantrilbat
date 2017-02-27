/**
 * file m1d_BoundaryCondition.h
 * 
 * Header for class BoundaryCondition and subclasses.
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */


#ifndef M1D_BC_BATTERY_H
#define M1D_BC_BATTERY_H

#include "m1d_BoundaryCondition.h"

namespace m1d {
//==================================================================================================================================
//!  Boundary condition to apply to the current equation that takes into
//!  account of the resistance of the current collector
class BC_anodeCC: public BoundaryCondition
{

public:

    BC_anodeCC(double thickness, double anodeCC_volts);

    BC_anodeCC(const BC_anodeCC& right);

    virtual ~BC_anodeCC();

    BC_anodeCC& operator=(const BC_anodeCC& right);

    virtual double valueAtTime(double time, double voltsAnode, int interval) const override;

protected:

    double anodeCC_volts_;

    double thickness_;
};
//==================================================================================================================================
//!  Boundary condition to apply to the current equation that takes into
//!  account of the resistance of the current collector
class BC_cathodeCC: public BoundaryCondition
{

public:

    BC_cathodeCC(double thickness, double extraResistance, double electrodeCrossSectionalArea,
                 double cathodeCC_volts);

    BC_cathodeCC(const BC_cathodeCC& right);

    virtual ~BC_cathodeCC();

    BC_cathodeCC& operator=(const BC_cathodeCC& right);
    
    //! Returns the current on a cross-sectional basis
    /*!
     *  @return  units = amps / m2
     */
    virtual double valueAtTime(double time, double voltsCathode, int interval) const override;

protected:

    double cathodeCC_volts_;

    double thickness_;

    //! Extra resistance attached to the battery (ohms)
    double extraResistance_;

    //! Cross sectional area of domain
    /*!
     *  (m2)
     */
    double electrodeCrossSectionalArea_;
};

//==================================================================================================================================
//!  Boundary condition to apply to the current equation that takes into
//!  account of the resistance of the current collector
class BC_cathodeCCLoad: public BoundaryCondition
{

public:

    BC_cathodeCCLoad(double thickness, double extraResistance, double electrodeCrossSectionalArea,
		     double cathodeCC_volts, double resistanceLoad, double voltageLoad);
    BC_cathodeCCLoad(const BC_cathodeCCLoad& right);
    virtual ~BC_cathodeCCLoad();
    BC_cathodeCCLoad& operator=(const BC_cathodeCCLoad& right);
    
    //! Returns the current on a cross-sectional basis
    /*!
     *  @return  units = amps / m2
     */
    virtual double valueAtTime(double time, double voltsCathode, int interval) const override;

protected:

    double cathodeCC_volts_;

    double thickness_;

    //! Extra resistance attached to the battery (ohms)
    double extraResistance_;

    //! Cross sectional area of domain
    /*!
     *  (m2)
     */
    double electrodeCrossSectionalArea_;

    double resistanceLoad_;
    double voltageLoad_;
};

//==================================================================================================================================

//!  Boundary condition to apply  a heat transfer flux formulation
class BC_heatTransfer: public BoundaryCondition
{
public:

    BC_heatTransfer(double tranCoeff, double TempRef, double electrodeCrossSectionalArea);
    BC_heatTransfer(const BC_heatTransfer &r);
    virtual ~BC_heatTransfer();
    BC_heatTransfer& operator=(const BC_heatTransfer& r);

    //! Returns the heat transfered out of the domain
    /*!
     *   Units = Watts / m2 / K
     */
    virtual double valueAtTime(double time, double tempCollector, int interval) const override;

protected:

    //! temperature of heat bath (Kelvin)
    double tempRef_;

    //! Heat transfer coefficient
    /*!
     *     Units = Watts / m2 / K
     */
    double tranCoeff_;

    //! cross sectional area
    /*!
     *  units = m2  Note -> this is not needed
     */
    double electrodeCrossSectionalArea_;
};

//==================================================================================================================================
} 
//----------------------------------------------------------------------------------------------------------------------------------
#endif // M1D_BOUNDARYCONDITION
