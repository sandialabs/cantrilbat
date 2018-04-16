/**
 * file m1d_m1d_BC_Battery.h
 * Header for class BoundaryCondition and subclasses associated with batteries
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

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! Boundary condition to apply to the current equation that takes into account of the resistance of the current collector
/*!
 *  Boundary condition that is applied on the anode current collector to take into account of the voltage loss in the anode
 *  current collector.
 *
 *  The valueAtTime() function returns the current flux exiting the domain from the anode wire as a function of the voltage at the
 *  anode current collector - anode material block interface.
 *
 *  To set up the boundary condition the thickness of the anode collector is needed along with specification of the electric
 *  potential at the anode, which is normally set to 0. Currently, the resistivity of the anode current collector is 
 *  fixed via a call to the copper resistivity function.
 *
 *  This boundary condition can be used whenever a function for the current flux boundary condition is needed.
 *  This is used as a Robin boundary condition where the current flux expression is needed as a function of the electric potential
 *  at the acc-anode boundary.
 */
class BC_anodeCC: public BoundaryCondition
{
public:

    //! Constructor
    /*!
     *  @param[in]           thickness           Thickness of the current collector in the region where the 
     *                                           anode current collector has the same cross-sectional area as the battery.
     *
     *  @param[in]           anodeCC_volts       Value of the electric potential at the anode wire.
     *                                             Defaults to 0.0;
     */
    BC_anodeCC(double thickness, double anodeCC_volts = 0.0);

    //! Copy constructor
    /*!
     *  @param[in]           right               object to be copied
     *
     *  @return                                  returns a reference to the current object
     */
    BC_anodeCC(const BC_anodeCC& right);

    //! Destructor
    virtual ~BC_anodeCC();

    //! Assignment operator
    /*!
     *  @param[in]           right               object to be copied
     *
     *  @return                                  returns a reference to the current object
     */
    BC_anodeCC& operator=(const BC_anodeCC& right);

    //! Returns the current density exiting the anode (i dot n_ext), which is nominally negative for the discharge cycle of a battery
    /*!
     *  @param[in]           time                Current time
     *  @param[in]           voltsAnode          Electric potential at the anode current collector - anode material block interface
     *  @param[in]           interval            unused.
     *
     *  @return                                  Returns the current density leaving the anode current collector into the anode wire
     *                                             Units:  Amps / m^2
     */
    virtual double valueAtTime(double time, double voltsAnode, int interval) const override;

protected:

    //! Value of the electric potential at the anode wire.
    /*!
     *    Units:   volts 
     */
    double anodeCC_volts_;

    //! Thickness of the current collector in the region where the acc has the same cross-sectional area.
    /*!
     *    Units:  m
     */
    double thickness_;
};

//==================================================================================================================================
//! Boundary condition to apply to the current equation that takes into
//! account of the resistance of the current collector
/*!
 *  Boundary condition that is applied on the cathode current collector to take into account of the voltage loss in the cathode
 *  current collector. The dependent variable is the electric potential at the cathode - cathode current collector boundary.
 *  The potential at the cathode current collector wire is considered fixed, and input as an initial condition.
 *
 *  The valueAtTime() function returns the current flux exiting the domain from the cathode wire as a function of the voltage at the
 *  cathode current collector - cathode material block interface.
 *
 *  To set up the boundary condition the thickness of the cathode collector is needed along with specification of the electric
 *  potential at the anode, which is normally set to 0. Currently, the resistivity of the anode current collector is 
 *  fixed via a call to the aluminum resistivity function.
 *
 *  This boundary condition can be used whenever a function for the current flux at the boundary in terms of the voltages are needed.
 *  This is used as a Robin boundary condition where the current flux expression is needed as a function of the electric potential
 *  at the ccc-cathode boundary.
 */
class BC_cathodeCC: public BoundaryCondition
{

public:

    //! Constructor
    /*!
     *  @param[in]           thickness           Thickness of the current collector in the region where the 
     *                                           anode current collector has the same cross-sectional area as the battery.
     *  @param[in]       extraWireResistance     Extra resistance when the cathode collector is considered to be a wire
     *  @param[in] electrodeWireCrossSectionalArea  Cross-sectional area of the wire when the cathode collector is considered to 
     *                                              be a wire with different area than the battery
     *
     *  @param[in]           cathodeCC_phi      Value of the electric potential at the cathode wire.
     *                                                There is no default
     */
    BC_cathodeCC(double thickness, double extraWireResistance, double electrodeWireCrossSectionalArea,
                 double cathodeCC_phi);

    //! Copy constructor
    /*!
     *  @param[in]           right               object to be copied
     *
     *  @return                                  returns a reference to the current object
     */
    BC_cathodeCC(const BC_cathodeCC& right);

    //! Destructor
    virtual ~BC_cathodeCC();

    //! Assignment operator
    /*!
     *  @param[in]           right               object to be copied
     *
     *  @return                                  returns a reference to the current object
     */
    BC_cathodeCC& operator=(const BC_cathodeCC& right);
    
    //! Returns the current density on a cross-sectional basis exiting the battery from the current collector wire
    /*!
     *  @param[in]           time                Value of the time
     *  @param[in]           voltsCathode        electric potential at the cathode - cathode current collector
     *  @param[in]           interval            Value of the time interval (unused here)
     *
     *  @return                                  Current exiting the battery from the cathode
     *                                             units = amps / m2
     */
    virtual double valueAtTime(double time, double voltsCathode, int interval) const override;

protected:

    //! Electric potential at the cathode current collector wire
    /*!
     *  This is considered as a fixed initial input
     *   Units: volts
     */
    double cathodeCC_phi_;

    //! Thickness of the current collector in the region where the ccc has the same cross-sectional area.
    //! as the battery
    /*!
     *    Units:  m
     */
    double thickness_;

    //! Extra resistance attached to the battery (ohms)
    /*!
     *  This reistance is applied when the cathode current collector is considered a wire
     *    Units = Ohms
     */
    double extraWireResistance_;

    //! Cross sectional area of the cathode wire when it is considered a wire and not an extension of the battery
    /*!
     *    Units: (m2)
     */
    double electrodeWireCrossSectionalArea_;
};
//==================================================================================================================================
//!  Boundary condition to apply to the current equation that takes into
//!  account of the resistance of the current collector
/*!
 *   Mixed boundary condition that is used to represent a battery under a load.
 */
class BC_cathodeCCLoad: public BoundaryCondition
{

public:

    //! Constructor
    /*!
     *  @param[in]           thickness           Thickness of the current collector in the region where the 
     *                                           anode current collector has the same cross-sectional area as the battery.
     *  @param[in]       extraWireResistance     Extra resistance when the cathode collector is considered to be a wire
     *  @param[in] electrodeCrossSectionalArea   Cross-sectional area of the battery
     *                                             Units: m2
     *
     *  @param[in]           cathodeCC_phi       Value of the electric potential at the cathode wire.
     *                                               NOT USED -> NOT NEEDED
     *                                                There is no default
     *  @param[in]           resistanceLoad      Resistance of the load that is in series with the terminals of the battery
     *                                             Units: ohms
     *  @param[in]           voltageLoad         Voltage of the load that is in series with the terminals of the battery
     *                                             Units: ohms                                         
     */
    BC_cathodeCCLoad(double thickness, double extraWireResistance, double electrodeCrossSectionalArea,
		     double cathodeCC_phi, double resistanceLoad, double voltageLoad);

    //! Copy constructor
    /*!
     *  @param[in]           right               object to be copied
     *
     *  @return                                  returns a reference to the current object
     */
    BC_cathodeCCLoad(const BC_cathodeCCLoad& right);


    //! Virtual destructor
    virtual ~BC_cathodeCCLoad();

    //! Assignment operator
    /*!
     *  @param[in]           right               object to be copied
     *
     *  @return                                  returns a reference to the current object
     */
    BC_cathodeCCLoad& operator=(const BC_cathodeCCLoad& right);
    
    //! Returns the current density at the cathode current collector on a cross-sectional of the battery basis
    //! given the value of the electric potential at the cathode - cathode current collector interface
    /*!
     *  @param[in]           time                Current value of the time
     *  @param[in]           phiCathode          electric potential at the current - current collector interface
     *  @param[in]           interval            Interval in the time coordinate.
     *
     *  @return                                  Current density flowing out of cathode current collector
     *                                             Units:  amps / m2
     */
    virtual double valueAtTime(double time, double phiCathode, int interval) const override;

protected:

    //! Thickness of the current collector in the region where the ccc has the same cross-sectional area.
    //! as the battery
    /*!
     *    Units:  m
     */
    double thickness_;

    //! Extra resistance attached to the battery 
    /*!
     *  Units: ohms
     */
    double extraResistance_;

    //! Cross sectional area of domain
    /*!
     *  Units:  m2
     */
    double electrodeCrossSectionalArea_;

    //! Resistance of the load
    /*!
     *   Units:  ohms
     */
    double resistanceLoad_;

 
    //! Voltage drop of the load
    /*!
     *   Units:  volts
     */
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

    //! Assignment operator
    /*!
     *  @param[in]           right               object to be copied
     *
     *  @return                                  returns a reference to the current object
     */
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
