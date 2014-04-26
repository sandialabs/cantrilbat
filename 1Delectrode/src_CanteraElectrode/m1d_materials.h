#ifndef _M1D_MATERIALS_H
#define _M1D_MATERIALS_H


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
    
}
#endif
