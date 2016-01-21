/*
 * @File GFCEO_Electrode_Integrator.h coupling 
 */
/*
 * Copywrite 2015 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _GFCEO_ELECTRODE_H
#define _GFCEO_ELECTRODE_H


#include "Electrode.h"
#include "cantera/numerics/ResidJacEval.h"
//-----------------------------------------------------------------------------------------------------------------------------------
namespace Cantera
{

//===================================================================================================================================
//! This class is a derived class used to carry out fully coupled simulations
/*!
 * Complete problem statement
 *
 */
class GFCEO_Electrode : public Cantera::ResidJacEval
{
public:

    // ---------------------------------------------------------------------------------------------
    // ----------------------- BASIC SETUP ROUTINES  -----------------------------------------------
    // ---------------------------------------------------------------------------------------------

    //! Constructor
    /*!
     *  @param[in]      atol                     Default value for the absolute tolerance
     */
    GFCEO_Electrode(doublereal atol = 1.0E-13);

    //! Destructor
    virtual ~GFCEO_Electrode();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    GFCEO_Electrode(const GFCEO_Electrode& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    GFCEO_Electrode& operator=(const GFCEO_Electrode& right);

};

} // End namespace Cantera
//-----------------------------------------------------------------------------------------------------------------------------------
#endif

