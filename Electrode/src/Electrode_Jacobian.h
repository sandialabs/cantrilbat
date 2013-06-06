/**
 * Predict the stability of a single phase that is currently zereod.
 *
 *  This routine fills in an estimate for the solution
 *  Return 1 if the phases are stable and 0 if they are not
 *
 */
/*
 * $Id: Electrode_PhaseStability.h 496 2013-01-07 21:15:37Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_EQUILIBRIUM_H
#define _ELECTRODE_EQUILIBRIUM_H

#include "cantera/equilibrium.h"
#include "cantera/thermo/FixedChemPotSSTP.h"

#include "Electrode.h"

#include <string>
#include <vector>

namespace Cantera {

//! Class which determines the stability of phases due to kinetics
/*!
 *  Note, I believe I can make this class simpler
 *
 */
class Electrode_Jacobian {

public:
    //! Constructor
    /*!
     *
     * @param elect  Electrode object pertaining to the stability problem
     *               The electrode object assigns this object as a "friend"
     */
    Electrode_Jacobian(Electrode* elect);

    //! Destructor
    virtual ~Electrode_Jacobian();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    Electrode_Jacobian(const Electrode_Jacobian& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    Electrode_Jacobian& operator=(const Electrode_Jacobian& right);



protected:

    //! This is a reference to the friend object, where we will pull most of the data from.
    Electrode* ee_;


public:
    //! Print level for input to vcs routines within this object
    /*!
     *    This is a public member so that it can be manipulated
     */
    int printLvl_;

};

}
#endif
/*****************************************************************************/

