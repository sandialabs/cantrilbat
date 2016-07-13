/*
 * $Id: Electrode_SurfaceRegion.h 298 2012-08-08 20:15:48Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_SURFACEREGION_H
#define _ELECTRODE_SURFACEREGION_H

#include "Electrode_Integrator.h"

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{


//! This class is a derived class used to model phase-change electrodes
/*!
 *  The class is an intermediary, support class. It's main purpose is to house the member data
 *  needed to handle the disretization of one region in the radial direction.
 *  Diffusion of material through the region is allowed.
 *  The regions can expand and contract through surface reactions.
 *  The surface reactions all occur on interfacial kinetics objects.
 *
 *  Child classes will actually complete each of the physical problems.
 *
 */
class Electrode_SurfaceRegion : public Electrode_Integrator
{
public:
    //! Constructor
    Electrode_SurfaceRegion();

    //! Destructor
    virtual  ~Electrode_SurfaceRegion();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    Electrode_SurfaceRegion(const Electrode_SurfaceRegion& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    Electrode_SurfaceRegion& operator=(const Electrode_SurfaceRegion& right);

    //! Return the type of electrode
    /*!
     *  Returns the enum type of the electrode. This is used in the factory routine.
     *
     *  @return Returns an enum type, called   Electrode_Types_Enum
     */
    virtual Electrode_Types_Enum electrodeType() const;

    //! create the electrode model
    int electrode_model_create(ELECTRODE_KEY_INPUT* ei);

    //! Specify initial conditions from an input file
    /*!
     *    @param ei pointer to a parent Structure containing the input file. The structure will be dynamically
     *              cast to a type  ELECTRODE_SurfaceRegion_KEY_INPUT type within the program.
     *              It is a fatal error not to be able to do the dynamic cast.
     *
     *    @param  Return flag. Returns a zero if everything is ok. Anything else is a fatal error.
     */
    virtual int setInitialConditions(ELECTRODE_KEY_INPUT* ei);

    //! Initialize the sizes
    void init_sizes();

    //! Get all of the reaction rates and parameters from Cantera
    void extractInfoJustBorn(std::vector<int>& justBornMultiSpecies);


    //! Print conditions of the electrode for the current integration step to stdout
    /*!
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values
     */
    virtual void printElectrode(int pSrc = 1, bool subTimeStep = false);

    //! Print condition of a phase in the electrode
    /*!
     *  @param iPhase        Print the phase
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values
     */
    virtual void printElectrodePhase(int iPhase, int pSrc = 1, bool subTimeStep = false);

    //! The internal state of the electrode must be kept for the initial and
    //! final times of an integration step.
    /*!
     *  (virtual function from Electrode)
     *  This function advances the initial state to the final state that was calculated
     *  in the last integration step.
     *
     * @param Tinitial   This is the New initial time. This time is compared against the "old"
     *                   final time, to see if there is any problem.
     * @param doResetAlways  Do the reset always, even if the Tinitial value is equal to t_init_init_
     */
    virtual void  resetStartingCondition(double Tinitial, bool doResetAlways = false);


    //! Take the state (i.e., the final state) within the Electrode_Model and push it down
    //! to the ThermoPhase objects and propogate it to all other aspects of the final state
    /*!
     *  (virtual function from Electrode)
     *  This virtual function should be written so that child routines are not called by parent routines.
     *
     *  We take the values of spMoles_final_[] and propagate them down to the ThermoPhase
     *  objects in the electrode.
     *
     *  We also take the state of the electrode as described by the mole numbers and mole fractions
     *  and calculate the geometrical details involved with the model. This includes the radii and
     *  thicknesses of the regions and the existence of the annular regions with their associated boolean flags.
     *
     *  All of these properties are defined for the _final_ state.
     *
     *  Thus, this is the main routine that reconciles all of the state information within the object.
     *  At the end of this routine, all aspects of the final state are consistent with each other.
     *
     *  prerequisites: The object must have been already created.
     *
     */
    virtual void updateState();




protected:

    //! Pointer to the reacting surface domain that is located on the surface
    ReactingSurDomain* rsd_;

    //! Surface index of the surface
    int indexSurfaceRegion_;

    //! Number of Surface species defined on the surface
    int numSurfaceSpecies_;

    //! Number of Surface phases defined on the surface
    int numSurfacePhases_;

    //! Define the number of species that are defined to have radially distributed distributions
    //! in the domain to the left
    /*!
     *
     */
    int numKRSpeciesLeft_;


    //! global index of the surface node
    int indexRNode_;

    //! Node position of the surface - final_final
    doublereal rnodePos_final_final_;

    //! Node position of the surface - final
    doublereal rnodePos_final_;

    //! Node position of the surface - init
    doublereal rnodePos_init_;

    //! Node position of the surface- init_init
    doublereal rnodePos_init_init_;

    //! Reference radius of the node at the surface
    doublereal rRefPos_final_final_;

    //! Reference radius at the node at the surface - local final value
    doublereal rRefPos_final_;

    //! Reference radius at the node at the surface - local init value
    doublereal rRefPos_init_;

    //! Reference radius at the node at the surface - global init value
    doublereal rRefPos_init_init_;


};

}


#endif
/*****************************************************************************/
