/**
 *  @file Electrode_SurfaceRegion.h
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

namespace Zuzax
{
//==================================================================================================================================
//! This class is a derived class used to model phase-change electrodes
/*!
 *  The class is an intermediary, support class. It's main purpose is to house the member data
 *  needed to handle a surface in the middle of two regions that are distributed radially.
 *  The regions can expand and contract through surface reactions.
 *  The surface reactions all occur on interfacial kinetics objects.
 *
 *  (NOT COMPLETED)
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
     * @param[in]            right               Object to be copied
     */
    Electrode_SurfaceRegion(const Electrode_SurfaceRegion& right);

    //! Assignment operator
    /*!
     *  @param[in]           right               object to be copied
     *
     *  @return                                  Returns a reference to a current object
     */
    Electrode_SurfaceRegion& operator=(const Electrode_SurfaceRegion& right);

    //! Duplicator function
    /*!
     *  Duplicate the current Electrode object, returning a base Electrode pointer
     *
     *  @return                                   Returns a duplicate of the current object as a base class pointer
     */
    virtual Electrode* duplMyselfAsElectrode() const override;

    //! Return the type of electrode
    /*!
     *  Returns the enum type of the electrode. This is used in the factory routine.
     *
     *  @return                                  Returns an enum type, called Electrode_Types_Enum
     */
    virtual Electrode_Types_Enum electrodeType() const override;

    //! Create the electrode model
    /*!
     *  (virtual from Electrode  - onion Out)
     *
     *    This is one of the most important routines. It sets up the electrode's internal structures
     *    After the call to this routine, the electrode should be internally ready to be integrated and reacted.
     *    It takes its input from an ELECTRODE_KEY_INPUT object which specifies the setup of the electrode
     *    object and the initial state of that object.
     *    The routine works like an onion Out initialization. The parent object is initialized before the
     *    child. This means the child object first calls the parent, before it does its own initializations.
     *
     *    There are some virtual member functions that won't work until this routine is called. That's because
     *    the data structures won't be set up for base and child Electrode objects until this is called.
     *
     *  @param[in]           ei                  BASE ELECTRODE_KEY_INPUT pointer object. Note, it must have the correct child class
     *                                           for the child electrode object.
     *
     *  @return                                  Returns zero if successful, and -1 if not successful.
     */
    virtual int electrode_model_create(ELECTRODE_KEY_INPUT* ei) override;

    //! Specify initial conditions from an input file
    /*!
     *    @param[in]         ei                  pointer to a parent Structure containing the input file. The structure will be dynamically
     *                                           cast to a type  ELECTRODE_SurfaceRegion_KEY_INPUT type within the program.
     *                                           It is a fatal error not to be able to do the dynamic cast.
     *
     *    @return                                Returns a zero if everything is ok. Anything else is a fatal error.
     */
    virtual int setInitialConditions(ELECTRODE_KEY_INPUT* ei) override;

    //! Initialize the sizes
    void init_sizes();

    //! Calculate the number of equations that will be solved during the nonlinear solver step.
    /*!
     *  (virtual from Electrode_Integrator)
     *  All classes which inherit from this routine must have a class that determines this value.
     *
     *  @return                                  Returns the number of unknowns in the nonlinear problem and time-stepping problem.
     */
    virtual size_t nEquations_calc() const override;

    //! Get all of the reaction rates and parameters from Zuzax
    /*!
     *  @param[out]          justBornMultiSpecies Vector of phases which are just born
     */
    void extractInfoJustBorn(std::vector<size_t>& justBornMultiSpecies);

    //! Print conditions of the electrode for the current integration step to stdout
    /*!
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values
     */
    virtual void printElectrode(int pSrc = 1, bool subTimeStep = false) override;

    //! Print condition of a phase in the electrode
    /*!
     *  @param[in]           iph                 Print the phase
     *  @param[in]           pSrc                Print Source terms that have occurred during the step from the initial_initial
     *                                           to the final_final time. The default is to print out the source terms
     *  @param[in]           subTimeStep         Print out conditions from the most recent subTimeStep and not the global
     *                                           time step. The default is to print out the global values
     */
    virtual void printElectrodePhase(size_t iph, int pSrc = 1, bool subTimeStep = false) override;

    //! The internal state of the electrode must be kept for the initial and
    //! final times of an integration step.
    /*!
     *  (virtual function from Electrode)
     *  This function advances the initial state to the final state that was calculated
     *  in the last integration step.
     *
     *  @param[in]            Tinitial            This is the New initial time. This time is compared against the "old"
     *                                           final time, to see if there is any problem.
     *  @param[in]            doAdvancementAlways      Do the advancement always, even if the Tinitial value is equal to t_init_init_
     *
     *  @return                                  Returns true if the time step is reset to t_init_init.
     */
    virtual bool resetStartingCondition(double Tinitial, bool doAdvancementAlways = false) override;

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
    virtual void updateState() override;

protected:

    //! Pointer to the reacting surface domain that is located on this surface
    ReactingSurDomain* rsd_;

    //! Surface index of the surface phase
    size_t indexSurfaceRegion_;

    //! Number of Surface species defined on the surface
    size_t numSurfaceSpecies_;

    //! Number of Surface phases defined on the surface
    size_t numSurfacePhases_;

    //! Define the number of species that are defined to have radially distributed distributions
    //! in the domain to the left 
    /*!
     *  The left is defined as the inside of the particle.
     */
    size_t numKRSpeciesLeft_;

    //! Global index of the surface node
    size_t indexRNode_;

    //! Node position of the surface - final_final
    double rnodePos_final_final_;

    //! Node position of the surface - final
    double rnodePos_final_;

    //! Node position of the surface - init
    double rnodePos_init_;

    //! Node position of the surface- init_init
    double rnodePos_init_init_;

    //! Reference radius of the node at the surface
    double rRefPos_final_final_;

    //! Reference radius at the node at the surface - local final value
    double rRefPos_final_;

    //! Reference radius at the node at the surface - local init value
    double rRefPos_init_;

    //! Reference radius at the node at the surface - global init value
    double rRefPos_init_init_;

};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
