/**
 *  @file Electrode_RadialDiffRegions.h
 *
 *
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program  Electrode_RadialDiffRegions.h  
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_RADIALDIFFREGIONS_H
#define _ELECTRODE_RADIALDIFFREGIONS_H

#include "Electrode_Integrator.h"

#include "Electrode_RadialRegion.h"
#include "Electrode_SurfaceRegion.h"

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{

class Electrode_RadialRegion;
class Electrode_SurfaceRegion;

//==================================================================================================================================
//! Extra Input for models with distributed radial diffusion regions
/*!
 *   This model is used for Electrode_SimpleDiff
 */
class ELECTRODE_RadialDiffRegions_KEY_INPUT : public ELECTRODE_KEY_INPUT
{
public:

    //! Constructor
    ELECTRODE_RadialDiffRegions_KEY_INPUT(int printLvl = 0);

    //! Destructor
    virtual ~ELECTRODE_RadialDiffRegions_KEY_INPUT();

    ELECTRODE_RadialDiffRegions_KEY_INPUT(const ELECTRODE_RadialDiffRegions_KEY_INPUT& right);

    ELECTRODE_RadialDiffRegions_KEY_INPUT& operator=(const ELECTRODE_RadialDiffRegions_KEY_INPUT& right);

    //!  First pass through the child setup system
    /*!
     *    Typically we will fill in all vectors that depend on the value of numRegions_ in this
     *    pass.
     *  @param cf    Pointer to the BlockEntry record
     */
    void setup_input_child1(BEInput::BlockEntry* cf);

    //!  Second pass through the child setup system
    /*!
     *    Typically we will fill in all vectors that depend on the value of numRegions_ in this
     *    pass.
     *  @param cf    Pointer to the BlockEntry record
     */
    void setup_input_child2(BEInput::BlockEntry* cf);

    virtual void post_input_child2(BEInput::BlockEntry* cf);

    virtual void setup_input_child3(BEInput::BlockEntry* cf);

    virtual void post_input_child3(BEInput::BlockEntry* cf);

    //   MEMBER DATA

    //! Number of regions in the model
    /*!
     *  The number of regions refers to the distribution of phases within the radially symmetric
     *  model. Within each region a certain number of phases are distributed radially. They are
     *  connected to other regions by interfacial surface regions.
     */
    int numRegions_;
 
    //! Number of region domains actually entered by the user.
    int numRegionsEntered_;

    //! Solid state diffusion model identification
    /*!
     *   0 = There is no diffusion model (default)
     *   1 = There is solid state diffusion limitations. The model assumes pseudo-steady state
     *       The limiting reacting surface location is determined by the  locationOfReactingSurface_ model.
     *       The details of the model is given in the writeup.
     */
    int solidDiffusionModel_;

    //!  Model for the formulation of the diffusive flux
    /*!
     *   0 = diffusive flux include the activity coefficient
     *   1 = diffusive flux doesn't include the activity coefficient
     */
    int diffusiveFluxModel_;

    //! Number of radial cells in each region
    /*!
     *  Number of cells in each region. Note we use a control volume formulation. So, two cells would mean
     *  a region with only two nodes, one on each exterior surface.
     *
     *  Length = numRegions_
     *  Default = 5
     *  units = none.
     */
    std::vector<int> numRadialCellsRegions_;

    //! Diffusion coefficient in the outer region of a two region spherical model
    /*!
     *  We identify the solids by zones, which are related to regions. Generally, when the extent of reaction is
     *  in the zeroeth region, we are using diffusionCoeffRegions_[1], i.e., the first zone, because that is
     *  zone 1 in the model. Zone 0 would be diffusion in the inner core, whose value we never
     *  actually need. It's actually conceptually simpler this way.
     *
     *  Length = numRegions_ + 1
     *  Units = m**2 s-1
     */
    std::vector<double> diffusionCoeffRegions_;

    //! Unused atm
    std::vector<double> rxnPerturbRegions_;

    //! Vector of  classes that are used to store the information about the regions
    std::vector<ELECTRODE_RadialRegion_KEY_INPUT> rregions_;
};

//==================================================================================================================================
//! This class is a derived class used to model phase - change electrodes
/*!
 *  The class is an intermediary, support class. It's main purpose is to house the member data
 *  needed to handle the disretization of one or more regions in the radial direction.
 *  Diffusion of material through each region is allowed.
 *  The regions can expand and contract through surface reactions.
 *  The surface reactions all occur on interfacial kinetics objects.
 *
 *  Child classes will actually complete each of the physical problems.
 *
 */
class Electrode_RadialDiffRegions : public Electrode_Integrator
{
public:
    //! Constructor
    Electrode_RadialDiffRegions();

    //! Destructor
    virtual  ~Electrode_RadialDiffRegions();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    Electrode_RadialDiffRegions(const Electrode_RadialDiffRegions& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    Electrode_RadialDiffRegions& operator=(const Electrode_RadialDiffRegions& right);

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
     *  @return Returns an enum type, called   Electrode_Types_Enum
     */
    virtual Electrode_Types_Enum electrodeType() const;

    //! create the electrode model
    int electrode_model_create(ELECTRODE_KEY_INPUT* ei);


    //! Specify initial conditions from an input file
    /*!
     *    @param[in]         ei                  pointer to a parent Structure containing the input file. The structure will be dynamically
     *                                           cast to a type  ELECTRODE_RadialDiffRegions_KEY_INPUT type within the program.
     *                                           It is a fatal error not to be able to do the dynamic cast.
     *
     *    @return                                Returns a zero if everything is ok. Anything else is a fatal error.
     */
    virtual int setInitialConditions(ELECTRODE_KEY_INPUT* ei) override;

    //! Initialize the sizes
    void init_sizes();

    //! Initialize the distribution of species evenly across the radius
    /*!
     *  This routine takes the spMoles values and spreads the values evenly in the radial
     *  direction so that the radial concencentrations are constant
     */
    void initializeAsEvenDistribution();

    //! Get all of the reaction rates and parameters from Cantera
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
     *  @param iph           Print the phase
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values
     */
    virtual void printElectrodePhase(size_t iph, int pSrc = 1, bool subTimeStep = false) override;

    //! The internal state of the electrode must be kept for the initial and
    //! final times of an integration step.
    /*
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

    // ------------------------------------------ D A T A -------------------------------------------------------------
protected:
    //! Number of Radial diffusion Regions
    int numRadialRegions_;

    //! List of Radial Regions
    std::vector<Electrode_RadialRegion*> RadialRegionList_;

    //! List of Surface Regions
    std::vector<Electrode_SurfaceRegion*> SurfaceRegionList_;

    //! Molar density of the solid phase in each cell under reference conditions
    /*!
     *  This is the molar volume of the first species in the mechanism at the
     *  starting temperature and pressure.
     *
     *  units are kmol m-3.
     */
    double MolarVolume_Ref_;

    //! Node position of the mesh - final_final
    std::vector<double> rnodePos_final_final_;

    //! Node position of the mesh - final
    std::vector<double> rnodePos_final_;

    //! Node position of the mesh - init
    std::vector<double> rnodePos_init_;

    //! Node position of the mesh - init_init
    std::vector<double> rnodePos_init_init_;

    //! Reference radius at the right cell boundary - global final value
    std::vector<double> rRefPos_final_final_;

    //! Reference radius at the right cell boundary - local final value
    std::vector<double> rRefPos_final_;

    //! Reference radius at the right cell boundary - local init value
    std::vector<double> rRefPos_init_;

    //! Reference radius at the right cell boundary - global init value
    std::vector<double> rRefPos_init_init_;

    //!  Spline system for the nodal equations
    /*!
     *   These factors are the fraction of the exterior node radius that the
     *   current node possesses.
     */
    std::vector<double> fracNodePos_;

};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------

#endif

