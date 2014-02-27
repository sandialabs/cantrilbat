/*
 * $Id: Electrode_RadialRegion.h 298 2012-08-08 20:15:48Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_RADIALREGION_H
#define _ELECTRODE_RADIALREGION_H


#include "Electrode.h"
#include "Electrode_Integrator.h"
#include "cantera/integrators.h"
#include "cantera/numerics/ResidJacEval.h"


namespace Cantera
{

class ELECTRODE_RadialRegion_KEY_INPUT : public ELECTRODE_KEY_INPUT
{
public:

    //! Constructor
    ELECTRODE_RadialRegion_KEY_INPUT(int printLvl = 0);

    //! Destructor
    virtual ~ELECTRODE_RadialRegion_KEY_INPUT();

    ELECTRODE_RadialRegion_KEY_INPUT(const ELECTRODE_RadialRegion_KEY_INPUT& right);

    ELECTRODE_RadialRegion_KEY_INPUT& operator=(const ELECTRODE_RadialRegion_KEY_INPUT& right);

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

    //! Index of the region in the model
    /*!
     * 
     */
    int indexRegion_;

    //! Phase indeces of the solid phases comprising the species that are radially distributed
    //! within the current region
    /*!
     *   This is a very important indexing array.
     *   There are numSPhase_ of these.
     *
     *  Length =  numSPhase_
     *  access method = phaseIndeciseKRsolidPhases_[iSPphase_] = iPhaseIndex;
     *
     *  The iPhaseIndex is the index within the Electrode object.
     */
    std::vector<int> phaseIndeciseKRsolidPhases_;

    //! Solid state diffusion model identification
    /*!
     *   0 = There is no diffusion model (default)
     *   1 = There is solid state diffusion limitations. The model assumes pseudo-steady state
     *       The limiting reacting surface location is determined by the  locationOfReactingSurface_ model.
     *       The details of the model is given in the writeup.
     */
    int solidDiffusionModel_;

    //!  String name of the phase that is represented by this region
    /*!
     *   Note will expand to multiple phases here soon
     */
    std::string phaseName_;


    //! Diffusion coefficient in the outer region of a two region spherical model
    /*!
     *  We identify the solids by Zones, which are related to regions. Generally, when the extent of reaction is
     *  in the zeroeth region, we are using diffusionCoeffRegions_[1], i.e., the first zone, because that is
     *  zone 1 in the model. Zone 0 would be diffusion in the inner core, whose value we never
     *  actually need. It's actually conceptually simpler this way.
     *
     *  Length = numRegions_ + 1
     *  Units = m**2 s-1
     */

    std::vector<double> diffusionCoeffSpecies_;   

    //! Default diffusion coefficients for species within the region
    /*!
     *  Defaults to 1.0E-12 m2/s = 1.0E-8 cm2/s
     */
    double defaultDiffusionCoeff_;
};

    class Electrode_RadialDiffRegions;

//! This class is a derived class used to model phase-change or intercalating electrodes
/*!
 *  The class is an intermediary, support class. It's main purpose is to house the member data
 *  needed to handle the disretization of one region in the radial direction.
 *  Diffusion of material through the region is allowed.
 *  The regions can expand and contract through surface reactions.
 *  The surface reactions all occur on interfacial kinetics objects.
 *
 *  Each region is separated from the next by a surface reaction object
 *
 *  Child classes will actually complete each of the physical problems.
 *
 *
 *   Unknowns in this problem are:
 *
 *           Rlattice
 *           concKRSpecies_Cell_final_[]
 *
 *  Because it's a supporting class, some of the electrode functionality must be turned off.
 *
 *  The basic idea is that the object inherits from the full Electrode object, which contains all
 *  of the phases. However, in the particular region, there exists only a subset of the phases.
 *  Also the region lies between a starting radius and an ending radius.
 *
 */
class Electrode_RadialRegion : public Cantera::Electrode_Integrator
{
public:
    //! Constructor
    Electrode_RadialRegion();

    //! Destructor
    virtual  ~Electrode_RadialRegion();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    Electrode_RadialRegion(const Electrode_RadialRegion& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    Electrode_RadialRegion& operator=(const Electrode_RadialRegion& right);

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
     *              cast to a type  ELECTRODE_RadialRegion_KEY_INPUT type within the program.
     *              It is a fatal error not to be able to do the dynamic cast.
     *
     *    @param  Return flag. Returns a zero if everything is ok. Anything else is a fatal error.
     */
    virtual int setInitialConditions(ELECTRODE_KEY_INPUT* ei);

    //! Initialize the sizes
    void init_sizes();

    //! Initialize the distribution of species evenly across the radius
    /*!
     *  This routine takes the spMoles values and spreads the values evenly in the radial
     *  direction so that the radial concencentrations are constant
     */
    void initializeAsEvenDistribution();


    //! Get all of the reaction rates and parameters from Cantera
    void extractInfoJustBorn(std::vector<int>& justBornMultiSpecies);

    int calcResid(double* const resid,  const ResidEval_Type_Enum evalType) ;


    //! Set the internal initial intermediate and initial global state from the internal final state
    /*!
     *  (virtual function)
     *
     *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
     *  routine as well.
     *
     * @param setInitInit   Boolean indicating whether you should set the init_init state as well
     */
    virtual void setInitStateFromFinal(bool setInitInit = false);

    //! Set the internal final intermediate state from the internal init state
    /*!
     *  (virtual function from Electrode)
     *
     *  Set the final state from the init state. This is commonly called during a failed time step
     */
    virtual void setFinalStateFromInit();

    //! Set the internal initial intermediate from the internal initial global state
    /*!
     *  Set the intial state from the init init state. We also can set the final state from this
     *  routine as well.
     *
     *  The final_final is not touched.
     *
     * @param setFinal   Boolean indicating whether you should set the final as well
     */
    virtual void setInitStateFromInitInit(bool setFinal = false);


    //! Set the internal initial intermediate and initial global state from the internal final state
    /*!
     *  (virtual function from Electrode)
     *
     *  Set the intial state from the final state. We also can set the init_init state from this
     *  routine as well.
     *
     * @param setInitInit   Boolean indicating whether you should set the init_init state as well. When
     *                      we do this we set the final state as well.
     */
    virtual void setInitStateFromFinalFinal(bool setInitInit = false);

    //! Set the internal final global state from the internal final intermediate state
    /*!
     *  (virtual function from Electrode)
     *
     *  Set the final_final state from the final state. This is commonly called at the end of successful base integration
     */
    virtual void setFinalFinalStateFromFinal();

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
     *  This routine does not calculate the global summed state. It does add contributions
     *  into spMoles_final_ and that's it.
     *
     *  We take the unknowns from the nonlinear problem and calculate all the values
     *
     *  All of these properties are defined for the _final_ state.
     *
     *  Thus, this is the main routine that reconciles all of the state information within the object.
     *  At the end of this routine, all aspects of the final state are consistent with each other.
     *
     *  Fundamental Variables:
     *       concKRSpecies_Cell_final_[]
     *
     *  Quantitities filled by this routine
     *     Distributed Quantities:
     *        cellBoundR_final_[iCell]
     *        spMoles_KRsolid_Cell_final_
     *        spMf_KRsolid_Cell_final_
     *        concTot_SPhase_Cell_final_
     *        partialMolarVolKRSpecies_Cell_final_
     *     Global Quantities:
     *         spMoles_final_   (additions)
     *
     *
     *  @param zeroGlobals   Zero the globals before adding to spMoles.
     */
    void updateStateDistrib(bool zeroGlobals);

    //! Generate a grid of radial points
    /*!
     *  @param numPoints   Number of grid points
     *  @param radiusInner  Inner radius
     *  @param radiusOuter  Outer Radius
     *  @param geomFactor   Geometric variation of the radii
     *                          Values greater than one increase density at the outer
     *                          Values less than one increase density at the inner
     *  @param radialMesh   Returned mesh
     */
    void radialGridGenerate(int numPoints, doublereal radiusInner, doublereal radiusOuter,
                            doublereal geomFactor,
                            std::vector<double>& radialMesh) const;

protected:

    //! Define the number of species that are defined to have radially distributed distributions
    //! within this solid region
    /*!
     *   Note, for this object there is only one radial distribution. This concept will be enhanced
     *   in later formulations.
     *
     *   There are a few arrays which have numKRSpecies_ as their inner loop dimension.
     */
    int numKRSpecies_;

    //!  Number of cells involved with the radial distribution, including the 1/2 end cells
    int numRCells_;

    //! Number of phases which have radial distributions
    int numSPhases_;

    //! Total number of equations defined at each node of the radial mesh
    int numEqnsCell_;

    //! Phase indeces of the solid phases comprising the species that are radially distributed
    //! within the current region
    /*!
     *   This is a very important indexing array.
     *   There are numSPhase_ of these.
     *
     *  Length =  numSPhase_
     *  access method = phaseIndeciseKRsolidPhases_[iSPphase_] = iPhaseIndex;
     *
     *  The iPhaseIndex is the index within the Electrode object.
     */
    std::vector<int> phaseIndeciseKRsolidPhases_;

    //! Species Index for the solid species that are distributed radially through the particle
    /*!
     *  The main object is a vector of species representing the homogenized values of the mole
     *  numbers of all species and phases that contribute to the Electrode.
     *  This vector contains the mapping between this index and the species that are radially
     *  distributed throughout the particle.
     *
     *    KRsolid_speciesList_[iKRSpecies] = kElectrodeObject
     */
    std::vector<int> KRsolid_speciesList_;



    //! Vector of solid species moles defined on the grid of the spherical particle, local final state
    /*!
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies.
     *
     *     Length = numKRSpecies_ * numRCells_
     *     access method:  spMoles_KRsolid_Cell_final_[numKRSpecies_ * iCell + iKRSpecies]
     *     units = kmol m-3.
     */
    std::vector<double> spMoles_KRsolid_Cell_final_;

    //! Vector of solid species moles defined on the grid of the spherical particle, local init state
    /*!
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies.
     *
     *     Length = numKRSpecies_ * numRCells_
     *     access method:  spMoles_KRsolid_Cell_init_[numKRSpecies_ * iCell + iKRSpecies]
     *     units = kmol m-3.
     */
    std::vector<double> spMoles_KRsolid_Cell_init_;


    //! Vector of solid species moles defined on the grid of the spherical particle, global final state
    /*!
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies.
     *
     *     Length = numKRSpecies_ * numRCells_
     *     access method:  spMoles_KRsolid_Cell_final_final_[numKRSpecies_ * iCell + iKRSpecies]
     *     units = kmol m-3.
     */
    std::vector<double> spMoles_KRsolid_Cell_final_final_;

    //! Vector of solid species moles defined on the grid of the spherical particle, global init state
    /*!
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies.
     *
     *     Length = numKRSpecies_ * numRCells_
     *     access method:  spMoles_KRsolid_Cell_init_init_[numKRSpecies_ * iCell + iKRSpecies]
     *     units = kmol m-3.
     */
    std::vector<double> spMoles_KRsolid_Cell_init_init_;

    //! Vector of solid species mole fractions defined on the grid of the spherical particle, local final state
    /*!
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies.
     *
     *     Length = numKRSpecies_ * numRCells_
     *     access method:  spMoles_KRsolid_Cell_final_[numKRSpecies_ * iCell + iKRSpecies]
     *     units = none
     */
    std::vector<double> spMf_KRsolid_Cell_final_;

    //! Vector of solid species mole fractions defined on the grid of the spherical particle, local init state
    /*!
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies.
     *
     *     Length = numKRSpecies_ * numRCells_
     *     access method:  spMoles_KRsolid_Cell_init_[numKRSpecies_ * iCell + iKRSpecies]
     *     units = none
     */
    std::vector<double> spMf_KRsolid_Cell_init_;

    //! Vector of solid species mole fractions defined on the grid of the spherical particle, global final state
    /*!
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies.
     *
     *     Length = numKRSpecies_ * numRCells_
     *     access method:  spMoles_KRsolid_Cell_final_[numKRSpecies_ * iCell + iKRSpecies]
     *     units = none
     */
    std::vector<double> spMf_KRsolid_Cell_final_final_;

    //! Vector of solid species mole fractions defined on the grid of the spherical particle, global init state
    /*!
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies.
     *
     *     Length = numKRSpecies_ * numRCells_
     *     access method:  spMoles_KRsolid_Cell_init_[numKRSpecies_ * iCell + iKRSpecies]
     *     units = none
     */
    std::vector<double> spMf_KRsolid_Cell_init_init_;

    //! Total concentration of each of the solid phases that are distributed - local final state
    /*!
     *     Length = numSPhase_ * numRCells_
     *     access method:  concTot_SPhase_Cell_final_[numSPhase * iCell + iSPhase]
     *     units = kmol m-3.
     */
    std::vector<double> concTot_SPhase_Cell_final_;

    //! Total concentration of each of the solid phases that are distributed - local init state
    /*!
     *     Length = numSPhase_ * numRCells_
     *     access method:  concTot_SPhase_Cell_init_[numSPhase * iCell + iSPhase]
     *     units = kmol m-3.
     */
    std::vector<double> concTot_SPhase_Cell_init_;

    //! Total concentration of each of the solid phases that are distributed - global final state
    /*!
     *     Length = numSPhase_ * numRCells_
     *     access method:  concTot_SPhase_Cell_final_[numSPhase * iCell + iSPhase]
     *     units = kmol m-3.
     */
    std::vector<double> concTot_SPhase_Cell_final_final_;

    //! Total concentration of each of the solid phases that are distributed - global init state
    /*!
     *     Length = numSPhase_ * numRCells_
     *     access method:  concTot_SPhase_Cell_init_init_[numSPhase * iCell + iSPhase]
     *     units = kmol m-3.
     */
    std::vector<double> concTot_SPhase_Cell_init_init_;

    //! Total concentration of each of the solid phase species that are distributed - final local state
    /*!
     *     Length = numKRSpecies_ * numRCells_
     *     access method:  concTot_SPhase_Cell_init_[numKRSpecies_ * iCell + iKRSpecies]
     *     units = kmol m-3.
     */
    std::vector<double> concKRSpecies_Cell_final_;

    //! total concentration of the solid phase species that are distributed - init loacl state
    /*!
     *     Length = numKRSpecies_ * numRCells_
     *     access method:  concTot_SPhase_Cell_init_[numKRSpecies_ * iCell + iKRSpecies]
     *     units = kmol m-3.
     */
    std::vector<double> concKRSpecies_Cell_init_;

    //! total concentration of the solid phase species that are distributed - final global state
    /*!
     *     Length = numKRSpecies_ * numRCells_
     *     access method:  concTot_SPhase_Cell_init_[numKRSpecies_ * iCell + iKRSpecies]
     *     units = kmol m-3.
     */
    std::vector<double> concKRSpecies_Cell_final_final_;

    //! total concentration of the solid phase species that are distributed - init state of global
    /*!
     *     Length = numKRSpecies_ * numRCells_
     *     access method:  concTot_SPhase_Cell_init_[numKRSpecies_ * iCell + iKRSpecies]
     *     units = kmol m-3.
     */
    std::vector<double> concKRSpecies_Cell_init_init_;

    //! Molar Volume of the reference lattice - final state of local step
    /*!
     *      Length numRCells_
     *      access method: molarVolume_refLat_Cell_final_[iCell];
     *      units = m3 kmol-1
     */
    std::vector<doublereal> molarVolume_refLat_Cell_final_;

    //! Molar Volume of the reference lattice - init state of local step
    /*!
     *      Length numRCells_
     *      access method: molarVolume_refLat_Cell_init_[iCell];
     *      units = m3 kmol-1
     */
    std::vector<doublereal> molarVolume_refLat_Cell_init_;

    //! Molar Volume of the reference lattice - final state of global step
    /*!
     *      Length numRCells_
     *      access method: molarVolume_refLat_Cell_final_final_[iCell];
     *      units = m3 kmol-1
     */
    std::vector<doublereal> molarVolume_refLat_Cell_final_final_;

    //! Molar Volume of the reference lattice - init state of global step
    /*!
     *      Length numRCells_
     *      access method: molarVolume_refLat_Cell_init_init_[iCell];
     *      units = m3 kmol-1
     */
    std::vector<doublereal> molarVolume_refLat_Cell_init_init_;


    //! Molar density of the solid phase in each cell under reference conditions
    /*!
     *  This is the molar volume of the first species in the mechanism at the
     *  starting temperature and pressure.
     *
     *  units are kmol m-3.
     */
    doublereal MolarVolume_refLat_Ref_;


    //! Node position of the mesh - final
    /*!
     *  Value of the node associated with the cell
     *  Vector of length number of cells =  numRCells_
     *  Index is the cell number
     */
    std::vector<doublereal> rnodePos_final_;

    //! Node position of the mesh - init
    /*!
     *  Value of the node associated with the cell
     *  Vector of length number of cells =  numRCells_
     *  Index is the cell number
     */
    std::vector<doublereal> rnodePos_init_;

    //! Node position of the mesh - final_final
    /*!
     *  Value of the node associated with the cell
     *  Vector of length number of cells =  numRCells_
     *  Index is the cell number
     */
    std::vector<doublereal> rnodePos_final_final_;

    //! Node position of the mesh - init_init
    /*!
     *  Value of the node associated with the cell
     *  Vector of length number of cells =  numRCells_
     *  Index is the cell number
     */
    std::vector<doublereal> rnodePos_init_init_;

    // Change the following to cellBoundRRefPos
    //! Reference radius at the right cell boundary - local final value
    std::vector<doublereal> rRefPos_final_;

    //! Reference radius at the right cell boundary - local init value
    std::vector<doublereal> rRefPos_init_;


    //! Reference radius at the right cell boundary - global final value
    std::vector<doublereal> rRefPos_final_final_;

    //! Reference radius at the right cell boundary - global init value
    std::vector<doublereal> rRefPos_init_init_;


    //! Value of the Radius' of the right cell boundaries at the end of the local step
    /*!
     *  Vector of length number of cells =  numRCells_
     *  Index is the cell number
     *  units = meters (note this is usually a micron sized quantity)
     */
    std::vector<doublereal> cellBoundR_final_;

    //! Value of the Radius' of the right cell boundaries at the beginning of the local step
    /*!
     *  Vector of length number of cells =  numRCells_
     *  Index is the cell number
     *  units = meters (note this is usually a micron sized quantity)
     */
    std::vector<doublereal> cellBoundR_init_;

    //! Value of the Radius' of the right cell boundaries at the end of the global step
    /*!
     *  Vector of length number of cells =  numRCells_
     *  Index is the cell number
     *  units = meters (note this is usually a micron sized quantity)
     */
    std::vector<doublereal> cellBoundR_final_final_;

    //! Value of the Radius' of the right cell boundaries at the beginning of the global step
    /*!
     *  Vector of length number of cells =  numRCells_
     *  Index is the cell number
     *  units = meters (note this is usually a micron sized quantity)
     */
    std::vector<doublereal> cellBoundR_init_init_;


    //!  Partial molar volume of all of the solid species located in all of the cells at the end of the local step
    /*!
     *   Vector of molar volumes for all solid species in all cells (KRSpecies, iCell)
     *   noUnits are m3 / kmol
     */
    std::vector<doublereal> partialMolarVolKRSpecies_Cell_final_;

    //!  Partial molar volume of all of the solid species located in all of the cells at the beginning of the local step
    /*!
     *     Length = numKRSpecies_ * numRCells_
     *     access method:   partialMolarVolKRSpecies_Cell_final_[numKRSpecies_ * iCell + iKRSpecies]
     *   Vector of molar volumes for all solid species in all cells (KRSpecies, iCell)
     *   Units are m3 / kmol
     */
    std::vector<doublereal> partialMolarVolKRSpecies_Cell_init_;

    //!  Partial molar volume of all of the solid species located in all of the cells at the end of the global step
    /*!
     *   Vector of molar volumes for all solid species in all cells (KRSpecies, iCell)
     *   noUnits are m3 / kmol
     */
    std::vector<doublereal> partialMolarVolKRSpecies_Cell_final_final_;

    //!  Partial molar volume of all of the solid species located in all of the cells at the beginning of the global step
    /*!
     *     Length = numKRSpecies_ * numRCells_
     *     access method:   partialMolarVolKRSpecies_Cell_final_[numKRSpecies_ * iCell + iKRSpecies]
     *   Vector of molar volumes for all solid species in all cells (KRSpecies, iCell)
     *   Units are m3 / kmol
     */
    std::vector<doublereal> partialMolarVolKRSpecies_Cell_init_init_;


    //!  Spline system for the nodal equations
    /*!
     *   These factors are the fraction of the total lattice growth or shrinkage that the
     *   current cell will take from the LHS of the domain. This does not change as a function of time.
     *   Therefore, the RHS will have a value of one, and this function will be monotonically
     *   increasing from left to right.
     */
    std::vector<doublereal> fractionVolExpansion_Cell_;

    //!  Molar creation rate of species in the electrode object due to the Exterior right surface
    /*!
     *   Length =  m_NumTotSpecies;
     *   access method:    DspMoles_RightSurf_final_[iSpecies]
     *   Vector of creation rates of species in the electrode object
     *
     *    units (kmol sec-1);
     */
    std::vector<doublereal> DspMoles_RightSurf_final_;


    //! Radius of the left boundary - reference value.
    doublereal radiusLeft_Ref_;

    //! Radius of the left boundary - value at the end of the local step.
    doublereal radiusLeft_final_;

    //! Radius of the left boundary - value at the beginning of the local step.
    doublereal radiusLeft_init_;

    //! Radius of the left boundary - value at the end of the local step.
    doublereal radiusLeft_final_final_;

    //! Radius of the left boundary - value at the beginning of the local step.
    doublereal radiusLeft_init_init_;

    //! Diffusion coefficient of species that are distributed
    /*!
     *   Here we are assuming that they are constant
     *
     *         Length = numKRSpecies_
     *         units = m2/s
     */
    std::vector<doublereal> Diff_Coeff_KRSolid_;

    //! Production Rate of phases at the right boundary
    /*!
     *   Length =  m_NumTotPhases
     *   access method:    DphMoles_RightSurf_final_[iPhase]
     *   Vector of creation rates of phase in the electrode object
     *
     *    units (kmol sec-1);
     */
    std::vector<doublereal> DphMoles_RightSurf_final_;

    //! we identify the phases here as being the right surface, which may be the exterior surface
    /*!
     *  The phase will be in direct contact with the electrolyte or it may not be
     */
    int surfIndex_RightSurf_;

    //! Vector of activity coefficients for all KR species at all nodes
    /*!
     *    Length = numKRSpecies_
     *    access method:   actCoeff_[iKRSpecies]
     *    This is dimensionless , intermediate number that is calculated at every point.
     */
    std::vector<doublereal> actCoeff_Cell_final_;


     Electrode_RadialDiffRegions* ee_;
};

}


#endif
/*****************************************************************************/
