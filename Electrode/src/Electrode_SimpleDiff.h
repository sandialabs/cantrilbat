/*
 * $Id: Electrode_SimpleDiff.h 571 2013-03-26 16:44:21Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_SIMPLEDIFF_H
#define _ELECTRODE_SIMPLEDIFF_H



#include "Electrode.h"
#include "Electrode_Integrator.h"
#include "cantera/integrators.h"
#include "cantera/numerics/ResidJacEval.h"

namespace Cantera
{

#define ELECTRODETYPE_ANODE   0
#define ELECTRODETYPE_CATHODE 1


//! This class is a derived class used to model phase - change electrodes
/*!
 * Complete problem statement
 *
 */
class Electrode_SimpleDiff : public Electrode_Integrator
{
public:
    //! Constructor
    Electrode_SimpleDiff();

    //! Destructor
    virtual  ~Electrode_SimpleDiff();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    Electrode_SimpleDiff(const Electrode_SimpleDiff& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    Electrode_SimpleDiff& operator=(const Electrode_SimpleDiff& right);

    //! Return the type of electrode
    /*!
     *  Returns the enum type of the electrode. This is used in the factory routine.
     *
     *  @return Returns an enum type, called   Electrode_Types_Enum
     */
    virtual Electrode_Types_Enum electrodeType() const;

  
    //! Add additional Keylines for child electrode objects, and then read them in
    /*!
     *   (virtual function from Electrode)
     *   (overload virtual function - replaces parent members)
     *
     *   This function will replace the ELECTRODE_KEY_INPUT structure with an expanded
     *   child member structure ELECTRODE_RadialRegion_KEY_INPUT containing the extra information.
     *
     *   If the command file has been read before, it will then reparse the command file
     *   storring the new information in the  ELECTRODE_RadialRegion_KEY_INPUT structure.
     *
     *    @param ei_ptr  Handle to the ELECTRODE_KEY_INPUT base pointer. This handle may change
     *                   as the child class of  ELECTRODE_KEY_INPUT gets malloced.
     *
     *    @return  0 successful but no change in ei
     *             1 Successful and ei has changed
     *            -1 unsuccessful fatal error of some kind.
     */
    virtual int electrode_input_child(ELECTRODE_KEY_INPUT** ei_ptr);
    

    //!  Setup the electrode for first time use
    /*!
     * (virtual from Electrode  - onion Out)
     *
     *    This is one of the most important routines. It sets up the electrode's internal structures
     *    After the call to this routine, the electrode should be internally ready to be integrated
     *    and reacted.
     *    It takes its input from an ELECTRODE_KEY_INPUT object which specifies the setup of the electrode
     *    object and the initial state of that object.
     *    The routine works like an onion Out initialization. The parent object is initialized before the
     *    child. This means the child object first calls the parent, before it does its own initializations.
     *
     *    There are some virtual member functions that won't work until this routine is called. That's because
     *    the data structures won't be set up for base and child Electrode objects until this is called.
     *
     *  @param ei   BASE ELECTRODE_KEY_INPUT pointer object. Note, it must have the correct child class
     *              for the child electrode object.
     *
     *  @return  Returns zero if successful, and -1 if not successful.
     */
    int electrode_model_create(ELECTRODE_KEY_INPUT* ei);

    //!  Set the electrode initial conditions from the input file.
    /*!
     *   (virtual from Electrode)
     *   (This is a serial virtual function or an overload function)
     *
     *    This is one of the most important routines. It sets up the initial conditions of the electrode
     *    from the input file. The electrode itself has been set up from a call to electrode_model_create().
     *    After the call to this routine, the electrode should be internally ready to be integrated and reacted.
     *    It takes its input from an ELECTRODE_KEY_INPUT object which specifies the setup of the electrode
     *    object and the initial state of that object.
     *
     *    The routine works like an onion initialization. The parent object is initialized before the
     *    child. This means the child object first calls the parent, before it does its own initializations.
     *
     * @param ei    ELECTRODE_KEY_INPUT pointer object
     *
     *  @return  Returns zero if successful, and -1 if not successful.
     */
    virtual int setInitialConditions(ELECTRODE_KEY_INPUT* ei);


    //! Initialize the sizes
    void init_sizes();

    //! Initilize the grid
    void init_grid();

    //! Initialize the distribution of species evenly across the radius
    /*!
     *  This routine takes the spMoles values and spreads the values evenly in the radial
     *  direction so that the radial concencentrations are constant
     */
    void initializeAsEvenDistribution();

 

    // --------------------------------------  QUERY VOLUMES -----------------------------------------------------------
    //!    Return the total volume of solid material
    /*!
     *  (virtual from Electrode.h)
     *
     *       This is the main routine for calculating the
     *       We leave out the solnPhase_ volume from the calculation
     *       units = m**3
     */
    virtual double SolidVol() const;

    //!  Returns the total moles in the electrode phases of the electrode
    /*!
     *  @return total moles (kmol)
     */
    virtual double SolidTotalMoles() const;

    // --------------------------------------------- SURFACE AREAS -------------------------------------------------------



    void calcRate(double deltaT);




    //!  Print one set of field variables
    /*!
     *  @param indentSpaces  Number of spacies to indent the whole group
     *  @param 
     */
    void showOneField(const std::string &title, int indentSpaces, const double * const radialValues, int numRadialVals, 
		      const double * const vals, const std::vector<std::string> &varNames, int numFields);

    void  showOneFieldInitFinal(const std::string &title, int indentSpaces, const double * const radialValues, int numRadialVals, 
				const double * const vals_init,  const double * const vals_final,
				const std::vector<std::string> &varNames, int numFields);

    void showOneResid(const std::string &title, int indentSpaces, const double * const radialValues, int numRadialVals, 
		      int numFields1, int iTerm, const double * const val_init,  const double * const val_final,
		      int numEqnsCell, int iEqn,const double * const resid_error,  const double * const solnError_tol,
		      const double * const residual);

    //! Print the solution to the current step to output
    void showSolution(int indentSpaces);
    void showResidual(int indentSpaces,  const double * const residual);



    //! Extract information from reaction mechanisms
    /*!
     *   (inherited from Electrode_Integrator)
     */
    virtual void extractInfo();

    //! Collect mole change information
    /*!
     *   We take the ROP_inner_[] and ROP_outer_[] rates of progress, combine them with the surface area calculation,
     *   and the stoichiometric coefficients to calculate the DspMoles_final_[], which is production rate for
     *   all species in the electrode due to surface reactions
     *
     *   (inherited from Electrode_Integrator)
     */
    virtual void updateSpeciesMoleChangeFinal();


    //------------------------------------------------------------------------------------------------------------------
    // -------------------------------  SetState Functions -------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------

    //! Set the internal initial intermediate and initial global state from the internal final state
    /*!
     *  (virtual function from Electrode.h)
     *
     *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
     *  routine as well.
     *
     * @param setInitInit   Boolean indicating whether you should set the init_init state as well
     */
    virtual void setInitStateFromFinal(bool setInitInit = false);


    //! Set the internal initial intermediate and initial global state from the internal final_final state
    /*!
     *  (virtual function from Electrude.h)
     *
     *  Set the init and init_init state from the final_final state.
     */
    virtual void setInitInitStateFromFinalFinal();


    //! Set the internal final intermediate state from the internal init state
    /*!
     *  (virtual function from Electrode)
     *
     *  Set the final state from the init state. This is commonly called during a failed time step
     */
    virtual void setFinalStateFromInit();

    //! Set the internal initial intermediate from the internal initial global state
    /*!
     *  (virtual function from Electrode.h)
     *
     *  Set the intial state from the init init state. We also can set the final state from this
     *  routine as well.
     *
     *  The final_final is not touched.
     *
     * @param setFinal   Boolean indicating whether you should set the final as well
     */
    virtual void setInitStateFromInitInit(bool setFinal = false);

    //! Set the internal final global state from the internal final intermediate state
    /*!
     *  (virtual function from Electrode)
     *
     *  Set the final_final state from the final state. This is commonly called at the end of successful base integration
     */
    virtual void setFinalFinalStateFromFinal();

    //--------------------------------------------------------------------------------------------------
    // -----------------------  STATE and PRINTING FUNCTIONS ----------------------------------------
    //--------------------------------------------------------------------------------------------------

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

    virtual void updatePhaseNumbers(int);

    //! Take the state (i.e., the final state) within the Electrode_SimpleDiff and push it up
    //! to the zero-dimensional parent object
    /*!
     * 
     *  update:
     *             spMoles_final_ [] -> sum solid phase species
     *             spMf_final_[]  -> Use exterior cell values
     *
     */
    void updateState_OneToZeroDimensions();

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

    //! Predict the solution
    /*!
     * Ok at this point we have a time step deltalimiTsubcycle_
     * and initial conditions consisting of phaseMoles_init_ and spMF_init_.
     * We now calculate predicted solution components from these conditions.
     *
     * @return   Returns the success of the operation
     *                 1  A predicted solution is achieved
     *                 2  A predicted solution with a multispecies phase pop is acheived
     *                 0  A predicted solution is not achieved, but go ahead anyway
     *                -1  The predictor suggests that the time step be reduced and a retry occur.
     */
    int predictSolnResid();

    //! Predict the solution
    /*!
     * Ok at this point we have a time step deltalimiTsubcycle_
     * and initial conditions consisting of phaseMoles_init_ and spMF_init_.
     * We now calculate predicted solution components from these conditions.
     *
     *  Calls predictSolnResid(). If it works that good, if not reduct time step and call again.
     *
     * @return   Returns the success of the operation
     *                 1  A predicted solution is achieved
     *                 2  A predicted solution with a multispecies phase pop is acheived
     *                 0  A predicted solution is not achieved, but go ahead anyway
     *                -1  The predictor suggests that the time step be reduced and a retry occur.
     */
    virtual int predictSoln();

    //! Unpack the soln vector
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *  This function unpacks the solution vector into  phaseMoles_final_,  spMoles_final_, and spMf_final_[]
     */
    virtual void unpackNonlinSolnVector(const double* const y);

    //! Set the base tolerances for the nonlinear solver within the integrator
    /*!
     *   The tolerances are based on controlling the integrated electron source term
     *   for the electrode over the integration interval.  The integrated source term
     *   has units of kmol.
     *
     *   Because the electron is only one molar quantity within a bunch of molar quantities,
     *   this requirement will entail that we control the source terms of all species within the
     *   electrode to the tolerance requirements of the electron source term.
     *
     *   @param rtolResid  Relative tolerance allowed for the electron source term over the interval.
     *                     This is a unitless quantity
     */
    virtual void setNLSGlobalSrcTermTolerances(double rtolResid);

    //!   Set the Residual absolute error tolerances
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *   Set the absolute error tolerances for the residuals for the nonlinear solvers. This is called at the top
     *   of the integrator() routine.
     *
     *   Calculates residAtolNLS_[]
     *   Calculates atolNLS_[]
     */
    void setResidAtolNLS();

    //! Evaluate the residual function
    /*!
     * (virtual from NonlinearSolver)
     *
     * @param t             Time                    (input)
     * @param delta_t       The current value of the time step (input)
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input, do not modify)
     * @param resid         Value of the residual that is computed (output)
     * @param evalType      Type of the residual being computed (defaults to Base_ResidEval)
     * @param id_x          Index of the variable that is being numerically differenced to find
     *                      the jacobian (defaults to -1, which indicates that no variable is being
     *                      differenced or that the residual doesn't take this issue into account)
     * @param delta_x       Value of the delta used in the numerical differencing
     *
     * @return
     */
    virtual int evalResidNJ(const doublereal t, const doublereal delta_t,
                            const doublereal* const y,
                            const doublereal* const ydot,
                            doublereal* const resid,
                            const ResidEval_Type_Enum evalType = Base_ResidEval,
                            const int id_x = -1,
                            const doublereal delta_x = 0.0);

    //!  Residual calculation for the solution of the Nonlinear integration problem
    /*!
     * @param t             Time                    (input)
     * @param delta_t       The current value of the time step (input)
     * @param y             Solution vector (input, do not modify)
     * @param ySolnDot      Rate of change of solution vector. (input, do not modify)
     * @param resid         Value of the residual that is computed (output)
     * @param evalType      Type of the residual being computed (defaults to Base_ResidEval)
     * @param id_x          Index of the variable that is being numerically differenced to find
     *                      the jacobian (defaults to -1, which indicates that no variable is being
     *                      differenced or that the residual doesn't take this issue into account)
     * @param delta_x       Value of the delta used in the numerical differencing
     */
    int integrateResid(const doublereal t, const doublereal delta_t,
		       const doublereal* const y, const doublereal* const ySolnDot,
		       doublereal* const resid,
		       const ResidEval_Type_Enum evalType, const int id_x,
		       const doublereal delta_x);

    //! Main internal routine to calculate the residual
    /*!
     *  This routine calculates the functional at the current stepsize, deltaTsubcycle_.
     *  A new stepsize, deltaTsubcycleCalc_, is calculated within this routine for changes
     *  in topology of type LimitingEventDeltaTsubcycle_
     *
     *  This routine calcules yval_retn, which is the calculated value of the residual for the
     *  nonlinear function
     *
     *   resid[i] = y[i] - yval_retn[i]
     *
     *  The formulation of the solution vector is as follows. The solution vector will consist of the following form
     *
     *     y =   phaseMoles_final[iph]    for iph = phaseIndexSolidPhases[0]
     *           phaseMoles_final[iph]    for iph = phaseIndexSolidPhases[1]
     *           . . .
     *           Xmol[k = 0]              for iph = phaseIndexSolidPhases[0]
     *           Xmol[k = nSpecies()-1]   for iph = phaseIndexSolidPhases[0]
     *           Xmol[k = 0]              for iph = phaseIndexSolidPhases[0]
     *           Xmol[k = nSpecies()-1]   for iph = phaseIndexSolidPhases[0]
     *           ...
     *           Xmol[k = 0]              for iph = phaseIndexSolidPhases[1]
     *           Xmol[k = nSpecies()-1]   for iph = phaseIndexSolidPhases[1]

     *
     *                                                             Unknown                           Index
     * --------------------------------------------------------------------------------------------------------------
     *         Residual (Time)                                     deltaSubcycleCalc_                   0
     *                                                                                            1
     *         Loop over cells                                                            0 <=  iCell < numRCells_
     *                                                                                     j = numEqnsCell_ * iCell
     *                    
     *            Residual (Reference/lattice Position)            rRefPos_final_[iCell];        (1+j)
     *            Residual (Mesh Position)                         rnodePos_final_[iCell]        (j+1) + 1
     *            Loop over distributed Phases
     *            Residual (Concentration _ k=0)                  concTot_SPhase_Cell_final_[iCell * numSPhases_ + jPh]
     *              . . .
     *            Residual (Concentration _ k=Ns-1)               concKRSpecies_Cell_final_[iCell * numKRSpecies_ + iKRSpecies]
     *  --------------------------------------------------------------------------------------------------------------
     *
     *  @param resid   Calculated residual vector whose form is described above
     */
    int  calcResid(double* const resid, const ResidEval_Type_Enum evalType);



protected:

    //! Type of the electrode, 0 for anode, 1 for cathode
    /*!
     *  0 - anode
     *  1 - cathode
     *
     *  The open circuit plateaus are ramped differently.
     *   For anodes, the open circuit plateaus are ramped going upwards
     *       A big negative voltage means that the DoD will stay or trend to 0
     *       A big positive voltage means that the DoD will stay or trend to 1.
     *
     *   For cathodes, the open circuit plateaus are ramped going downwards.
     *       A big positive voltage means that the DoD will stay or trend to 0
     *       A big negative voltage means that the DoD will stay or trend to 1.
     */
    int electrodeType_;

    //! Define the number of species that are defined to have radially distributed distributions
    //! within the solid
    /*!
     *   Note, for this object there is only one radial distribution. This concept will be enhanced
     *   in later formulations.
     *
     *   There are a few arrays which have numKRSpecies_ as their inner loop dimension.
     */
    int numKRSpecies_;

    //!  Number of cells involved with the radial distribution, including the 1/2 end cells
    int numRCells_;

    //! Number of phases which have radial distributions of their species
    int numSPhases_;

    //! Total number of equations defined at each node of the radial mesh
    int numEqnsCell_;

    std::vector<ThermoPhase *> thermoSPhase_List_;

    //! Vector of solid species defined on the grid of the spherical particle
    /*!
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies;
     *   final state value
     */
    std::vector<double> spMoles_KRsolid_Cell_final_;

    //! Vector of solid species defined on the grid of the spherical particle
    /*!
     *
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies;
     *   init state value
     */
    std::vector<double> spMoles_KRsolid_Cell_init_;

    //! Vector of solid species defined on the grid of the spherical particle
    /*!
     *
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies;
     *   final_final state value
     */
    std::vector<double> spMoles_KRsolid_Cell_final_final_;

    //! Vector of solid species defined on the grid of the spherical particle
    /*!
     *
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies;
     *   init_init state value
     */
    std::vector<double> spMoles_KRsolid_Cell_init_init_;

    //! Species Index for the solid species that are distributed radially through the particle
    /*!
     *  The main object is a vector of species representing the homogenized values of the mole
     *  numbers of all species and phases that contribute to the Electrode.
     *  This vector contains the mapping between this index and the species that are radially
     *  distributed throughout the particle.
     *
     *    KRsolid_speciesList_[KRsolid] = kElectrodeObject
     *
     *      where
     *           KRsolid = index of the radially distributed species
     *           kElectrodeObject = PhaseList index of the species
     */
    std::vector<int> KRsolid_speciesList_;

    //! Species Names for the solid species that are distributed radially through the particle
    std::vector<std::string> KRsolid_speciesNames_;

    //! Phase indeces of the solid phases comprising the species that are radially distributed
    /*!
     *  There are numSPhases_ of these
     *
     *  The phaseIndex is the index within the Electrode object
     *
     *  
     *  iPh = phaseIndeciseKRsolidPhases_[distribPhIndex];
     *
     *
     */
    std::vector<int> phaseIndeciseKRsolidPhases_;

    //! Returns the distributed phase index given the regular phase Index.
    /*
     *  
     *  distribPhIndex = phaseIndeciseKRsolidPhases_[iPh];
     *
     *  If distribPhIndex is -1, the phase isn't distributed.
     */
    std::vector<int> distribPhIndexKRsolidPhases_;

    //! Number of species in each of the radially distributed phases
    /*!
     *  There are numSPhases_ of these
     *
     *  The phaseIndex is the index within the Electrode object
     */
    std::vector<int> numSpeciesInKRSolidPhases_;

    std::vector<int> kstartKRSolidPhases_;

    //! Phase Names of the solid phases comprising the species that are radially distributed
    /*!
     *  There are numSPhases_ of these
     */
    std::vector<std::string>  KRsolid_phaseNames_;

    //! Phase indeces of the solid phases comprising the species that are not radially distributed
    /*!
     *  There are numNonSPhase_ of these
     *
     *  The phaseIndex is the index within the Electrode object
     */
    std::vector<int> phaseIndeciseNonKRsolidPhases_;

    int numNonSPhases_;

    //! Total concentration of each of the solid phases that are distributed - global final state
    /*!
     *       concTot_SPhase_Cell_final_(iSPhase, iCell)
     */
    std::vector<double> concTot_SPhase_Cell_final_final_;

    //! Total concentration of each of the solid phases that are distributed - local final state
    /*!
     *       concTot_SPhase_Cell_final_(iSPhase, iCell)
     */
    std::vector<double> concTot_SPhase_Cell_final_;

    //! Total concentration of each of the solid phases that are distributed - local init state
    /*!
     *       concTot_SPhase_Cell_final_(iSPhase, iCell)
     */
    std::vector<double> concTot_SPhase_Cell_init_;

    //! Total concentration of each of the solid phases that are distributed - global init state
    /*!
     *       concTot_SPhase_Cell_final_(iSPhase, iCell)
     */
    std::vector<double> concTot_SPhase_Cell_init_init_;

    //! total concentration of the solid phases that are distributed - init state
    /*!
     *      concKRSpecies_Cell_init_[numKRSpecies_ * iCell + iKRSpecies]
     */
    std::vector<double> concKRSpecies_Cell_init_;

    //! total concentration of the solid phases that are distributed - init state
    /*!
     *     concKRSpecies_Cell_final_[numKRSpecies_ * iCell + iKRSpecies]
     */
    std::vector<double> concKRSpecies_Cell_final_;

    //! total concentration of the solid phases that are distributed - init state
    /*!
     *     concKRSpecies_Cell_init_init_[numKRSpecies_ * iCell + iKRSpecies]
     */
    std::vector<double> concKRSpecies_Cell_init_init_;

    //! total concentration of the solid phases that are distributed - init state
    /*!
     *     concKRSpecies_Cell_final_final_[numKRSpecies_ * iCell + iKRSpecies]
     */
    std::vector<double> concKRSpecies_Cell_final_final_;

    //! Mole fraction of the solid phase species that are distributed = final state
    /*!
     *
     */
    std::vector<double> spMf_KRSpecies_Cell_final_;

    //! Mole fraction of the solid phase species that are distributed = final state
    /*!
     *
     */
    std::vector<double> spMf_KRSpecies_Cell_init_;

    //! Molar density of the solid phase in each cell under reference conditions
    /*!
     *  This is the molar volume of the first species in the mechanism at the
     *  starting temperature and pressure.
     *
     *  units are kmol m-3.
     */
    doublereal MolarVolume_Ref_;

    //! Node position of the mesh - final_final
    std::vector<doublereal> rnodePos_final_final_;

    //! Node position of the mesh - final
    std::vector<doublereal> rnodePos_final_;

    //! Node position of the mesh - init
    std::vector<doublereal> rnodePos_init_;

    //! Node position of the mesh - init_init
    std::vector<doublereal> rnodePos_init_init_;

    //! Reference radius at the right cell boundary - global final value
    std::vector<doublereal> rLatticeCBR_final_final_;

    //! Reference radius at the right cell boundary - local final value
    std::vector<doublereal> rLatticeCBR_final_;

    //! Reference radius at the right cell boundary - local init value
    std::vector<doublereal> rLatticeCBR_init_;

    //! Reference radius at the right cell boundary - global init value
    std::vector<doublereal> rLatticeCBR_init_init_;

    std::vector<doublereal> rLatticeCBR_ref_;

    //! Lattice velocity of the right cell boundary - global final value
    std::vector<doublereal> vLatticeCBR_cell_;

    std::vector<doublereal> cellBoundR_final_;
    std::vector<doublereal> cellBoundL_final_;

    //! Volume of each cell on a per particle basis
    /*!
     *      Length numRCells:
     *      units = m3
     */
    std::vector<doublereal> volPP_Cell_final_;

    //!  Spline system for the nodal equations
    /*!
     *   These factors are the fraction of the exterior node radius that the
     *   current node possesses.
     */
    // std::vector<doublereal> fracNodePos_;

    //!  Spline System for the nodal points
    /*!
     *  Fraction of the domain's volume which is inside the current node position
     */
    std::vector<doublereal> fracVolNodePos_;

    //!  Partial molar volume of all of the solid species located in all of the cells
    /*!
     *   Vector of molar volumes for all solid species in all cells (KRSpecies, iCell)
     *   noUnits are m3 / kmol
     */
    std::vector<doublereal> partialMolarVolKRSpecies_Cell_final_;


    //! Rate of progress of the surface reactions
    /*!
     *   length = maxNumRxns;
     */
    std::vector<double> ROP_;


    //!  Molar creation rate of species in the electrode object due to the Exterior surface
    /*!
     *   This is a quantity over all species in the PhaseList
     *
     *    units (kmol sec-1);
     */
    std::vector<doublereal> DspMoles_final_;

    //! Domain boundary at the inner radius.
    /*!
     *   Frequently this will be zero. default is zero.
     */
    doublereal m_rbot0_;

    std::vector<doublereal> Diff_Coeff_KRSolid_;

    //! Molar creation rate of phases in the electrode object.
    /*!
     *  This is calculated as a sum of species creation rates summed
     *  up over all cells in the domain
     *
     *    units = kmol / s
     */
    std::vector<doublereal> DphMolesSrc_final_;

    //! we identify the phases here as being the exterior surface
    /*!
     *  The  phase will be in direct contact with the electrolyte.
     */
    int surfIndexExteriorSurface_;

    //!  Value of the total flux at the outer edge - kmol m-2 s-1
    doublereal NTflux_final_;


    //! Local value of the diffusion coefficient
    doublereal DiffCoeff_;

    //! Local value of the diffusion coefficient
    doublereal DiffCoeff_default_;

    //! Vector of activity coefficients for all KR species at all nodes
    /*!
     *    This is calculated at the final state
     */
    std::vector<doublereal> actCoeff_Cell_final_;


    //! phase id of the phase which will die at the shortest time
    int phaseID_TimeDeathMin_;

    //! Cell id of the shortest time to a phase death
    int cellID_TimeDeathMin_;

   //! This integer describes if the system is current on a Region boundary at the start of a subgrid
    //! integration step
    /*!  
     *  We define a region boundary here iff all cells are on that boundary
     *  If it isn't, we set it to -1. If it is, we set it to the region boundary index
     */
    int onRegionBoundary_init_;

    //! This integer describes if the system is current on a Region boundary at the end of a subgrid
    //! integration step
    /*!
     *   We define a region boundary here iff all cells are on that boundary
     *  If it isn't, we set it to -1. If it is, we set it to the region boundary index
     */
    int onRegionBoundary_final_;

    //! This integer describes if the system is current on a Region boundary at the start of a global
    //! integration step
    /*!
     *  If it isn't, we set it to -1. If it is, we set it to the region boundary index
     */
    int onRegionBoundary_init_init_;

    //! This integer describes if the system is current on a Region boundary at the end of a global
    //! integration step
    /*!
     *  If it isn't, we set it to -1. If it is, we set it to the region boundary index
     */
    int onRegionBoundary_final_final_;

    //! Absolute tolerance for nonlinear residual
    double atolBaseResid_;

    //!  This Boolean is true when we are at a plateau boundary and the voltage is in between the
    //!  top value and the bottom value.
    int goNowhere_;



};

}


#endif
/*****************************************************************************/
