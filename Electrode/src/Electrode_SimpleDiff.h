/**
 *  @file Electrode_SimpleDiff.h
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_SIMPLEDIFF_H
#define _ELECTRODE_SIMPLEDIFF_H

#include "Electrode_Integrator.h"

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{

//! anode type of electrode
#define ELECTRODETYPE_ANODE   0
//! cathode type of electrode
#define ELECTRODETYPE_CATHODE 1

class EState_RadialDistrib;

//==================================================================================================================================
//! This class is a derived class used to model intercalating electrodes where solid phase diffusion of ionic or solid-phase
//! intermediate species is a rate limiting or at least a slow step
/*!
 *  This is a distributed radial model with a simple diffusion mechanism, with linkages to a simple surface reaction mechanism
 *  at the exterior of the radial particle.
 *
 *  It is assumed throughout the model that there is one and only one ReactingSurDomain object and that it is located
 *  at the exterior of the particle, between the particle and the electrolyte.
 *
 *  The particle can swell due to reaction or diffusion.
 *
 *  Various diffusion mechanisms are allowed.
 *
 */
class Electrode_SimpleDiff : public Electrode_Integrator
{
public:
    //! Constructor
    Electrode_SimpleDiff();

    //! Destructor
    virtual ~Electrode_SimpleDiff();

    //! Copy Constructor
    /*!
     *  @param right Object to be copied
     */
    Electrode_SimpleDiff(const Electrode_SimpleDiff& right);

    //! Assignment operator
    /*!
     *  @param[in]           right               object to be copied
     *
     *  @return                                  Returns a reference to the current object
     */
    Electrode_SimpleDiff& operator=(const Electrode_SimpleDiff& right);

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
  
    //! Add additional Keylines for child electrode objects, and then read them in
    /*!
     *  (virtual function from Electrode)
     *  (overload virtual function - replaces parent members)
     *
     *  This function will replace the ELECTRODE_KEY_INPUT structure with an expanded
     *  child member structure ELECTRODE_RadialRegion_KEY_INPUT containing the extra information.
     *
     *  If the command file has been read before, it will then reparse the command file
     *  storring the new information in the  ELECTRODE_RadialRegion_KEY_INPUT structure.
     *
     *  @param[in]           ei_ptr              Handle to the ELECTRODE_KEY_INPUT base pointer. This handle may change
     *                                           as the child class of ELECTRODE_KEY_INPUT gets malloced.
     *
     *  @return                                   0 successful but no change in ei
     *                                            1 Successful and ei has changed
     *                                           -1 unsuccessful fatal error of some kind.
     */
    virtual int electrode_input_child(ELECTRODE_KEY_INPUT** ei_ptr) override;

    //! Setup the electrode for first time use
    /*!
     *  (virtual from Electrode  - onion Out)
     * 
     *  This is one of the most important routines. It sets up the electrode's internal structures
     *  After the call to this routine, the electrode should be internally ready to be integrated and reacted.
     *  It takes its input from an ELECTRODE_KEY_INPUT object which specifies the setup of the electrode
     *  object and the initial state of that object.
     *  The routine works like an onion Out initialization. The parent object is initialized before the
     *  child. This means the child object first calls the parent, before it does its own initializations.
     *
     *  There are some virtual member functions that won't work until this routine is called. That's because
     *  the data structures won't be set up for base and child Electrode objects until this is called.
     *
     *  @param[in]           ei                  BASE ELECTRODE_KEY_INPUT pointer object. Note, it must have the correct child class
     *                                           for the child electrode object.
     *
     *  @return                                  Returns zero if successful, and -1 if not successful.
     */
    virtual int electrode_model_create(ELECTRODE_KEY_INPUT* ei) override;

    //! Set the electrode initial conditions from the input file.
    /*!
     *  (virtual from Electrode)
     *  (This is a serial virtual function or an overload function)
     *
     *  This is one of the most important routines. It sets up the initial conditions of the electrode
     *  from the input file. The electrode itself has been set up from a call to electrode_model_create().
     *  After the call to this routine, the electrode should be internally ready to be integrated and reacted.
     *  It takes its input from an ELECTRODE_KEY_INPUT object which specifies the setup of the electrode
     *  object and the initial state of that object.
     * 
     *  The routine works like an onion initialization. The parent object is initialized before the
     *  child. This means the child object first calls the parent, before it does its own initializations.
     *
     *  @param[in]           ei                  ELECTRODE_KEY_INPUT pointer object
     *
     *  @return                                  Returns zero if successful, and -1 if not successful.
     */
    virtual int setInitialConditions(ELECTRODE_KEY_INPUT* ei) override;

    //! Create an object that saves the electrode state and can print out an XML solution to file
    /*!
     *  The pointer to the malloced object is saved in the internal variable eState_final_ .
     *  Because this pointer is now non-null, the state of the electrode will be saved at each step
     *  and the electrode object has a restart capability.
     *  If the pointer is null, no restart information is generated
     *
     *  @param[in]           force               Force the creation of a new eState object. Defaults to false.
     *  @return                                  Returns zero if successful, and -1 if not successful.
     */
    virtual int electrode_stateSave_create(bool force = false) override;

    //! Resize the solid phase and electrolyte mole numbers within the object
    /*!
     *  (virtual from Electrode)
     * 
     *  This routine uses particleDiameter_ , particleNumberToFollow_, and porosity_ to recalculate
     *  all the mole numbers in the electrode. This is done by rescaling all of the numbers.
     *  At the end of the process, the total volume of the electrode object is
     *
     *    grossVol = SolidVol() / ( 1.0 - porosity_)
     *
     *  where the SolidVol() is equal to
     *
     *   SolidVol() =  particleNumberToFollow_  Pi *  particleDiameter_**3 / 6.0;
     */
    virtual void resizeMoleNumbersToGeometry() override;

    //! Calculate the number of equations that will be solved during the nonlinear solver step.
    /*!
     *  (virtual from Electrode_Integrator)
     *  All classes which inherit from this routine must have a class that determines this value.
     *
     *  @return                                  Returns the number of unknowns in the nonlinear problem and time-stepping problem.
     */
    virtual size_t nEquations_calc() const override;

    //! Initialize the sizes
    void init_sizes();

    //! Initilize the grid
    void init_grid();

    //! Initialize the distribution of species evenly across the radius
    /*!
     *  This routine takes the spMoles values and spreads the values evenly in the radial
     *  direction so that the radial concencentrations and mole fractions are constant.
     */
    void initializeAsEvenDistribution();

    //-------------------------------------------------------------------------------------------------------------------
    // ----------------------------------- QUERY VOLUMES ----------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------
 
    //! Return the total volume of solid material (m**3)
    /*!
     *  (virtual from Electrode.h)
     *
     *  This is the main routine for calculating the solid volume of the electode.
     *  We leave out the solnPhase_ volume from the calculation.
     *
     *  @return                                  Returns the solid volume of the electrode (m**3)
     */
    virtual double SolidVol() const override;

    //-------------------------------------------------------------------------------------------------------------------
    // --------------------------------------  QUERY HEAT CAPACITY  -----------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------

    //! Returns the total Heat Capacity of the Material in the Solid Electrode at constant volume
    /*!
     *  (virtual from Electrode)
     *  This is an extensive quantity.
     *
     *  @return                                  Returns the heat capacity of the solid phases (Joule K-1)
     */
    virtual double SolidHeatCapacityCV() const override;

    //! Returns the total enthalpy of the Material in the Solid Electrode
    /*!
     *  (virtual from Electrode)
     *  This is an extensive quantity.
     *
     *  @return                                  Returns the total enthalpy (Joules)
     */
    virtual double SolidEnthalpy() const override;

    //-------------------------------------------------------------------------------------------------------------------
    // --------------------------------------------- SURFACE AREAS -------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------

    //! Extract information from reaction mechanisms
    /*!
     *  (virtual from Electrode_Integrator)
     *  Calculate the reaction rates on the one exterior surface
     */
    virtual void extractInfo() override;

    //! Collect mole change information
    /*!
     *  (virtual from Electrode_Integrator)
     * 
     *  We take the ROP_inner_[] and ROP_outer_[] rates of progress, combine them with the surface area calculation,
     *  and the stoichiometric coefficients to calculate the DspMoles_final_[], which is production rate for
     *  all species in the electrode due to surface reactions
     */
    virtual void updateSpeciesMoleChangeFinal() override;

    //-------------------------------------------------------------------------------------------------------------------
    // --------------------------------------- PRINTING UTILITIES -------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------

    //! Print a group of radially varying field variables to stdout
    /*!
     *  The fields are printed out 5 to a line, with each preceded by the radial coordinate value.
     *
     *  @param[in]           title               Title of the printout (can be multiple lines)
     *  @param[in]           indentSpaces        Number of spacies to indent the whole group
     *  @param[in]           radialValues        Value of the radial values of the nodes -> usually starts with 0 and 
     *                                           goes to the edge of the sphere.
     *  @param[in]           numRadialVals       Number of radial nodes
     *  @param[in]           vals                Vector of values to be printed out (cell number is the outer loop)
     *                                           vals[iCell * numFields + iField]
     *  @param[in]           varNames            string names of the fields
     *  @param[in]           numFields           Number of variables to be printed out (this is the inner loop).
     */
    void showOneField(const std::string &title, int indentSpaces, const double* const radialValues, int numRadialVals, 
		      const double* const vals, const std::vector<std::string> &varNames, int numFields);

    //! Print a group of radially varying field variables to stdout. Print the init vs the final for each
    /*!
     *  The fields are printed out 4 to a line, with each preceded by the radial coordinate value. Final and init values are 
     *  printed side-by-side
     *
     *  @param[in]           title               Title of the printout (can be multiple lines)
     *  @param[in]           indentSpaces        Number of spacies to indent the whole group
     *  @param[in]           radialValues        Value of the radial values of the nodes -> usually starts with 0 and 
     *                                           goes to the edge of the sphere.
     *  @param[in]           numRadialVals       Number of radial nodes
     *  @param[in]           vals_init           Vector of init values to be printed out (cell number is the outer loop)
     *                                           vals[iCell * numFields + iField]
     *  @param[in]           vals_final          Vector of final values to be printed out (cell number is the outer loop)
     *                                           vals[iCell * numFields + iField]
     *  @param[in]           varNames            string names of the fields
     *  @param[in]           numFields           Number of variables to be printed out (this is the inner loop).
     */
    void showOneFieldInitFinal(const std::string &title, int indentSpaces, const double * const radialValues, int numRadialVals, 
			const double* const vals_init,  const double * const vals_final,
			const std::vector<std::string> &varNames, int numFields);

    //! Print out the residual and associated tolerances for one field variable that is distributed over the radial coordinate
    /*!
     *  iField and iEqns can be different values
     *
     *  @param[in]           title               Title of the printout (can be multiple lines)
     *  @param[in]           indentSpaces        Number of spacies to indent the whole group
     *  @param[in]           radialValues        Value of the radial values of the nodes -> usually starts with 0 and 
     *                                           goes to the edge of the sphere.
     *  @param[in]           numRadialVals       Number of radial nodes
     *  @param[in]           numFields           Number of fields that are included in the val_init[] and val_final[] array
     *                                           (this is the inner loop).
     *  @param[in]           iField              Actual field to print out the residuals for (0 <= iField < numFields)
     *  @param[in]           val_init            Vector of values to be printed out (cell number is the outer loop)
     *                                           vals[iCell * numFields + iField]
     *  @param[in]           val_final           Vector of values to be printed out (cell number is the outer loop)
     *                                           vals[iCell * numFields + iField]
     *  @param[in]           numEqnsCell         Number of equations per cell included in the residual array
     *                                           (this is the inner loop).
     *  @param[in]           iEqn                Actual Equation index to print out for the residuals for (0 <= iEqn < numEqnsCell)
     *  @param[in]           resid_error         Vector of residual errors (length = numRadialVals * numEqnsCell) 
     *                                           (can be zero, in which case, it's calculated within the routine)
     *  @param[in]           solnError_tol       Solution error tolerance. (length = numRadialVals * numEqnsCell)
     *                                           Can be zero. If existing, its value is printed as AbsResErrorTol 
     *  @param[in]           residual            value of the residual (length = numRadialVals * numEqnsCell)
     */
    void showOneResid(const std::string &title, int indentSpaces, const double * const radialValues, size_t numRadialVals, 
		      size_t numFields, size_t iField, const double * const val_init, const double * const val_final,
		      size_t numEqnsCell, size_t iEqn, const double * const resid_error, const double * const solnError_tol,
		      const double * const residual);

    //! Print the solution for the current step to standard output
    /*!
     *  @param[in]           indentSpaces        Number of spaces to indent each line
     */
    void showSolution(int indentSpaces);

    //! Print the satisfaction of the residual for the current subtimestep integration to the standard output
    /*!
     *  @param[in]           indentSpaces        Number of spaces to indent each line
     *  @param[in]           residual            Residual vector from the residual calculation routine.
     */
    void showResidual(int indentSpaces, const double* const residual);

    //-------------------------------------------------------------------------------------------------------------------
    // -------------------------------- QUERY AND SET MOLE NUMBERS-------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------

    // -----------------------------------------------------------------------------------------------------------------
    // ------------------------------------ CALCULATE INSTANTANEOUS RATES ----------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------
    //     -- These are now called "ProductionRate"  or SourceTerms values

    //! Overpotential term for the heat generation from a single surface for the current global time step (J / s)
    /*!
     *  @param[in]           isk                 Surface index 
     * 
     *  @return                                  Returns the thermal energy overpotential term ( units = J / s) 
     */
    virtual double thermalEnergySourceTerm_overpotential(size_t isk) override;

    //! Reversible Entropy term leading to  heat generation from a single surface for the current global time step (J / s)
    /*!
     *  @param[in]           isk                 Surface index
     * 
     *  @return                                  Returns the thermal energy source term for the global step ( units = J / s) 
     */
    virtual double thermalEnergySourceTerm_reversibleEntropy(size_t isk) override;

    //! Enthalpy Release term leading to heat generation from a single surface for the current global time step (J/s)
    /*!  
     *  (virtual from Electrode.h)
     *  @param[in]           isk                 Surface index
     * 
     *  @return                                  Returns the thermal energy source term for the global step ( units = J / s) 
     */
    virtual double thermalEnergySourceTerm_EnthalpyFormulation(size_t isk) override;

    // -----------------------------------------------------------------------------------------------------------------
    // ------------------------------------ CARRY OUT INTEGRATION OF EQUATIONS -----------------------------------------
    // -----------------------------------------------------------------------------------------------------------------

    // Pack the nonlinear solver proplem
    /*
     *  formulate the nonlinear solver problem to be solved.
     *     Fields to be filled in
     *             yvalNLS_
     *             ylowNLS_
     *             yhighNLS_
     *             deltaBoundsMagnitudesNLS_
     */
    virtual void initialPackSolver_nonlinFunction();

    //! Pack the solution vector
    /*!
     *  @param[in]           y                   solution vector to be filled
     */
    void packNonlinSolnVector(double* const y) const;

    //! Check existing or formulate an initial value for the solution vector
    /*!
     *  @param[in]           doOthers            Fill the other yvalNLS values with the value
     */
    virtual void check_yvalNLS_init(bool doOthers) override;

    //!   Calculate the integrated source terms and do other items now that we have a completed time step
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *  Calculate source terms on completion of a step. At this point we have solved the nonlinear problem
     *  for the current step, and we are calculating post-processed quantities like source terms.
     */
    virtual void calcSrcTermsOnCompletedStep() override;

    //! Determine the species with the largest mole fraction
    /*!
     *    Fill in the array phaseMFBig_[iph] for all phases in the Electrode.
     *    This is currently called once at the start of the problem, from setResidAtolNLS()
     */
    virtual void determineBigMoleFractions() override;

private:    

    //! Carry out work that needs to be done when the solution components are changed
    /*!
     *  Solution components that must be modified when the special species changes:
     *
     *    solnDot_init_
     *
     *  @param[in]           iCell               Value of the cell that is changing
     *  @param[in]           jRPh                Value of the distributed phase that is being changed
     *  @param[in]           oldK                Old value of the special species
     *  @param[in]           bigK                New value of the special species
     */
    void changeBigMF(size_t iCell, size_t jRPh, size_t oldK, size_t bigK);


public:
    //!  Gather the predicted solution values and the predicted integrated source terms
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *  Both the predicted solution values and the predicted integrated source terms are used in the time step control
     */
    virtual void gatherIntegratedSrcPrediction() override;

    //!  Calculate the norm of the difference between the predicted answer and the final converged answer for the current time step
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *   The norm calculated by this routine is used to determine whether the time step is accurate enough. Two norms are taken,
     *   one from the predictedSolution routine and the other from the solnDot prediction. The lesser of the deviation of the
     *   predictions from the final answer is used as the final error predictor 
     *
     *  @param[in]           yvalNLS             Converged final answer for the solution unknowns, yval, from the nonlinear solver
     *                                           for the current time step.
     *
     *  @return                                  Returns the norm of the difference. Normally this is the weighted L0 norm 
     *                                           of the difference between predictor and the corrector.
     *                                           The lesser of the deviation in the two norms is now taken as the answer.
     */
    virtual double predictorCorrectorWeightedSolnNorm(const std::vector<double>& yvalNLS) override;

    //! Calculate the vector of predicted errors in the source terms that this integrator is responsible for
    /*!
     *  (virtual from Electrode_Integrator)
     *   In the base implementation we assume that the there are just one source term, the electron  source term.
     *   However, this will be wrong in almost all cases. The number of source terms is unrelated to the number of 
     *   unknowns in the nonlinear problem. Source terms will have units associated with them.
     *   For example the integrated source term for electrons will have units of kmol.
     */
    virtual void predictorCorrectorGlobalSrcTermErrorVector() override;

    //! Print table representing prediction vs. corrector information
    /*!
     *  @param[in]           yval                Vector of corrector values
     *  @param[in]           pnormSrc            Norm of the predictor-corrector comparison for the source vector.
     *  @param[in]           pnormSoln           Norm of the predictor-corrector comparison for the solution vector.
     */
    virtual void predictorCorrectorPrint(const std::vector<double>& yval, double pnormSrc, double pnormSoln) const override;


    //------------------------------------------------------------------------------------------------------------------
    // ---------------------------- INTEGRATED SOURCE TERM QUERIES -----------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------

    //! Energy released during a single local time step (joules)
    /*!
     *  (virtual from Electrode.h)
     *  Total Energy released within the electrode during a local time step due to all processes
     *
     *  @return                                  Return the energy released (joules)
     */
    virtual double thermalEnergySourceTerm_EnthalpyFormulation_SingleStep() override;

    //! Reversible Entropy release during a single step (joules)
    /*!
     *  (virtual from Electrode.h)
     *  Energy released within the electrode during a local time step due to reversible entropy generation
     *
     *  @return                                  Return the energy released (joules)
     */
    virtual double thermalEnergySourceTerm_ReversibleEntropy_SingleStep() override;

    //! Irreversible thermal energy release during a single step (joules)
    /*!
     *  (virtual from Electrode.h)
     *  Energy released within the electrode during a local time step due to the overpotential 
     *
     *  @return                                  Return the energy released (joules)
     */
    virtual double thermalEnergySourceTerm_Overpotential_SingleStep() override;

    // -----------------------------------------------------------------------------------------------------------------
    // ---------------------------- SOLUTION OF NONLINEAR TIME DEPENDENT SYSTEM  ---------------------------------------
    // -----------------------------------------------------------------------------------------------------------------


    //------------------------------------------------------------------------------------------------------------------
    // -------------------------------  SetState Functions -------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------

    //! Sets the state of the Electrode object given an EState object
    /*!
     *   (virtual function from Electrode)
     *   This sets all of the states within the object to the same state.
     *   It is an error to call this function during a pending step where there can be a difference between t_init and t_final.
     *
     *   @param[in]  es          Const reference to the EState object. Must be an EState_RadialDistrib child
     *                           or else it will throw an error. However, there is an option to 
     *                           read EState objects with less information. 
     */
    virtual void setState_EState(const EState& es) override;


    //! Set the internal initial intermediate and initial global state from the internal final state
    /*!
     *  (virtual function from Electrode.h)
     *
     *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
     *  routine as well.
     *
     *  @param[in]           setInitInit         Boolean indicating whether you should set the init_init state as well
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
     * @param[in]            setFinal            Boolean indicating whether you should set the final as well
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
    virtual void printElectrode(int pSrc = 1, bool subTimeStep = false) override;

    //! Print condition of a phase in the electrode
    /*!
     *  @param[in]           iph                 Print the phase
     *  @param[in]           pSrc                Print Source terms that have occurred during the step from the initial_initial
     *                                           to the final_final time.
     *                                           The default is to print out the source terms
     *  @param[in]           subTimeStep         Print out conditions from the most recent subTimeStep and not the global
     *                                           time step. The default is to print out the global values
     */
    virtual void printElectrodePhase(size_t iph, int pSrc = 1, bool subTimeStep = false) override;

    //! Reset the Internal state of the electrode at the start of a new global integration step or a redo of the current global step
    /*!
     *  (virtual from Electrode)
     *  This function advances the initial state to the final state that was calculated in the last integration step.
     *
     *  @param[in]           Tinitial            This is the New initial time. This time is compared against the "old"
     *                                           final time, to see if there is any problem. They should be the same.
     *                                           If Tinitial == t_init_init, we redo the time step.
     *
     *  @param[in]           doAdvancementAlways Always do the advancement, no matter what. Normally, Tinitial is checked against the 
     *                                           current t_init_init value. If they are the same, then we redo the time step.
     *                                           However, if  doAdvancementAlways is true, we advance the solution unknowns to the 
     *                                           final_final values produced in the last global step no matter what.
     *                                           Defaults to false.
     *
     *  @return                                  Returns true if the time step is reset to t_init_init.
     */
    virtual bool resetStartingCondition(double Tinitial, bool doAdvancementAlways = false) override;

    //! Update the stat of a single phase
    /*!
     *  @param[in]           iph                 Index of the phase within the %PhaseList object
     */    
    void updateState_Phase(size_t iph);

    //! Take the state (i.e., the final state) within the Electrode_SimpleDiff and push it up
    //! to the zero-dimensional parent object
    /*!
     * 
     *  update:
     *             spMoles_final_ [] -> sum soli/d phase species
     *             spMf_final_[]  -> Use exterior cell values
     *
     */
    void updateState_OneToZeroDimensions();

    //!  Set the state of the ThermoPhases to the exterior surface
    /*!
     *   We can set the ThermoPhases using the common temperature and pressure and the mole fraction
     *   to the exterior cell.
     */
    void setState_exteriorSurface();

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

    //! Check total moles of stuff against geometry of particle
    /*!
     *  error exit if inconsistency
     *  works on _final_ state only
     */
    void checkGeometry() const;

    //! Check routine to see if we can account for the mass loss
    /*!
     *   Algorithm is to check for mass loss. If there is some, then add moles back into far last cell.
     *   Note, this works because we have an accounting of all possible sources for
     *
     *  @param[in]           doErr               If true we skip some printing if we find significant discrepancies
     *                                           Defaults to false
     */
    void checkMoles_final_init(bool doErr = false);

    //!  Here we fix the mole balance locally from init state to final state
    /*!
     *  Algorithm is to check for mole loss for each distributed species in the problem.
     *  If there is some, then add moles back into far last cell to ensure complete conservation of moles/mass.
     *  Note, this works because we have an accounting of all possible sources for losses, and because
     *  It has been experimentally verified that the if you crank the nonlinear solver's convergence requirements down far enough
     *  then mass conservation up to round-off error occurs.
     *
     *  HKM experience. This is a local fix for the local subgrid time step.
     *                  Therefore, it must be done at every local step irrespective of a cutoff.
     *                  If this is not done, then the global mole balances may start to drift off of the target.
     *                  When this is done at every step, then this fix is equivalent to fixCapacityBalances_final();
     *
     * @return   Returns the largest mole balance deficit
     */
    double fixMoleBalance_final_init();

    //! Experimental routine to enforce a balance on the electrode capacity given previous balance information
    /*! 
     *  This routine equalizes the capacity 
     */
    virtual void fixCapacityBalances_final();

    //!
    virtual void check_final_state();

    //! Calculate the diffusive flux of all distributed species at the right cell boundary of the cell, iCell
    /*!
     *  Algorithm assumes that species 0 is special. It's usually called the vacancy species. Think of it as the vacency
     *  species. We sum up the diffusive fluxes of all the other species. Then, the diffusive flux of the vacency is calculated
     *  as the negative of that sum. What we are doing is ensuring that the sum of the diffusive flux of all species is equal
     *  to zero. This is necessary if we are to ensure that there are no lattice movements.
     *
     *  The diffusive flux is the based on the gradient of the activity concentration rather than the concentration. 
     *  This difference is significant in many battery systems.
     *
     *  @param[out]          fluxRCB             diffusive flux for all distributed species at the right cell boundary
     *  @param[in]           iCell               Cell number
     *  @param[in]           finalState          If true we use the _final_ values of the concentrations.
     *                                           if false, we use the _init_ values of the concentrations; this is used
     *                                           to make predictions of the flux. 
     */
    void diffusiveFluxRCB(double* const fluxRCB, size_t iCell, bool finalState) const;  

    //! Predict the solution
    /*!
     *  Ok at this point we have a time step deltalimiTsubcycle_
     *  and initial conditions consisting of phaseMoles_init_ and spMF_init_.
     *  We now calculate predicted solution components from these conditions.
     *
     * @return                                   Returns the success of the operation
     *                                            1  A predicted solution is achieved
     *                                            2  A predicted solution with a multispecies phase pop is acheived
     *                                            0  A predicted solution is not achieved, but go ahead anyway
     *                                           -1  The predictor suggests that the time step be reduced and a retry occur.
     */
    int predictSolnResid();

    //! Predict the solution at the next time step
    /*!
     *  Ok at this point we have a time step deltalimiTsubcycle_
     *  and initial conditions consisting of phaseMoles_init_ and spMF_init_.
     *  We now calculate predicted solution components from these conditions.
     *
     *  Calls predictSolnResid(). If it works that good, if not reduce the time step and call the procedure again.
     *
     *  @return                                  Returns the success of the operation
     *                                            1  A predicted solution is achieved
     *                                            2  A predicted solution with a multispecies phase pop is acheived
     *                                            0  A predicted solution is not achieved, but go ahead anyway
     *                                           -1  The predictor suggests that the time step be reduced and a retry occur.
     */
    virtual int predictSoln() override;

    //! Predict the solution using the time derivative of the solution from the last step
    /*!
     *  (virtual from Electrode_Integrator)
     * 
     *  Predicts the solution at the final time from the current derivative of the solution at the initial time.
     *
     *   @return                                 Returns the success of the operation
     *                                           -  1  A predicted solution is achieved
     *                                           -  2  A predicted solution with a multispecies phase pop is acheived
     *                                           -  0  A predicted solution is not achieved, but go ahead anyway
     *                                           - -1  The predictor suggests that the time step be reduced and a retry occur.
     */
    virtual int predictSolnDot() override;
  
    //! Unpack the solution vector on return from the time stepper
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *  This function unpacks the solution vector into  phaseMoles_final_,  spMoles_final_, and spMf_final_[]
     *
     *  @param[in]           ySoln               Solution vector as returned from the time stepper.
     *                                           What the solution actually refers to depends on the individual Electrode objects.
     *                                           y[0] is always the current  deltaTsubcycleCalc_ value
     *
     *  @return                                  Returns the success of the operation
     *                                             1 success
     *                                             0 failure - bounds problems
     */
    virtual int unpackNonlinSolnVector(const double* const ySoln) override;

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
     *   @param[in]          rtolResid           Relative tolerance allowed for the electron source term over the interval.
     *                                           This is a unitless quantity.
     */
    virtual void setNLSGlobalSrcTermTolerances(double rtolResid) override;

    //!   Set the Residual absolute error tolerances
    /*!
     *  (virtual from Electrode_Integrator)
     *   Set the absolute error tolerances for the residuals for the nonlinear solvers. This is called at the top
     *   of the integrator() routine.
     *
     *   Calculates residAtolNLS_[]
     *   Calculates atolNLS_[]
     */
    virtual void setResidAtolNLS() override;

    //! Evaluate the residual function
    /*!
     * (virtual from NonlinearSolver)
     *
     * @param[in]            t                   Time                    (input)
     * @param[in]            delta_t             The current value of the time step (input)
     * @param[in]            y                   Solution vector (input, do not modify)
     * @param[in]            ydot                Rate of change of solution vector. (input, do not modify)
     * @param[out]           resid               Value of the residual that is computed (output)
     * @param[in]            evalType            Type of the residual being computed (defaults to Base_ResidEval)
     * @param[in]            id_x                Index of the variable that is being numerically differenced to find
     *                                           the jacobian (defaults to -1, which indicates that no variable is being
     *                                           differenced or that the residual doesn't take this issue into account)
     * @param[in]            delta_x             Value of the delta used in the numerical differencing
     *
     * @return                                   Returns an integer that gets fed back through evalResidNJ() to the
     *                                           nonlinear solver. Anything other than a 1 causes an immediate failure
     *                                           of the nonlinear solver to occur.
     */
    virtual int evalResidNJ(const double t, const double delta_t, const double* const y, const double* const ydot,
                            double* const resid, const ResidEval_Type evalType = ResidEval_Type::Base_ResidEval,
                            const int id_x = -1, const double delta_x = 0.0) override;

    //!  Residual calculation for the solution of the Nonlinear integration problem
    /*!
     *    Given tfinal and delta_t, and given y and ydot which are estimates of the solution
     *    and solution derivative at tfinal, this function calculates the residual equations.
     *    It is the residual function used in the nonlinear solver that relaxes the equations at each time step.
     *
     *    This is typically called from evalResidNJ(), which is called directly from the
     *    nonlinear solver. However, we expose this routine so that the residual can be queried given all of the inputs.
     *
     *  @param[in]           tfinal              Time                    (input)
     *  @param[in]           delta_t             The current value of the time step (input)
     *  @param[in]           y                   Solution vector (input, do not modify)
     *  @param[in]           ySolnDot            Rate of change of solution vector. (input, do not modify)
     *  @param[out]          resid               Value of the residual that is computed (output)
     *  @param[in]           evalType            Type of the residual being computed (defaults to Base_ResidEval)
     *  @param[in]           id_x                Index of the variable that is being numerically differenced to find
     *                                           the jacobian (defaults to -1, which indicates that no variable is being
     *                                           differenced or that the residual doesn't take this issue into account)
     *  @param[in]           delta_x             Value of the delta used in the numerical differencing
     *
     *  @return                                  Returns an integer that gets fed back through evalResidNJ() to the
     *                                           nonlinear solver. Anything other than a 1 causes an immediate failure
     *                                           of the nonlinear solver to occur.
     */
    virtual int integrateResid(const double tfinal, const double delta_t, const double* const y, const double* const ySolnDot,
                               double* const resid, const ResidEval_Type evalType, const int id_x,
                               const double delta_x) override;

    //! Additional bounds checking
    /*!
     *  Do any problem specific bounds checking.
     *
     *  @param[in]           t                   Time                    (input)
     *  @param[in]           ybase               Solution vector (input, output)
     *  @param[in]           step                Proposed step in the solution that will be cropped
     *  @param[in]           dampIn              Existing damping coefficient already applied
     *
     *  @return                                  Returns the damping coefficient needed for this problem. Will be less
     *                                             than or equal to the input, dampIn
     */
    virtual doublevalue boundsCheckAddn(const doublevalue t, const doublevalue* const ybase, const doublevalue* const step, 
                                        const doublevalue dampIn) override;

    //! Main internal routine to calculate the residual
    /*!
     *  This routine calculates the functional at the current stepsize, deltaTsubcycle_.
     *  A new stepsize, deltaTsubcycleCalc_, is calculated within this routine for changes
     *  in topology of type LimitingEventDeltaTsubcycle_
     *
     *  This routine calculates resid[i], which is the calculated value of the residual for the nonlinear function:
     *
     *   resid[i] = y[i] - yval_retn[i]
     *
     *  The formulation of the residual and the corresponding independent unknown is described below:
     *
     *   Residual entry:                                            Unknown                           Index
     * --------------------------------------------------------------------------------------------------------------
     *         Residual (Time)                                     deltaSubcycleCalc_                   0
     *                                                                                            1
     *         Loop over cells                                                            0 <=  iCell < numRCells_
     *                                                                                     j = numEqnsCell_ * iCell
     *                    
     *            Loop over distributed Phases
     *            Residual (Concentration _ k=0)                  totMoles_SPhase__Cell_final_[iCell * numSPhases_ + jPh]
     *              . . .
     *            Residual (Concentration _ k=Ns-1)               spMoles_KRSolid_final_[iCell * numKRSpecies_ + iKRSpecies]
     *  --------------------------------------------------------------------------------------------------------------
     *
     *   Previously, we should have called:
     *        unpackNonlinSolnVector()
     *        updateState()
     *        extractInfo(); 
     *        updateSpeciesMoleChangeFinal();
     *
     *  @param[in]           resid               Residual vector. Length = neq_
     *
     *  @param[in]           evalType            Type of the residual being computed (defaults to Base_ResidEval)
     *                                           Other types are JacBase_ResidEval, JacDelta_ResidEval, and Base_ShowSolution.
     *
     *  @return                                  Returns a value of 1 if everything went well.
     *                                           Returns negative numbers to indicate types of failures.
     */
    virtual int calcResid(double* const resid, const ResidEval_Type evalType) override;

    //! Returns the equilibrium OCV for the selected ReactingSurfaceDomain and current conditions
    /*!
     *  (virtual from Electrode)
     *  This routine uses a root finder to find the voltage at which there
     *  is zero net electron production.  It leaves the object unchanged. However, it
     *  does change the voltage of the phases during the calculation, so this is a non const function.
     *
     * @param[in]       isk                            Reacting surface domain id
     * @param[in]       comparedToReferenceElectrode   Boolean indicating whether voltage is referenced to the solution at
     *                                                 the current conditions (false) or compared to the voltage wrt the 
     *                                                 reference electrode (true). The later is akin to using the standard 
     *                                                 state thermo functions for the electrolyte species.
     *
     *   @return                                       Returns the OCV (volts)
     */
    virtual double openCircuitVoltage(size_t isk, bool comparedToReferenceElectrode = false) override;

    //! Get the open circuit potential at the mixture averaged conditions of the electrode
    /*!
     *  (virtual from Electrode)
     *  This routine creates a mixture averaged condition of the electrode (eliminating any diffusion process)
     *  before calculating the OCV. The result is the same as openCircuitVoltage() in the base class.
     *
     *  @param[in]     isk                           Reacting surface domain id
     *  @param[in]     comparedToReferenceElectrode  Boolean, if true compare to the reference electrode. Defaults to false.  
     *
     *  @return                                      Returns the OCV (volts)
     */
    virtual double openCircuitVoltage_MixtureAveraged(size_t isk, bool comparedToReferenceElectrode = false) override;

#ifdef DEBUG_THERMAL
 
    double netElectrons();
#endif

protected:

    // ------------------------------------------------------- D A T A -------------------------------------------------------------

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

    //! Define the number of species that are defined to have radially distributed distributions within the solid
    /*!
     *   Note, for this object there is only one radial distribution. This concept will be enhanced
     *   in later formulations.
     *
     *   There are a few arrays which have numKRSpecies_ as their inner loop dimension.
     */
    size_t numKRSpecies_;

    //!  Number of cells involved with the radial distribution, including the 1/2 end cells
    size_t numRCells_;

    //! Number of phases which have radial distributions of their species
    size_t numSPhases_;

    //! Total number of equations defined at each node of the radial mesh
    size_t numEqnsCell_;


    //! Pointers to the ThermoPhase objects that make up the radially distributed portion of this Electrode object
    /*!
     *  No attempt is made to malloc a ThermoPhase object for each cell within the radial distribution
     *  Length: numSPhases_
     */
    std::vector<ThermoPhase *> thermoSPhase_List_;

    //! Vector of solid species moles defined on the grid of the spherical particle
    /*!
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies;
     *   final state value
     *
     *   Length: numRCells_ * numKRSpecies_
     *   units: kmol 
     */
    std::vector<double> spMoles_KRsolid_Cell_final_;

    //! Vector of solid species defined on the grid of the spherical particle
    /*!
     *
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies;
     *   init state value
     *
     *   Length: numRCells_ * numKRSpecies_
     *   units: kmol 
     */
    std::vector<double> spMoles_KRsolid_Cell_init_;

    //! Vector of solid species defined on the grid of the spherical particle
    /*!
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies;
     *   final_final state value
     *
     *   Length: numRCells_ * numKRSpecies_
     *   units: kmol 
     */
    std::vector<double> spMoles_KRsolid_Cell_final_final_;

    //! Vector of solid species defined on the grid of the spherical particle
    /*!
     *
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies;
     *   init_init state value
     *
     *   Length: numRCells_ * numKRSpecies_
     *   units: kmol 
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
     *
     *  Length: numKRSpecies_
     */
    std::vector<size_t> KRsolid_speciesList_;

    //! Species Names for the solid species that are distributed radially through the particle
    std::vector<std::string> KRsolid_speciesNames_;

    //! Phase indeces of the solid phases comprising the species that are radially distributed
    /*!
     *  The phase index value is the index within the Electrode object's PhaseList
     *  
     *  iPh = phaseIndeciseKRsolidPhases_[distribPhIndex];
     *
     *  Length: numSPhases_
     */
    std::vector<size_t> phaseIndeciseKRsolidPhases_;

    //! Returns the distributed phase index given the regular phase index.
    /*!
     *  distribPhIndex = phaseIndeciseKRsolidPhases_[iPh];
     *
     *  If distribPhIndex is -1, the phase isn't distributed.
     *
     *  Length: m_NumTotPhases
     */
    std::vector<size_t> distribPhIndexKRsolidPhases_;

    //! Number of species in each of the radially distributed phases
    /*!
     *  The index is the distributed phase index number
     *
     *  Length: numSPhases_
     */
    std::vector<size_t> numSpeciesInKRSolidPhases_;

    //!  Index of species starting position for each distributed phase within the distributed species vectors
    /*!
     *  Length: numSPhases_
     */
    std::vector<size_t> kstartKRSolidPhases_;

    //! Phase Names of the solid phases comprising the species that are radially distributed
    /*!
     *  Length: numSPhases_
     */
    std::vector<std::string> KRsolid_phaseNames_;

    //! Big Mole Fractions within each cell within each distributed phaes
    /*!
     *
     *   Usage: phaseMFBig_SPhase_Cell_[numSPhases_ * iCell + iJRPh]
     *   Length: numRCells_ * numSPhases_
     *   units: unitless
     */
    std::vector<size_t> phaseMFBig_SPhase_Cell_;

    //! Phase indeces of the solid phases comprising the species that are not radially distributed
    /*!
     *  The phaseIndex is the index within the Electrode object
     *
     *  Length: numNonSPhases_
     */
    std::vector<size_t> phaseIndeciseNonKRsolidPhases_;

    //! Number of phases that are not distributed
    size_t numNonSPhases_;

    //! Total concentration of each of the solid phases that are distributed - local final state
    /*!
     *   concTot_SPhase_Cell_final_[numSPhases_ * iCell + iJRPh]
     *   Length: numRCells_ * numSPhases_
     *   units: kmol / m3
     */
    std::vector<double> concTot_SPhase_Cell_final_;

    //! Total concentration of each of the solid phases that are distributed - local init state
    /*!
     *   concTot_SPhase_Cell_init_[numSPhases_ * iCell + iJRPh]
     *   Length: numRCells_ * numSPhases_
     *   units: kmol / m3
     */
    std::vector<double> concTot_SPhase_Cell_init_;

    //! Total concentration of each of the solid phases that are distributed - global final state
    /*!
     *   concTot_SPhase_Cell_final_final_[numSPhases_ * iCell + iJRPh]
     *   Length: numRCells_ * numSPhases_
     *   units: kmol / m3
     */
    std::vector<double> concTot_SPhase_Cell_final_final_;

    //! Total concentration of each of the solid phases that are distributed - global init state
    /*!
     *   concTot_SPhase_Cell_init_init_[numSPhases_ * iCell + iJRPh]
     *   Length: numRCells_ * numSPhases_
     *   units: kmol / m3
     */
    std::vector<double> concTot_SPhase_Cell_init_init_;

    //! Total moles of each of the solid phases that are distributed - local final state
    /*!
     *   phaseMoles_KRsolid_Cell_final_[numSPhases_ * iCell + iJRPh]
     *   Length: numRCells_ * numSPhases_
     *   units: kmol in all of the particles
     */
    std::vector<double> phaseMoles_KRsolid_Cell_final_;

    //! Total moles of each of the solid phases that are distributed - local init state
    /*!
     *   phaseMoles_KRsolid_Cell_init_[numSPhases_ * iCell + iJRPh]
     *   Length: numRCells_ * numSPhases_
     *   units: kmol in all of the particles
     */
    std::vector<double> phaseMoles_KRsolid_Cell_init_;

    //! Total moles of each of the solid phases that are distributed - global final state
    /*!
     *   phaseMoles_KRsolid_Cell_final_final_[numSPhases_ * iCell + iJRPh]
     *   Length: numRCells_ * numSPhases_
     *   units: kmol in all of the particles
     */
    std::vector<double> phaseMoles_KRsolid_Cell_final_final_;

    //! Total moles of each of the solid phases that are distributed per cell - global init state
    /*!
     *   phaseMoles_KRsolid_init_init_[numSPhases_ * iCell + iJRPh]
     *   Length: numRCells_ * numSPhases_
     *   units: kmol in all of the particles
     */
    std::vector<double> phaseMoles_KRsolid_Cell_init_init_;

    //! total concentration of the solid phases that are distributed - init state
    /*!
     *   concKRSpecies_Cell_init_[numKRSpecies_ * iCell + iKRSpecies]
     *   Length: numRCells_ * numKRSpecies_
     *   units: kmol/m3 
     */
    std::vector<double> concKRSpecies_Cell_init_;

    //! total concentration of the solid phases that are distributed - init state
    /*!
     *   concKRSpecies_Cell_final_[numKRSpecies_ * iCell + iKRSpecies]
     *   Length: numRCells_ * numKRSpecies_
     *   units: kmol/m3 
     */
    std::vector<double> concKRSpecies_Cell_final_;

    //! total concentration of the solid phases that are distributed - init state
    /*!
     *   concKRSpecies_Cell_init_init_[numKRSpecies_ * iCell + iKRSpecies]
     *   Length: numRCells_ * numKRSpecies_
     *   units: kmol/m3 
     */
    std::vector<double> concKRSpecies_Cell_init_init_;

    //! total concentration of the solid phases that are distributed - init state
    /*!
     *   concKRSpecies_Cell_final_final_[numKRSpecies_ * iCell + iKRSpecies]
     *   Length: numRCells_ * numKRSpecies_
     *   units: kmol/m3 
     */
    std::vector<double> concKRSpecies_Cell_final_final_;

    //! Mole fraction of the solid phase species that are distributed = final state
    /*!
     *   spMf_KRSpecies_Cell_final_[numKRSpecies_ * iCell + iKRSpecies]
     *   Length: numRCells_ * numKRSpecies_
     *   units:  mole fractions
     */
    std::vector<double> spMf_KRSpecies_Cell_final_;

    //! Mole fraction of the solid phase species that are distributed - final state of subtimestep integration
    /*!
     *   spMf_KRSpecies_Cell_init_[numKRSpecies_ * iCell + iKRSpecies]
     *   Length: numRCells_ * numKRSpecies_
     *   units:  mole fractions
     */
    std::vector<double> spMf_KRSpecies_Cell_init_;

    //! Mole fraction of the solid phase species that are distributed = final state
    /*!
     *   spMf_KRSpecies_Cell_final_final_[numKRSpecies_ * iCell + iKRSpecies]
     *   Length: numRCells_ * numKRSpecies_
     *   units:  mole fractions
     */
    std::vector<double> spMf_KRSpecies_Cell_final_final_;

    //! Mole fraction of the solid phase species that are distributed = final state
    /*!
     *   spMf_KRSpecies_Cell_init_init_[numKRSpecies_ * iCell + iKRSpecies]
     *   Length: numRCells_ * numKRSpecies_
     *   units:  mole fractions
     */
    std::vector<double> spMf_KRSpecies_Cell_init_init_;

    //! Node position of the mesh - final_final
    /*!
     *  Length: numRCells_
     *  units:  m  (particle radius coordinate)
     */
    std::vector<double> rnodePos_final_final_;

    //! Node position of the mesh - final
    /*!
     *  Length: numRCells_
     *  units:  m  (particle radius coordinate)
     */
    std::vector<double> rnodePos_final_;

    //! Node position of the mesh - init
    /*!
     *  Length: numRCells_
     *  units:  m  (particle radius coordinate)
     */
    std::vector<double> rnodePos_init_;

    //! Node position of the mesh - init_init
    /*!
     *  Length: numRCells_
     *  units:  m  (particle radius coordinate)
     */
    std::vector<double> rnodePos_init_init_;

    //! Location of the right boundary of each cell at the end time of subtimestep integration
    /*!
     *  Length: numRCells_
     *  units:  m  (particle radius coordinate)
     */
    std::vector<double> cellBoundR_final_;
  
    //! Location of the right boundary of each cell at the init time of subtimestep integration
    /*!
     *  Length: numRCells_
     *  units:  m  (particle radius coordinate)
     */
    std::vector<double> cellBoundR_init_;

    //! Location of the right boundary of each cell at the final time of global timestep integration
    /*!
     *  Length: numRCells_
     *  units:  m  (particle radius coordinate)
     */
    std::vector<double> cellBoundR_final_final_;

    //! Location of the right boundary of each cell at the init time of global timestep integration
    /*!
     *  Length: numRCells_
     *  units:  m  (particle radius coordinate)
     */
    std::vector<double> cellBoundR_init_init_;

    //! Location of the left boundary of each cell at the final time of subtimestep integration
    /*!
     *  Length: numRCells_
     *  units:  m  (particle radius coordinate)
     */
    std::vector<double> cellBoundL_final_;

    //! Location of the left boundary of each cell at the init time of subtimestep integration
    /*!
     *  Length: numRCells_
     *  units:  m  (particle radius coordinate)
     */
    std::vector<double> cellBoundL_init_;

    //! Location of the left boundary of each cell at the init time of global timestep integration
    /*!
     *  Length: numRCells_
     *  units:  m  (particle radius coordinate)
     */
    std::vector<double> cellBoundL_init_init_;

    //! Location of the left boundary of each cell at the final time of global timestep integration
    /*!
     *  Length: numRCells_
     *  units:  m  (particle radius coordinate)
     */
    std::vector<double> cellBoundL_final_final_;

    //! Volume of each cell on a per particle basis
    /*!
     *  Length: numRCells_
     *  units:  m3 (particle volume basis)
     */
    std::vector<double> volPP_Cell_final_;

    //!  Spline System for the nodal points
    /*!
     *  Fraction of the domain's volume which is inside the current node position
     *  Length: numRCells_
     *  units:  unitless
     */
    std::vector<double> fracVolNodePos_;

    //!  Partial molar volume of all of the solid species located in all of the cells
    /*!
     *   Vector of molar volumes for all solid species in all cells (KRSpecies, iCell)
     *   Length: numKRSpecies_ *  numRCells_
     *   Units:  m3 / kmol
     */
    std::vector<double> partialMolarVolKRSpecies_Cell_final_;

    //!  Partial molar Heat Capacity  of all of the solid species located in all of the cells
    /*!
     *   Vector of partial molar heat capacity const press (KRSpecies, iCell)
     *   Length: numKRSpecies_ *  numRCells_
     *   Units:  Joules/(kmol K)
     */
    mutable std::vector<double> partialMolarCpKRSpecies_Cell_final_;

    //!  Partial molar Enthalpy of all of the solid species located in all of the cells
    /*!
     *   Vector of partial molar Enthalpy  (KRSpecies, iCell)
     *   Length: numKRSpecies_ *  numRCells_
     *   Units:  Joules/(kmol)
     */
    mutable std::vector<double> partialMolarEnthKRSpecies_Cell_final_;

    //!  Partial molar Enthalpy of all of the solid species located in all of the cells
    /*!
     *   Vector of partial molar Enthalpy  (KRSpecies, iCell)
     *   Length: numKRSpecies_ *  numRCells_
     *   Units:  Joules/(kmol)
     */
    mutable std::vector<double> partialMolarEnthKRSpecies_Cell_init_;

    //!  Partial molar chemical potential of all of the solid species located in all of the cells
    /*!
     *   Vector of partial molar Enthalpy  (KRSpecies, iCell)
     *   Length: numKRSpecies_ *  numRCells_
     *   Units:  Joules/(kmol)
     */
    mutable std::vector<double> chemPotKRSpecies_Cell_final_;

    //!  Partial molar chemical potential of all of the solid species located in all of the cells
    /*!
     *   Vector of partial molar Enthalpy  (KRSpecies, iCell)
     *   Length: numKRSpecies_ *  numRCells_
     *   Units:  Joules/(kmol)
     */
    mutable std::vector<double> chemPotKRSpecies_Cell_init_;

    //!  Partial molar Entropy of all of the solid species located in all of the cells
    /*!
     *   Vector of partial molar Entropy  (KRSpecies, iCell)
     *   Length: numKRSpecies_ *  numRCells_
     *   Units:  Joules/(kmol K)
     */
    mutable std::vector<double> partialMolarEntropyKRSpecies_Cell_final_;

    //!  Partial molar Entropy of all of the solid species located in all of the cells
    /*!
     *   Vector of partial molar Entropy  (KRSpecies, iCell)
     *   Length: numKRSpecies_ *  numRCells_
     *   Units: Joules/(kmol K)
     */
    mutable  std::vector<double> partialMolarEntropyKRSpecies_Cell_init_;

    //! Rate of progress of the surface reactions
    /*!
     *   Length: maxNumRxns;
     *   Units:  kmol/m2/s
     */
    std::vector<double> ROP_;

    //!  Molar creation rate of species in the electrode object due to the Exterior surface
    /*!
     *   This is a quantity over all species in the PhaseList object
     *
     *   Length: m_NumTotSpecies in PhaseList
     *   Indexing:  kGlob (indexing of species in the PhaseList)
     *   Units:  kmol sec-1;
     */
    std::vector<double> DspMoles_final_;

    //! Domain boundary at the inner radius.
    /*!
     *   Frequently this will be zero. default is zero.
     */
    double m_rbot0_;

    //! Solid phase diffusion coefficients of the species which are radially diffusing
    /*!
     *  Length: numKRSpecies_
     *  Units:  m2/s
     */
    std::vector<double> Diff_Coeff_KRSolid_;

    //! Molar creation rate of phases in the electrode object.
    /*!
     *  This is calculated as a sum of species creation rates summed up over all cells in the domain
     *
     *  Length: m_NumtotPhases
     *  Units:  kmol / s
     */
    std::vector<double> DphMolesSrc_final_;

    //! we identify the phases here as being the exterior surface
    /*!
     *  The phase will be in direct contact with the electrolyte.
     */
    int surfIndexExteriorSurface_;

    //!  Value of the total flux at the outer edge - kmol m-2 s-1
    double NTflux_final_;

    //!  Model for the formulation of the diffusive flux
    /*!
     *   0 = diffusive flux include the activity coefficient
     *   1 = diffusive flux doesn't include the activity coefficient
     */
    int diffusiveFluxModel_;

    //! Local value of the diffusion coefficient
    double DiffCoeff_;

    //! Local value of the diffusion coefficient
    double DiffCoeff_default_;

    //! Vector of activity coefficients for all KR species at all nodes
    /*!
     *    This is calculated at the final state
     *
     *   Length: numKRSpecies_ *  numRCells_
     *   Units:  unitless
     */
    std::vector<double> actCoeff_Cell_final_;

    //! Vector of activity coefficients for all KR species at all nodes
    /*!
     *    This is calculated at the init state
     *
     *   Length: numKRSpecies_ *  numRCells_
     *   Units:  unitless
     */
    std::vector<double> actCoeff_Cell_init_;

    //! Temporary storage for solution vector
    /*!
     *   Length: ResidEvalJac::neq_
     */
    std::vector<double> m_yvalNLS_Alt;

    //! phase id of the phase which will die at the shortest time
    int phaseID_TimeDeathMin_;

    //! Cell id of the shortest time to a phase death
    int cellID_TimeDeathMin_;

    //! This integer describes if the system is currently on a Region boundary at the start of a subgrid integration step
    /*!  
     *  We define a region boundary here iff all cells are on that boundary
     *  If it isn't, we set it to -1. If it is, we set it to the region boundary index
     */
    int onRegionBoundary_init_;

    //! This integer describes if the system is currently on a Region boundary at the end of a subgrid
    //! integration step
    /*!
     *   We define a region boundary here iff all cells are on that boundary
     *  If it isn't, we set it to -1. If it is, we set it to the region boundary index
     */
    int onRegionBoundary_final_;

    //! This integer describes if the system is currently on a Region boundary at the start of a global
    //! integration step
    /*!
     *  If it isn't, we set it to -1. If it is, we set it to the region boundary index
     */
    int onRegionBoundary_init_init_;

    //! This integer describes if the system is currently on a Region boundary at the end of a global
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

public:

    //! Total number of lattice sites within the electrode at the end of the current subgrid time step (kmol)
    /*!
     *  For many problems used by this class, this is a constant
     */
    double numLattices_final_;

    //! Total number of lattice sites within the electrode at the start of the current time step (kmol)
    /*!
     *  For many problems used by this class, this is a constant
     *  Units: kmol
     */
    double numLattices_init_;

    //! Total number of lattice sites within the electrode predicted for the end of the current subgrid time step
    /*!
     *  For many problems used by this class, this is a constant
     */
    double numLattices_pred_;

    //! Total number of lattice sites within each Cell at the start of the current electrode subcycle time step
    /*!
     *  Length: numRCells_
     *  Units:  kmol
     */
    std::vector<double> numLatticeCBR_init_;

    //! Total number of lattice sites within each Cell at the end of the current electrode subcycle time step
    /*!
     *  Length: numRCells_
     *  Units:  kmol
     */
    std::vector<double> numLatticeCBR_final_;

    friend class ZZCantera::EState;
    friend class ZZCantera::EState_RadialDistrib;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
