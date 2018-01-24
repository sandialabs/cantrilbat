/*
 *  @file Electrode_CSTR.h
 *     Headers for the declarations of the Electrode_CSTR class, used to model 
 *     Electrode processes in particles with no transport limits
 *     (see \ref electrode_mgr and class \link Zuzax::Electrode_CSTR Electrode_CSTR\endlink).
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_CSTR_H
#define _ELECTRODE_CSTR_H


#include "Electrode_input.h"
#include "Electrode_Integrator.h"

#include <string>
#include <vector>

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//! Definition for an anode type of electrode
#define ELECTRODETYPE_ANODE   0

//! Definition for a cathode type of electrode
#define ELECTRODETYPE_CATHODE 1

#define DEBUG_ELECTRODE_DEATH 1

//==================================================================================================================================
//! Extra input for the CSTR Model
/*!
 *
 */
class ELECTRODE_CSTR_KEY_INPUT : public ELECTRODE_KEY_INPUT
{
public:

    //! Constructor
    /*!
     *  @param[in]           printLvl            Print level. Defaults to zero.
     */
    ELECTRODE_CSTR_KEY_INPUT(int printLvl = 0);

    //! Destructor
    virtual ~ELECTRODE_CSTR_KEY_INPUT();

    //! First pass through the child setup system
    /*!
     *  Typically we will fill in all vectors that depend on the value of numRegions_ in this  pass.
     *
     *  @param[in]           cf                  Pointer to the BlockEntry record
     */
    void setup_input_child1(BEInput::BlockEntry* cf);

    //! Second pass through the child setup system
    /*!
     *  Typically we will fill in all vectors that depend on the value of numRegions_ in this  pass.
     *
     *  @param[in]           cf                  Pointer to the BlockEntry record
     */
    void setup_input_child2(BEInput::BlockEntry* cf);

    //! Resistance of the boundary region
    /*!
     *  Units:    ohms
     *  default: 0.0  ohms
     */
    double boundaryResistance_;

    //! Minimum particle radius for purposes of calculating the reaction surface area
    /*!
     *  Units:    meters
     *  default:  1.0E-7 m
     */
    double Radius_exterior_min_;
};
//==================================================================================================================================

//! Electrode_CSTR class is an electrode that models a solid electrode particle as a CSTR
/*!
 *  The base class is close to being a CSTR. This class will have to ensure
 *  that there are no plateaus that create complicated morphologies.
 *
 *  There will just be an interfacial reaction on the exterior of the particle.
 */
class Electrode_CSTR : public Electrode_Integrator
{
public:
    //! Constructor
    Electrode_CSTR();

    //! Destructor
    virtual ~Electrode_CSTR();

    //! Copy Constructor
    /*!
     * @param[in]            right               Object to be copied
     */
    Electrode_CSTR(const Electrode_CSTR& right);

    //! Assignment operator
    /*!
     *  @param[in]           right               object to be copied
     *
     *  @return                                  Returns a reference to the current object
     */
    Electrode_CSTR& operator=(const Electrode_CSTR& right);

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
     *  @return                                  Returns an enum type, called   Electrode_Types_Enum
     */
    virtual Electrode_Types_Enum electrodeType() const override;

    //! This function will extend the input options for the electrode to included child input
    //! and then return a child object of ELECTRODE_KEY_INPUT
    /*!
     *  @param[in]           ei                  Base input for the electrode
     *                                           On output it will be an ELECTRODE_CSTR_KEY_INPUT object
     *
     *  @return                                  Returns zero on success, -1 on failure
     */
    virtual int electrode_input_child(ELECTRODE_KEY_INPUT** ei) override;

    //!  Setup the electrode using the ELECTRODE_KEY_INPUT object that is read from an input file
    /*!
     *   (virtual from Electrode)
     *
     *   @param[in]          ei                  ELECTRODE_KEY_INPUT pointer object
     *
     *   @return                                 Returns 0 if successful, -1 if not.
     */
    virtual int electrode_model_create(ELECTRODE_KEY_INPUT* ei) override;

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
     *  @param[in]           eibase              ELECTRODE_KEY_INPUT pointer object
     *
     *  @return                                  Returns zero if successful, and -1 if not successful.
     */
    virtual int setInitialConditions(ELECTRODE_KEY_INPUT* eibase) override;

    //! local routine to resize arrays that this object is responsible for
    void init_sizes();

    //! Set the sizes of the electrode from the input parameters
    /*!
     *
     * @param electrodeArea   Area of the electrode
     * @param electrodeThickness  Width of the electrode
     * @param porosity        Volume of the electrolyte phase
     */
    virtual void setElectrodeSizeParams(double electrodeArea, double electrodeThickness, double porosity) override;

protected:
    //! Resize the solid phase and electrolyte mole numbers within the object
    /*!
     *  This routine uses particleDiameter_ , particleNumberToFollow_, and porosity_ to recalculate
     *  all the mole numbers in the electrode. This is done by rescaling all of the numbers.
     *  At the end of the process, the total volume of the electrode object is
     *
     *    grossVol = SolidVol() / ( 1.0 - porosity_)
     *
     *  where the SolidVol() is equal to
     *
     *   SolidVol() =  particleNumberToFollow_  Pi *  particleDiameter_**3 / 6.0;
     *
     */
    virtual void resizeMoleNumbersToGeometry() override;

public:

    //! Create and malloc the solvers
    /*!
     * (virtual from Electrode_Integrator)
     *
     * @return                                   returns 1 if ok
     */
    virtual int create_solvers() override;

    //!  Value of the standard state open circuit voltage for the standard state conditions for region xRegion at the _final_
    //!  conditions.
    /*!
     *  This is the standard state open circuit potential at the current _final_ conditions.
     *  The value is dependent on the region input. The region input, currently is identified with the surface reacting phase object.
     *  Therefore, this call is basically a wrapper around openCircuitVoltageSS(isk) which
     *  calculates the open circuit standard state voltage for the isk reacting surface.
     *
     *  Additions include extening the region values to include false values for DoD = 0 and 1 conditions.
     *
     *  @param[in]           xRegion            Value of the region. If -1, this is at the DoD = 0. If nR+1,
     *                                          this is at the DoD = 1.0 condition
     *
     *  @return                                 Returns the voltage of the standard state for that region.
     */
    double openCircuitVoltageSS_Region(int xRegion) const;

    //! Value of the open circuit voltage for region xRegion at the final_ conditions.
    /*!
     *  This is the open circuit potential at the current _final_ conditions.
     *  The value is dependent on the region input. The region input, currently is identified with the surface reacting phase object.
     *  Therefore, this call is basically a wrapper around openCircuitVoltage(isk) which
     *  calculates the open circuit voltage for the isk reacting surface.
     *
     *  Additions include extending the region values to include false values for DoD = 0 and 1 conditions.
     *
     *  @param[in]           xRegion             Value of the region. If -1, this is at the DoD = 0. If nR+1,
     *                                           this is at the DoD = 1.0 condition
     *
     *  @param[in]           comparetoReferenceElectrode  compare against the reference electrode voltage and don't 
     *                                                    just return the volage. Defaults to false
     *
     *  @return                                  Returns the voltage across the interface 
     *                                           (or the voltage compare to the reference electrode voltage) 
     */
    double openCircuitVoltage_Region(int xRegion, bool comparetoReferenceElectrode = false) const;

    //! Returns the total capacity of the electrode in Amp seconds = coulombs
    /*!
     *  Returns the capacity of the electrode in Amps seconds.
     *  This is the same as the number of coulombs that can be delivered at any voltage.
     *  Note, this number differs from the capacity of electrodes that is usually quoted for
     *  a battery. That number depends on the rate of discharge and also depends on the
     *  specification of a cutoff voltage. Here, we dispense with both of these specifications.
     *  So, it should be considered a theoretical capacity at zero current and minimal cutoff voltage
     *  considering the current state of the battery. The initial theoretical capacity given
     *  ideal conditions is given by capacityInitial().
     *  It will also include all plateaus that are defined by the electrode object.
     *
     *  @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     *
     *  @return returns the theoretical capacity of the electrode in Amp seconds = coulombs.
     */
    virtual double capacity(int platNum = -1) const override;

    //! Amount of charge that the electrode that has available to be discharged in coulombs
    /*!
     *   We report the number in terms of Amp seconds = coulombs
     *
     *   @param[in]     platNum                  Plateau number. Default is -1 which treats all plateaus as a single entity.
     *                                           If positive or zero, each plateau is treated as a separate entity.
     *   @param[in]     voltsMax                 Maximum voltage to search for capacity. Defaults to +50 volts
     *   @param[in]     voltsMin                 Minimum voltage to search for capacity. Defaults to -50 volts
     *
     *   @return                                 Returns the capacity left  in units of Amp sec = coulombs
     */
    virtual double capacityLeft(int platNum = -1, double voltsMax = 50.0, double voltsMin = -50.0) const override;

    //! Set the relative current capacity discharged per mole
    /*!
     * (virtual function from Electrode)
     * ( this is a serial virtual function - not to be used in creation)
     *
     * This is roughly equal to the total number of electrons that has been discharged
     * from a fully charged state divided by the total moles of solid species in the electrode
     *
     * In this child routine, we translate between relDischarged and RelExtentRxn. The difference between
     * these two concepts is that RelExtentRxn doesn't go between 0 and 1, but has limits based on
     * the relative mole fractions that can be put into the solid. One and only one plateau is also assumed.
     *
     *  @param relDischargedPerMole     Relative value of the discharge per mole. Always goes between 0 and number of electrons
     *                                   per active mole, num
     *                                  0 means that the electrode is fully charged, num means that it is fully discharged.
     *
     *  @param[in]           platNum             Plateau number. Default is -1 which treats all plateaus as a single entity and
     *                                           the relative discharged as a single combined fraction. If platNum is
     *                                           >= 0, then the discharge is relative to the current plateau.
     */
    virtual void setRelativeCapacityDischargedPerMole(double relDischargedPerMole, int platNum = -1) override;

    //! Extract information from reaction mechanisms
    /*!
     *   (inherited from Electrode_Integrator)
     */
    virtual void extractInfo() override;

    //! Calculate the production rate of species in the electrode at the final time of the time step
    /*!
     *   We take the ROP_inner_[] and ROP_outer_[] rates of progress, combine them with the surface area calculation,
     *   and the stoichiometric coefficients to calculate the DspMoles_final_[], which is production rate for
     *   all species in the electrode.
     */
    void updateSpeciesMoleChangeFinal();

    //! Take the state (i.e., the final state) within the Electrode object and push it down
    //! to the ThermoPhase Objects and other variables that are part of the Electrode object
    /*!
     *  We take the values of spMoles_final_[], the number of particles, and their size,
     *  which are the default specification of the state variables within the Electrode object,
     *  and propagate them down to the ThermoPhase objects in the electrode.
     *  We also calculate the volumetric properties of the Electrode, the phase moles,
     *  and the mole fractions.
     *
     *  All of these properties are defined for the _final_ state.
     *
     *
     *
     *  Summary: State Variables
     *              spMoles_final_[kGlobal]  : kGlobal in solid phase species
     *
     *           Independent Variables
     *              spMoles_final_[kGlobal]  :   kGlobal in electrolyte phase
     *              phaseVoltages_[iph] :        iph in solid phase electrodes
     *              Temperature
     *              Pressure
     *
     *   Dependent StateVariables
     *              phaseMoles_final_[iph]       iph in solid phase electrode
     *              spMf_final_[kGlobal]         kGlobal in solid phase species
     *              VolPM_[kGlobal]              kGlobal in solid phase species
     *              spElectroChemPot_[kGlobal]   kGlobal in solid phase species
     *              ThermoPhase[iph]             iph in solid phase electrode
     *        	    phaseMolarVolumes_[iph]      iph in solid phase electrode
     *              ElectrodeSolidVolume_
     *
     */
    virtual void updateState() override;

protected:
    //! Recalculate the surface areas of the surfaces for the final state
    /*!
     *  (virtual function from Electrode)
     *
     *  We used the internal variable locationOfReactingSurface_ to determine the behavior.
     *  A value of zero indicates that the surface 0 follows the reaction front as it goes from outer to inner as
     *  a function of the % though the plateau.
     *
     *  A value of locationOfReactingSurface_ = 1 indicates that the surface 0 follows the exterior surface of the particle
     *
     *  We also assume that the surface area is equal to the particle surface area multiplied by the numbers of particles.
     *
     *  Dependent StateVariables used:
     *
     *         Radius_exterior_final;
     *         particleNumberToFollow_
     *
     *  Dependent StateVariables calculated:
     *          surfaceAreaRS_final_[]
     */
    virtual void updateSurfaceAreas() override;

    //! This is used to set the phase information that is implicit but not set by a restart or an initialization
    /*!
     *  (virtual function from Electrode)
     *
     *  This is called immediately after the restart file's contents are loaded into the electrode object.
     *  We then call this function to calculate the internal flags. Then, we call updateState() to make sure
     *  all information about the state is self-consistent.
     *
     *  We recalculate the values of  onRegionBoundary_final_ and xRegion_final_  based on the current value of 
     *  RelativeExtentRxn_final_ and the values of RelativeExtentRxn_RegionBoundaries_[i]. 
     *  We compare against where we thought we were.
     *
     *  @param[in]           flagErrors          If true any changes in the current flags caused by a mismatch between the state
     *                                           and the values of the flags will cause an error exit.
     *
     *  @return                                  Returns whether there has been any discovered errors (currently ignored)
     */
    virtual bool stateToPhaseFlagsReconciliation(bool flagErrors) override;

public:
    //! Sets the state of the Electrode object given an EState object
    /*!
     *   (virtual function from Electrode)
     *   This sets all of the states within the object to the same state.
     *   It is an error to call this function during a pending step where there can be a difference between t_init and t_final.
     *
     *   @param[in]  es          const reference to the EState object.  Must be the correct EState object for the
     *                           current Electrode object, or else it will throw an error. However, there is an option to 
     *                           read EState objects with less information. 
     */
    virtual void setState_EState(const EState& es) override;

    //! Set the final state of the electrode using the relExtentRxn
    /*!
     *  (virtual from Electrode)
     *  This sets the state of the system, i.e., spmoles_final_[] for the solid phase
     *  components of the electrode using a single number.
     *
     *  It must be the case that  calcRelativeExtentRxn_final() and setState_relativeExtentRxn()
     *  are inverses of one another. Note, this means that if the state of the system has more than one rank,
     *  then the other ranks are unperturbed by the round trip.
     *
     *  @param[in]           relativeExtentRxn   Value of the relative extent of reaction
     */
    virtual void setState_relativeExtentRxn(double relativeExtentRxn) override;

protected:
    //! Predict the solution
    /*!
     *  Ok at this point we have a time step deltalimiTsubcycle_
     *  and initial conditions consisting of phaseMoles_init_ and spMF_init_.
     *  We now calculate predicted solution components from these conditions.
     *
     *  @return                                  Returns the success of the operation
     *                                             1  A predicted solution is achieved
     *                                             2  A predicted solution with a multispecies phase pop is acheived
     *                                             0  A predicted solution is not achieved, but go ahead anyway
     *                                            -1  The predictor suggests that the time step be reduced and a retry occur.
     */
    virtual int predictSoln() override;

    //! Predict the solution at t_final_ from the derivative of the solution at the current time t_init_
    /*!
     *  (virtual from Electrode_Integrator)
     *  Predicts the solution at the final time from the current derivative of the solution at the initial time.
     *  This override handles the cases where a phase is estimated to disappear. 
     *
     *  @return                                  Returns the success of the operation
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
     *                                            y[0] is always the current  deltaTsubcycleCalc_ value
     *
     *  @return                                  Returns the success of the operation
     *                                             1  success
     *                                             0  failure - bounds problems
     */
    virtual int unpackNonlinSolnVector(const double* const ySoln) override;

private:

    //! This routine doesn't look like it checks anything anymore
    /*!
     *  There were some print statements that got commented out
     */
    void checkStillOnRegionBoundary();

    //!  This routine makes sure that the reaction extent is confined to the xRegion_init_ region.
    /*!
     *  If the extent is above the top of the region it's put on that boundary, and onRegionBoundary_final_ is set
     *  If the extent is below the bottom of the region it's put on that boundary, and onRegionBoundary_final_ is set
     *  If the extent is within the region onRegionBoundary_final_ is set to -1.
     */
    void setOnRegionBoundary();

public:

    //! Pack the nonlinear solver proplem
    /*!
     *  formulate the nonlinear solver problem to be solved.
     *     Fields to be filled in
     *             yvalNLS_
     *             ylowNLS_
     *             yhighNLS_
     *             deltaBoundsMagnitudesNLS_
     */
    virtual void initialPackSolver_nonlinFunction() override;

    //! Pack the solution vector
    /*!
     *  Fill up the solution vector with the unknowns for the nonlinear solver for the time step.
     *
     *  @param[in]           y                   solution vector to be filled
     */
    void packNonlinSolnVector(double* const y) const;

    //!  Calculate the norm of the errors in the global source terms
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *   This routine make use of the source term error vector along with rtols and atols for the
     *   individual source terms to calculated a normalized error measure. This is the single number
     *   that the integration routine will try to control as it calculates a time stepping strategy.
     *
     *   @return  Returns a single nondimensional number representing the normalized error
     *            for the calculation of the source term
     */
    virtual double predictorCorrectorGlobalSrcTermErrorNorm() override;

    //! Print table representing prediction vs. corrector information
    /*!
     *  @param yval           Vector of corrector values
     *  @param pnormSrc       Norm of the predictor-corrector comparison for the source vector.
     *  @param pnormSoln      Norm of the predictor-corrector comparison for the solution vector.
     */
    virtual void predictorCorrectorPrint(const std::vector<double>& yval, double pnormSrc, double pnormSoln) const override;

    //! Possibly change the solution due to phase births and deaths.
    /*!
     *   (virtual from Electrode_Integrator)
     *  @return Returns a bool true if the step is acceptable, and false if it is unacceptable.
     */
    virtual bool changeSolnForBirthDeaths() override;

    //! Possibly change the solution due to phase births and deaths after phase has been accepted.
    /*!
     *   (virtual from Electrode_Integrator)
     *
     *  This routine is carried out after the step is deemed a success. Massaging of the solution
     *  must be carried out within strict tolerances.
     */
    virtual void manageBirthDeathSuccessfulStep() override;

    //! Error check on the routine step
    /*!
     *    (virtual from Electrode_Integrator)
     *
     *   Error checks go here. All errors are fatal exits.
     */
    virtual void check_final_state() override;

    //!  Gather the predicted solution values and the predicted integrated source terms
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *  Both the predicted solution values and the predicted integrated source terms are used
     *  in the time step control
     */
    virtual void gatherIntegratedSrcPrediction() override;

    // -----------------------------------------------------------------------------------------------------------------
    // ------------------------------------ CALCULATE INSTANTANEOUS RATES ----------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------
    //     -- These are now called "ProductionRate" values

    //! Calculate the instantaneous time derivative of the species vector as determined by all source terms
    /*!
     *  (virtual from Electrode)
     *
     *  This is the rate of change in the moles of species defined in the electrode  at t_final.
     *  This calculation does not necessarily use an interval of time to calculate anything.
     *
     *  @param spMoleDot   The end result in terms of the rate of change in moles of species in the
     *                     electrode. (kmol s-1)
     */
    virtual void speciesProductionRates(double* const spMoleDot) override;

    // -----------------------------------------------------------------------------------------------------------------------------
    // ----------------------------- CARRY OUT INTEGRATIONS ------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------------------

    //! Reset the Internal state of the electrode at the start of a new global integration step or a redo of the current global step
    /*!
     *  (virtual from Electrode)
     *  This function advances the initial state to the final state that was calculated
     *  in the last integration step. If the initial time is input, then the code doesn't advance or change anything.
     *
     *  @param[in]           Tinitial            This is the New initial time. This time is compared against the "old"
     *                                           final time, to see if there is any problem. They should be the same.
     *                                           If Tinitial == t_init_init, we redo the time step.
     *
     *  @param[in]           doAdvancementAlways Always do the reset, no matter what. Normally, Tinitial is checked against the 
     *                                           current t_init_init value. If they are the same, then we redo the time step.
     *                                           However, if  doAdvancementAlways is true, we advance the solution unknowns to the 
     *                                           final_final values produced in the last global step no matter what.
     *                                           Defaults to false.
     *
     *  @return                                  Returns true if the time step is reset to t_init_init.
     */
    virtual bool resetStartingCondition(double Tinitial, bool doAdvancementAlways = false) override;

    //! Check to see that the preceding step is a successful one
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *   We check to see if the preceding step is a successful one.
     *
     *  @return                                  Returns a bool true if the step is acceptable, and false if it is unacceptable.
     */
    virtual bool  checkSubIntegrationStepAcceptable() const override;

    //! Set the internal initial intermediate and initial global state from the internal final state
    /*!
     *  (virtual function)
     *
     *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
     *  routine as well.
     *
     * @param setInitInit   Boolean indicating whether you should set the init_init state as well
     */
    virtual void setInitStateFromFinal(bool setInitInit = false) override;

    //! Set the internal initial intermediate from the internal initial global state
    /*!
     *  Set the intial state from the init init state. We also can set the final state from this
     *  routine as well.
     *
     *  The final_final is not touched.
     *
     * @param setFinal   Boolean indicating whether you should set the final as well
     */
    virtual void setInitStateFromInitInit(bool setFinal = false) override;

    //! Set the internal initial intermediate and initial global state from the internal final_final state
    /*!
     *  (virtual function)
     *
     *  Set the intial  and init_int state and the final_final from the final state. We also can set the init_init state from this
     *  routine as well.
     *
     */
    virtual void setInitInitStateFromFinalFinal() override;

    //! Set the internal final intermediate and from the internal init state
    /*!
     *  (virtual function from Electrode)
     *
     *  Set the final state from the init state. This is commonly called during a failed time step
     *
     */
    virtual void setFinalStateFromInit() override;

    //! Set the internal final global state from the internal final intermediate state
    /*!
     *  (virtual function from Electrode)
     *
     *  Set the final_final state from the final state. This is commonly called at the end of successful base integration
     */
    virtual void setFinalFinalStateFromFinal() override;

    //  Return a vector of delta y's for calculation of the numerical Jacobian
    /*
     *   There is a default algorithm provided.
     *
     *        delta_y[i] = atol[i] + 1.0E-6 ysoln[i]
     *        delta_y[i] = atol[i] + MAX(1.0E-6 ysoln[i] * 0.01 * solnWeights[i])
     *
     * @param t             Time                    (input)
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input, do not modify)
     * @param delta_y       Value of the delta to be used in calculating the numerical jacobian
     * @param solnWeights   Value of the solution weights that are used in determining convergence (default = 0)
     *
     * @return Returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *            0  Means an unsuccessful operation
     */
    int calcDeltaSolnVariables(const double t, const double* const ySoln,
                               const double* const ySolnDot, double* const deltaYSoln,
                               const double* const solnWeights);

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
    virtual void setNLSGlobalSrcTermTolerances(double rtolResid) override;

    //   Calculate the integrated source terms and do other items now that we have a completed time step
    /*
     *  Calculate source terms on completion of a step. At this point we have solved the nonlinear problem
     *  for the current step, and we are calculating post-processed quantities like source terms.
     */
    virtual void calcSrcTermsOnCompletedStep() override;


    //! Determine the species with the largest mole fraction
    /*!
     *    Fill in the array phaseMFBig_[iph] for all phases in the Electrode.
     *    This is currently called once at the start of the problem, from setResidAtolNLS()
     */
    void determineBigMoleFractions();

    // -----------------------------------------------------------------------------------------------------------------

    //!  Residual calculation for the solution of the Nonlinear time-integration problem
    /*!
     *    (virtual from Electrode)
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
    virtual int integrateResid(const double tfinal, const double delta_t, const double* const y, 
                               const double* const ySolnDot, double* const resid, const ResidEval_Type evalType,
                               const int id_x, const double delta_x) override;

    //! Calculate the residual
    /*!
     *  (virtual function from Electrode_Integrator)
     *  This is the main routine for calculating the residual for the time step. All preliminary calculations have been
     *  carried out, and we are ready to assemble the residual vector.
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
     *  @return                                  Returns a value of 1 if everything went well
     *                                           Returns negative numbers to indicate types of failures
     */
    virtual int calcResid(double* const resid, const ResidEval_Type evalType) override;

    virtual int GFCEO_evalResidNJ(const double t, const double delta_t,
                            const double* const y,
                            const double* const ydot,
                            double* const resid,
                            const ResidEval_Type evalType = ResidEval_Type::Base_ResidEval,
                            const int id_x = -1,
                            const double delta_x = 0.0) override;

    virtual int GFCEO_calcResid(double* const resid, const ResidEval_Type evalType) override;


    void printElectrodeCapacityInfo(int pSrc, bool subTimeStep);

    //! Print condition of a phase in the electrode
    /*!
     *  @param iPhase        Print the phase
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values
     */
    virtual void printElectrodePhase(size_t iPhase, int pSrc = 1,  bool subTimeStep = false) override;

    // ---------------------------------------------------------------------------------------------
    // ---------------------------- SOLUTION OF NONLINEAR TIME DEPENDENT SYSTEM  --------------------
    // ---------------------------------------------------------------------------------------------

    //! Set the Residual and Solution absolute error tolerance vectors
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *   Set the absolute error tolerances fror the nonlinear solvers. This is called at the top
     *   of the integrator() routine.
     *
     *   Calculates atolResidNLS_[]
     *   Calculates atolNLS_[]
     */
    virtual void setResidAtolNLS() override;

    //! Return the number of equations in the equation system that is used to solve the ODE integration
    /*!
     *  (virtual from Electrode_Integrator)
     *  This will be
     *           1   for the deltaT
     *           nSolidPhases   Total moles in each solid phases
     *           Xmol           Mole fraction of species in each solid phase, except for largest mole fraction
     *
     *  @return                                  Returns the number of equations in the ODE integration system
     */
    virtual size_t nEquations_calc() const override;

    //! Evaluate the residual function
    /*!
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
     * @return                                   Returns 1 if everything is ok. Different values mean an error has occurred.
     */
    virtual int evalResidNJ(const double t, const double delta_t,
                            const double* const y,
                            const double* const ydot,
                            double* const resid,
                            const ResidEval_Type evalType = ResidEval_Type::Base_ResidEval,
                            const int id_x = -1,
                            const double delta_x = 0.0) override;

    //! Fill in the initial conditions at the given time value
    /*!
     *  (Virtual from ResidJacEval.h)
     *  Values for both the solution and the value of ydot may be provided.
     *
     *  @param[in]           t0                  Time   
     *  @param[out]          y                   Solution vector (output)
     *  @param[out]          ydot                Rate of change of solution vector. (output)
     *
     *  @return                                  Returns 1 if everything is ok. Error otherwise
     */
    virtual int getInitialConditionsWithDot(const double t0, double* const y, double* const ydot) override;

    //! Calculate the relative extent of reaction from the current state of the object
    /*!
     *  (virtual from Electrode.h)
     *
     *  Calculate the relative extent of reaction from the final state, spmoles_final.
     * 
     *  The relative extent of reaction is a dimensionless number that varies. It doesn't
     *  always vary between 0 and 1. Sometimes there are Li's that can be reacted or sites
     *  that can't be filled with Li.... At 0, the battery is fully charged. At ~1, the battery
     *  is fully discharged.
     *
     *  Here we evaluate this number as the number of electrons produced by the electrode divided by
     *  the normalization constant, which is the initial solid moles in the electrode of active material.
     *
     *  @return returns the relative extent of reaction (dimensionless).
     */
    virtual double calcRelativeExtentRxn_final() const override;

    //! Returns the capacity derivative wrt time in coulombs per second
    /*!
     *  @param[in]           platNum             Plateau number. Defaults to -1, which indicates all plateaus.
     *
     *  @return                                  Returns the derivative of the capacity with respect to time
     *                                           in units of coulombs / sec.
     */
    double capacityDot(int platNum = -1) const;

    //! Returns the capacity left derivative wrt time in coulombs per second
    /*!
     *  @param[in]           platNum             Plateau number. Defaults to -1, which indicates all plateaus.
     *
     *  @return                                  Returns the derivative of the capacity with respect to time
     *                                           in units of coulombs / sec.
     */
    double capacityLeftDot(int platNum = -1) const;

    //! Calculates the derivative of the reaction extent wrt time
    /*!
     *  @param[in]           platNum             Plateau number. Defaults to -1, which indicates all plateaus.
     *
     *  @return                                  Returns the derivative of the reaction extent with respect to time
     *                                           in units of one / sec.
     */
    double RxnExtentDot(int platNum = -1) const;

    //! Returns the capacity left in units of amp sec = coulombs
    /*!
     *  @param[in]           platNum             Plateau number. Defaults to -1, which indicates all plateaus.
     *
     *  @return                                  Returns the capacity left in units of coulombs
     */
    double capacityLeftRaw(int platNum = -1) const;
	
    //! Returns the capacity in units of amp sec = coulombs
    /*!
     *  @param[in]           platNum             Plateau number. Defaults to -1, which indicates all plateaus.
     *
     *  @return                                  Returns the capacity in units of coulombs
     */
    double capacityRaw(int platNum = -1) const; 

    // -------------------------------------       MEMBER DATA --------------------------------------------------------------------


    /* ----------------------------------------------------------------------------------------------
     *             DATA ASSOCIATED WITH REGION BOUNDARIES, RELATIVE EXTENTS, AND THEIR MARKERS
     * ---------------------------------------------------------------------------------------------- */

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
     *
     * @todo duplicate of base class capability -> refactor
     */
    int electrodeType_;

    //! Initial relative extent of reaction at the start of the current step
    double RelativeExtentRxn_init_;

    //! Initial relative extent of reaction at the start of the global step
    double RelativeExtentRxn_init_init_;

    //! Final relative extent of reaction at the end of the current step
    /*!
     *  A note about what relative extent is.
     *
     *    Relative extent is the amount of reaction per mole of reactant. It extends across
     *    multiple plateaus. For graphite, it would correspond to the amount of Lithium per 6
     *    moles of carbon, for example.
     *
     *  In this implementation within Electrode_CSTR, the relative extent is not part of the solution vector. It is calculated
     *  from the mole numbers and follows along with the solution as a tracking equation.
     */
    double RelativeExtentRxn_final_;

    //! final relative extent of reaction at the end of the global step
    double RelativeExtentRxn_final_final_;

    //! Normalization factor that is used to take the Electron production term and turn it into
    //! the relative extent of reaction, which is dimensionless
    /*!
     *  This is usually equal to the total moles of solid phase electrons possible to be created
     *  from the solid phase. We ignore here any issues having to do with boundaries and the
     *  inaccessibility of electrons. 
     *
     *  For anodes like zinc where the solid moles of the electrode decrease as a function of the
     *  extent of reaction, we make the following definition. This value will be equal to the
     *  number of initial moles of zinc, i.e., it will be closely associated with the
     *  value of capacityInitial() / Faraday.
     *
     *  Units: kmol
     */
    double RelativeExtentRxn_NormalizationFactor_;

public:
    //! Boundaries of the regions in terms of the relative extent of reaction
    /*!
     *  Currently, there is some confusion as to whether this must be set or not.
     */
    std::vector<double> RelativeExtentRxn_RegionBoundaries_;

protected:
    //! This integer describes if the system is current on a Region boundary at the start of a subgrid
    //! integration step
    /*!
     *  If it isn't, we set it to -1. If it is, we set it to the region boundary index
     */
    int onRegionBoundary_init_;

    //! This integer describes if the system is current on a Region boundary at the end of a subgrid
    //! integration step
    /*!
     *  If it isn't, we set it to -1. If it is, we set it to the region boundary index
     */
    int onRegionBoundary_final_;

    //! This integer describes if the system is current on a Region boundary at the start of a global
    //! integration step
    /*!
     *  If it isn't, we set it to -1. If it is, we set it to the region boundary index
     */
    int onRegionBoundary_init_init_;

    //! This integer describes if the system is current on a Region boundary at the end of a global integration step
    /*!
     *  If it isn't, we set it to -1. If it is, we set it to the region boundary index
     */
    int onRegionBoundary_final_final_;

    //! Initial region for the extent of reaction
    int xRegion_init_;

    //! Initial region for the extent of reaction of the global step
    int xRegion_init_init_;

    //! Final region for the extent of reaction
    int xRegion_final_;

    //! Final region for the extent of reaction of the global step
    int xRegion_final_final_;

    //!  This Boolean is true when we are at a plateau boundary and the voltage is in between the
    //!  top value and the bottom value.
    int goNowhere_;

    //! Time step to reach a collision with a region boundary.
    /*!
     *  This is used in the nonlinear solver.
     */
    double deltaT_RegionBoundaryCollision_;

    //! Absolute tolerance for nonlinear residual
    double atolBaseResid_;


    // ----------------------------------------------------------------------------------------

    /* ----------------------------------------------------------------------------------------------
     *             DATA ASSOCIATED WITH THE SOLUTION OF NONLINEAR EQUATIONS
     * ---------------------------------------------------------------------------------------------- */

    //! Vector of the rates of progress of the surface reactions
    /*!
     *  Length:   largest value of maxNumRxns in a reacting surface class
     *  Indexing: reaction index for the current ReactingSurface
     *  Units:    kmol / m2 / s
     */
    std::vector<double> ROP_;

    //! Molar source rate for the species vector of all species in the electrode object
    //! for the final time during a time step
    /*!
     *  Units:    kmol s-1
     */
    std::vector<double> DspMoles_final_;

    //! Source for the reaction extent
    /*!
     *  units are kmol s-1
     */
    double SrcDot_RxnExtent_final_;

    //! List of the the volume phases in the PhaseList object which are actually solid phases that are
    //! part of the electrode. 
    /*!
     *  Right now, everything that is not the metal phase or the electrolyte phase is a solid phase.
     */
    std::vector<size_t> phaseIndexSolidPhases_;

    //! The number of species in each of those solid phases that are part of the electrode
    /*!
     *  Length: phaseIndexSolidPhases_.size()
     */
    std::vector<size_t> numSpecInSolidPhases_;

    //! Change in the number of moles of species during the current step
    std::vector<double> deltaSpMoles_;

    //! Phase index which is predicted to die first during the current nonlinear problem
    size_t minPH_;

    //! Vector of species indexes for each phase in the problem which are the largest mole fraction
    //! in that phase
    std::vector<size_t> phaseMFBig_;

    //! Contains a boolean for phases which have died during the current subgrid step.
    /*!
     *  length: m_NumTotPhases
     *  Index:  PhaseList iphGlob
     *  Either has a value of 0 or 1
     */
    std::vector<int> justDied_;

    //! Lagged phase moles vector
    /*!
     *  Values are not updated for jacobian calculations, or for Base_LaggedSolutionComponents calls
     *  Length: Number phases in PhaseList =  m_NumTotPhases
     *  Index:  PhaseList iphGlob
     */
    std::vector<double> phaseMoles_final_lagged_;

    //! Delta phase moles vector
    /*!
     *  Changed in the number of moles in each phase during the local time step
     *  Length: m_NumTotPhases
     *  Index:  PhaseList iphGlob
     */
    std::vector<double> DphMoles_final_;

    //!  Total value of the solid moles at the start of the nonlinear iteration
    /*!
     *  units are kmol
     */
    double solidMoles_init_;

    //!  Minimum size of particle
    /*!
     *   This is used to calculate a minimum surface area for the particle
     */
    double Radius_exterior_min_;

#ifdef DEBUG_ELECTRODE_DEATH
public:
    double time_AtDeath;

    double deltaT_intermed_AtDeath;

    double deltaT_intermed_min_AtDeath;


    double phaseMoles_Init_AtDeath;

    double atol_AtDeath;

    size_t  counterNumberSubIntegrations_atDeath;

#endif
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

