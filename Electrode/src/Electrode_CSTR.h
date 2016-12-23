/*
 * $Id: Electrode_CSTR.h 571 2013-03-26 16:44:21Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_CSTR_H
#define _ELECTRODE_CSTR_H


//#include "cantera/numerics/NonlinearSolver.h"

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
};
//==================================================================================================================================

//! Electrode_CSTR class is an electrode that models a particle as a CSTR
/*!
 *  The base class is close to being a CSTR. This class will have to ensure
 *  that there are no plateaus that create complicated morphologies. There
 *  will just be an interfacial reaction on the exterior of the particle.
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
     * @param right Object to be copied
     */
    Electrode_CSTR(const Electrode_CSTR& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    Electrode_CSTR& operator=(const Electrode_CSTR& right);

    //! Return the type of electrode
    /*!
     *  Returns the enum type of the electrode. This is used in the factory routine.
     *
     *  @return Returns an enum type, called   Electrode_Types_Enum
     */
    virtual Electrode_Types_Enum electrodeType() const;

    //! Set the electrode ID information
    void setID(int domainNum, int cellNum);

    //! This function will extend the input options for the electrode to included child input
    //! and then return a child object of ELECTRODE_KEY_INPUT
    /*!
     *  @param ei   Base input for the electrode
     *              On output it will be an ELECTRODE_CSTR_KEY_INPUT object
     *
     *  @return  Returns zero on success, -1 on failure
     */
    virtual int electrode_input_child(ELECTRODE_KEY_INPUT** ei);

    //!  Setup the electrode using the ELECTRODE_KEY_INPUT object that is read from an input file
    /*!
     *   (virtual from Electrode)
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
     * @param ei    ELECTRODE_KEY_INPUT pointer object
     *
     *  @return  Returns zero if successful, and -1 if not successful.
     */
    virtual int setInitialConditions(ELECTRODE_KEY_INPUT* eibase);

    //! Create an object that saves the electrode state and can print out an XML solution to file
    /*!
     *  The pointer to the malloced object is saved in the internal variable eState_final_ .
     *  Because there is an object, the state of the electrode will be saved at each step.
     *
     *  @return  Returns zero if successful, and -1 if not successful.
     */
    virtual int electrode_stateSave_create();

    //! local routine to resize arrays that this object is responsible for
    void init_sizes();

    //! Set the sizes of the electrode from the input parameters
    /*!
     *
     * @param electrodeArea   Area of the electrode
     * @param electrodeThickness  Width of the electrode
     * @param porosity        Volume of the electrolyte phase
     */
    virtual void setElectrodeSizeParams(double electrodeArea, double electrodeThickness,
                                        double porosity);

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
    virtual void resizeMoleNumbersToGeometry();

public:


    //! Create and malloc the solvers
    /*!
     * @return   returns 1 if ok
     */
    virtual int create_solvers();



    //! Value of the standard state open circuit voltage for the standard state conditions for region xRegion
    //! at the final_ conditions.
    /*!
     *  This is the standard state open circuit potential at the current _final_ conditions.
     *  The value is dependent on the region input. The region input, currently is identified
     *  with the surface reacting phase object.
     *  Therefore, this call is basically a wrapper around openCircuitVoltageSS(isk) which
     *  calculates the open circuit standard state voltage for the isk reacting surface.
     *
     *  Additions include extening the region values to include false values for DoD = 0 and 1
     *  conditions.
     *
     *  @param xRegion  Region   Value of the region. If -1, this is at the DoD = 0. If nR+1,
     *                           this is at the DoD = 1.0 condition
     */
    double openCircuitVoltageSS_Region(int xRegion) const;

    //! Value of the open circuit voltage for region xRegion at the final_ conditions.
    /*!
     *  This is the open circuit potential at the current _final_ conditions.
     *  The value is dependent on the region input. The region input, currently is identified
     *  with the surface reacting phase object.
     *  Therefore, this call is basically a wrapper around openCircuitVoltage(isk) which
     *  calculates the open circuit voltage for the isk reacting surface.
     *
     *  Additions include extending the region values to include false values for DoD = 0 and 1
     *  conditions.
     *
     *  @param xRegion  Region   Value of the region. If -1, this is at the DoD = 0. If nR+1,
     *                           this is at the DoD = 1.0 condition
     */
    double openCircuitVoltage_Region(int xRegion, bool comparetoReferenceElectrode = false) const;

    //! Returns the total capacity of the electrode in Amp seconds
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
    virtual double capacity(int platNum = -1) const;

    //! Amount of charge that the electrode that has available to be discharged
    /*!
     *   We report the number in terms of Amp seconds = coulombs
     *
     *   @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     */
    virtual double capacityLeft(int platNum = -1, double voltsMax = 50.0, double voltsMin = -50.0) const;

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
     *  @param platNum            Plateau number. Default is -1 which treats all plateaus as a single entity and
     *                            the relative discharged as a single combined fraction. If platNum is
     *                            >= 0, then the discharge is relative to the current plateau.
     */
    virtual void setRelativeCapacityDischargedPerMole(double relDischargedPerMole, int platNum = -1);

    //! Extract information from reaction mechanisms
    /*!
     *   (inherited from Electrode_Integrator)
     */
    void extractInfo();


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
    virtual void updateState();

    //!  Recalculate the surface areas of the surfaces for the final state
    /*!
     *    (virtual function from Electrode)
     *
     *    We used the internal variable locationOfReactingSurface_ to determine the behavior.
     *    A value of zero indicates that the surface 0 follows the reaction front as it goes from outer to inner as
     *    a function of the % though the plateau.
     *    A value of locationOfReactingSurface_ = 1 indicates that the surface 0 follows the exterior surface of the particle
     *
     *    We also assume that the surface area is equal to the particle surface area multiplied by the numbers of particles.
     *
     *
     *    Dependent StateVariables Used
     *         Radius_exterior_final;
     *         particleNumberToFollow_
     *
     *    Dependent StateVariables Calculated
     *          surfaceAreaRS_final_[]
     */
    virtual void updateSurfaceAreas();

    //!
    virtual bool stateToPhaseFlagsReconciliation(bool flagErrors);

    //! Set the final state of the electrode using the relExtentRxn
    /*!
     *  This sets the state of the system, i.e., spmoles_final_[] for the solid phase
     *  components of the electrode using a single number.
     *
     *  It is virtual because there is no way to do this except by knowing about the system
     *
     *  It must be the case that  calcRelativeExtentRxn_final() and setState_relativeExtentRxn()
     *  are inverses of one another. Note, this means that if the state of the system has more than one rank,
     *  then the other ranks are unperturbed by the round trip.
     *
     *  @param relExtentRxn  input of the relative extent of reaction
     */
    virtual void setState_relativeExtentRxn(double relativeExtentRxn);

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
    virtual int predictSoln();

    //! Unpack the soln vector
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *  This function unpacks the solution vector into  phaseMoles_final_,  spMoles_final_, and spMf_final_[]
     */
    virtual void unpackNonlinSolnVector(const double* const y) override;
private:
    void checkStillOnRegionBoundary();
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
    virtual void initialPackSolver_nonlinFunction();


    //! Pack the solution vector
    /*!
     *  @param y  solution vector to be filled
     */
    void packNonlinSolnVector(double* const y) const ;

    //!  Calculate the norm of the difference between the predicted answer and the final converged answer
    //!  for the current time step
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *   The norm calculated by this routine is used to determine whether the time step is accurate enough.
     *
     *  @return    Returns the norm of the difference. Normally this is the L2 norm of the difference
     */
    virtual double predictorCorrectorWeightedSolnNorm(const std::vector<double>& yvalNLS);

    //! Calculate the vector of predicted errors in the source terms that this integrator is responsible for
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *    In the base implementation we assume that the there are just one source term, the electron
     *    source term.
     *    However, this will be wrong in almost all cases.
     *    The number of source terms is unrelated to the number of unknowns in the nonlinear problem.
     *    Source terms will have units associated with them.
     *    For example the integrated source term for electrons will have units of kmol
     */
    virtual void predictorCorrectorGlobalSrcTermErrorVector();

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
    virtual double predictorCorrectorGlobalSrcTermErrorNorm();

    //! Print table representing prediction vs. corrector information
    /*!
     *  @param yval           Vector of corrector values
     *  @param pnormSrc       Norm of the predictor-corrector comparison for the source vector.
     *  @param pnormSoln      Norm of the predictor-corrector comparison for the solution vector.
     */
    virtual void predictorCorrectorPrint(const std::vector<double>& yval,
                                         double pnormSrc, double pnormSoln) const;

    //! Possibly change the solution due to phase births and deaths.
    /*!
     *   (virtual from Electrode_Integrator)
     *  @return Returns a bool true if the step is acceptable, and false if it is unacceptable.
     */
    virtual bool changeSolnForBirthDeaths();

    //! Possibly change the solution due to phase births and deaths after phase has been accepted.
    /*!
     *   (virtual from Electrode_Integrator)
     *
     *  This routine is carried out after the step is deemed a success. Massaging of the solution
     *  must be carried out within strict tolerances.
     */
    virtual void manageBirthDeathSuccessfulStep();

    //! Error check on the routine step
    /*!
     *    (virtual from Electrode_Integrator)
     *
     *   Error checks go here. All errors are fatal exits.
     */
    virtual void check_final_state();

    //!  Gather the predicted solution values and the predicted integrated source terms
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *  Both the predicted solution values and the predicted integrated source terms are used
     *  in the time step control
     */
    virtual void gatherIntegratedSrcPrediction();

    double l0normM(const std::vector<double>& v1, const std::vector<double>& v2, int num,
                   const std::vector<double>& atolVec, const double rtol) const;


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
    virtual void speciesProductionRates(double* const spMoleDot);

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
     *  @param[in]           doResetAlways       Always do the reset, no matter what. Normally, Tinitial is checked against the 
     *                                           current t_init_init value. If they are the same, then we redo the time step.
     *                                           However, if  doResetAlways is true, we advance the solution unknowns to the 
     *                                           final_final values produced in the last global step no matter what.
     *                                           Defaults to false.
     */
    void  resetStartingCondition(double Tinitial, bool doTestsAlways = false) override;

    //! Check to see that the preceding step is a successful one
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *   We check to see if the preceding step is a successful one.
     *
     *  @return Returns a bool true if the step is acceptable, and false if it is unacceptable.
     */
    virtual bool  checkSubIntegrationStepAcceptable() const;

    //! Set the internal initial intermediate and initial global state from the internal final state
    /*!
     *  (non-virtual function onionize in-first)
     *
     *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
     *  routine as well.
     *
     * @param setInitInit   Boolean indicating whether you should set the init_init state as well
     */
    void setInitStateFromFinal_Oin(bool setInitInit = false);

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

    //! Set the internal initial intermediate and initial global state from the internal final_final state
    /*!
     *  (virtual function)
     *
     *  Set the intial  and init_int state and the final_final from the final state. We also can set the init_init state from this
     *  routine as well.
     *
     */
    virtual void setInitInitStateFromFinalFinal();

    //! Set the internal final intermediate and from the internal init state
    /*!
     *  (virtual function from Electrode)
     *
     *  Set the final state from the init state. This is commonly called during a failed time step
     *
     */
    virtual void setFinalStateFromInit();

    //! Set the internal final global state from the internal final intermediate state
    /*!
     *  (virtual function from Electrode)
     *
     *  Set the final_final state from the final state. This is commonly called at the end of successful base integration
     */
    virtual void setFinalFinalStateFromFinal();

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
    virtual void setNLSGlobalSrcTermTolerances(double rtolResid);

    //   Calculate the integrated source terms and do other items now that we have a completed time step
    /*
     *  Calculate source terms on completion of a step. At this point we have solved the nonlinear problem
     *  for the current step, and we are calculating post-processed quantities like source terms.
     */
    virtual void calcSrcTermsOnCompletedStep();


    //! Determine the species with the largest mole fraction
    /*!
     *    Fill in the array phaseMFBig_[iph] for all phases in the Electrode.
     *    This is currently called once at the start of the problem, from setResidAtolNLS()
     */
    void determineBigMoleFractions();

    // -----------------------------------------------------------------------------------------------------------------

    //!  Residual calculation for the solution of the Nonlinear integration problem
    /*!
     *    Given tfinal and delta_t, and given y and ydot which are estimates of the solution
     *    and solution derivative at tfinal, this function calculates the residual equations.
     *    It is the residual function used in the nonlinear solver that relaxes the equations
     *    at each time step.
     *
     *    This is typically called from evalResidNJ(), which is called directly from the
     *    nonlinear solver. However, we expose this routine so that the residual can be queried
     *    given all of the inputs.
     *
     * @param[in]  tfinal      Time                    (input)
     * @param[in]  delta_t     The current value of the time step (input)
     * @param[in]  y           Solution vector (input, do not modify)
     * @param[in]  ydot        Rate of change of solution vector. (input, do not modify)
     * @param[out] resid       Value of the residual that is computed (output)
     * @param[in]  evalType    Type of the residual being computed (defaults to Base_ResidEval)
     * @param[in]  id_x        Index of the variable that is being numerically differenced to find
     *                         the jacobian (defaults to -1, which indicates that no variable is being
     *                         differenced or that the residual doesn't take this issue into account)
     * @param[in]  delta_x     Value of the delta used in the numerical differencing
     *
     * @return                Returns an integer that gets fed back through evalResidNJ() to the
     *                        nonlinear solver. Anything other than a 1 causes an immediate failure
     *                        of the nonlinear solver to occur.
     */
    int integrateResid(const double tfinal, const double delta_t,
                       const double* const y, const double* const ySolnDot,
                       double* const resid,
                       const ResidEval_Type_Enum evalType, const int id_x,
                       const double delta_x);

    virtual int calcResid(double* const resid, const ResidEval_Type_Enum evalType);

    virtual int GFCEO_evalResidNJ(const double t, const double delta_t,
                            const double* const y,
                            const double* const ydot,
                            double* const resid,
                            const ResidEval_Type_Enum evalType = Base_ResidEval,
                            const int id_x = -1,
                            const double delta_x = 0.0);

    virtual int GFCEO_calcResid(double* const resid, const ResidEval_Type_Enum evalType);


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
                            const ResidEval_Type_Enum evalType = Base_ResidEval,
                            const int id_x = -1,
                            const double delta_x = 0.0) override;

    //! Fill in the initial conditions at the given time value
    /*!
     *  (Virtual from ResidJacEval.h)
     * Values for both the solution and the value of ydot may be provided.
     *
     * @param t0            Time                    (input)
     * @param y             Solution vector (output)
     * @param ydot          Rate of change of solution vector. (output)
     *
     * @return                                   Returns 1 if everything is ok. Error otherwise
     */
    virtual int getInitialConditions(const double t0, double* const y, double* const ydot) override;

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

    double capacityDot(int platNum = -1) const;
    double capacityLeftDot(int platNum = -1, double voltsMax = 50.0, double voltsMin = - 50.0) const;
    double RxnExtentDot(int platNum = -1, double voltsMax = 50.0, double voltsMin = -50.0) const;
    double capacityLeftRaw(int platNum = -1, double voltsMax = 50.0, double voltsMin = -50.0) const;
	
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
     *    moles of carbon.
     *
     *  In this implementation, the relative extent is not part of the solution vector. It is calculated
     *  from the mole numbers and follows along with the solution.
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
     *  units = kmol
     */
    double RelativeExtentRxn_NormalizationFactor_;

public:
    //! Boundaries of the regions in terms of the relative extent of reaction
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

    //! This integer describes if the system is current on a Region boundary at the end of a global
    //! integration step
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

    std::vector<double> ROP_;

    //! Molar source rate for the species vector of all species in the electrode object
    //! for the final time during a time step
    /*!
     *  units kmol s-1
     */
    std::vector<double> DspMoles_final_;

    //! Source for the reaction extent
    /*!
     *  units are kmol s-1
     */
    double SrcDot_RxnExtent_final_;

    //! List of the the volume phases in the PhaseList object  which are actually solid phases that are
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
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

