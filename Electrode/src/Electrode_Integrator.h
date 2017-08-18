/**
 *  @file Electrode_Integrator.h
 *     Headers for the declarations of the Electrode_Integrator class, used to perform
 *     subtimestep integrations on top of the Electrode class
 *     (see \ref electrode_mgr and class \link Zuzax::Electrode_Integrator Electrode_Integrator\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_INTEGRATOR_H
#define _ELECTRODE_INTEGRATOR_H

#include "Electrode.h"
#include "cantera/numerics/ResidJacEval.h"
//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
class NonlinearSolver;
class SquareMatrix;
//==================================================================================================================================
//! This class stores an overall summary of one time step taken in the integrator
/*!
 *  It contains the starting and ending times. It contains the number of nonlinear iterations,
 *  and a couple of key results like the voltage and the source term for electrons.
 *  The equals operator is overloaded for comparison of one of these objects with another
 */
class TimeStepHistory
{
public:

    //! Default constructor
    TimeStepHistory() :
        t_init_(0.0),
        t_final_(0.0),
        t_final_calc_(0.0),
	timeTypeSoln_(0),
        numNonLinSolves_(0),
        solnErrorNorm_(0.0),
        wdotErrorNorm_(0.0),
	volts_(0.0),
	srcTermStepElectrons_(0.0),
	currentStep_(0.0)
    {
    }

    //! Constructor that is used to initialize all of the fields
    /*!
     *    @param[in]         t_init              Initial time of the step
     *    @param[in]         t_final             Final time of the calc called for
     *    @param[in]         t_final_calc        Final time of calc actually computed
     *    @param[in]         timeTypeSoln        type of the solution (did we run into a discontinuous point or not)
     *    @param[in]         numNonLinSolves     Number of nonlinear solves
     *    @param[in]         solnErrorNorm       Error norm of the solution during integration
     *    @param[in]         wdotErrorNorm       Error norm on the integrated source terms that are the final product
     *                                           of this electrode object
     */ 
    TimeStepHistory(double t_init, double t_final, double t_final_calc, int timeTypeSoln, int numNonLinSolves,
		    double solnErrorNorm, double  wdotErrorNorm) :
	t_init_(t_init),
        t_final_(t_final),
        t_final_calc_(t_final_calc),
	timeTypeSoln_(timeTypeSoln),
        numNonLinSolves_(numNonLinSolves),
        solnErrorNorm_(solnErrorNorm),
        wdotErrorNorm_(wdotErrorNorm),
	volts_(0.0),
	srcTermStepElectrons_(0.0),
	currentStep_(0.0)
    {
    }

    //! Copy constructor
    /*!
     *  @param[in]           right               object to be copied
     */ 
    TimeStepHistory(const TimeStepHistory& right) :
        t_init_(right.t_init_),
        t_final_(right.t_final_),
        t_final_calc_(right.t_final_calc_),
        timeTypeSoln_(right.timeTypeSoln_),
	numNonLinSolves_(right.numNonLinSolves_),
        solnErrorNorm_(right.solnErrorNorm_),
        wdotErrorNorm_(right.wdotErrorNorm_) ,
	volts_(right.volts_),
	srcTermStepElectrons_(right.srcTermStepElectrons_),
	currentStep_(right.currentStep_)
    {
    }

    //! Assignment operator
    /*!
     *  @param[in]           right               object to be copied
     *  @return                                  Returns a reference to the current object
     */ 
    TimeStepHistory& operator=(const TimeStepHistory& right) 
    { 
	if (this == &right) {
	    return *this;
	}
	t_init_             = right.t_init_;
	t_final_            = right.t_final_;
	t_final_calc_       = right.t_final_calc_;
	numNonLinSolves_    = right.numNonLinSolves_;
	timeTypeSoln_       = right.timeTypeSoln_;
	solnErrorNorm_      = right.solnErrorNorm_;
	wdotErrorNorm_      = right.wdotErrorNorm_;
	volts_              = right.volts_;
	srcTermStepElectrons_ = right.srcTermStepElectrons_;
	currentStep_        = right.currentStep_;
	return *this;
    }

    //!  Overload of the "equals to" operator for the TimeStepHistory object
    /*!
     *  Compares the current time step history with a different object. 
     *  A comparison of times is made with an absolute tolerance of 1.0E-13 (hard coded)
     *  The number of nonlinear solves is also checked.
     *  
     *  @param[in]           other               Reference to the TimeStepHistory object that will be compared against
     *  @return                                  Returns true if the two objects have the equivalent time step history
     */
    bool operator==(const TimeStepHistory& other) const
    {
	if (fabs(t_init_ - other.t_init_) > 1.0E-13) {
	    return false;
	}
	if (fabs(t_final_ - other.t_final_) > 1.0E-13) {
	    if (timeTypeSoln_ != 2) {
		return false;
	    }
	}
	if (fabs(t_final_calc_ - other.t_final_calc_) > 1.0E-13) {
	    return false;
	}
	if (numNonLinSolves_ != other.numNonLinSolves_) {
	    return false;
	}
	if (timeTypeSoln_ != other.timeTypeSoln_) {
	    return false;
	}

	return true;
    }

    //! Overload of the "not equals to" operator for the TimeStepHistory object
    /*!
     *  Compares the current time step history with a different object. 
     *  
     *  @param[in]           other               Reference to the TimeStepHistory object that will be compared against
     *  @return                                  Returns true if the two objects don't have the equivalent time step history
     */
    bool operator!=(const TimeStepHistory& other) const
    {
	return !(*this == other);
    }

    // ------------------------------------------------ D A T A -------------------------------------------------------

    //! Initial time step
    double t_init_;

    //! Final time step before special events
    double t_final_;

    //! Actually calculated final time
    /*!
     *  Can be different if there was a special event such as a phase death
     */
    double t_final_calc_;

    //!  Type of the solution
    /*!
     *      0 normal solution type where the time step is determined from
     *        a predictor-corrector step
     *      1 normal solution type where the time step is determined from
     *        a fixed discretization of the global time step
     *      2 Special time step where the time step is determined by
     *        a special event causing loss of source term continuity (phase death or phase birth)
     */
    int timeTypeSoln_;

    //! Number of nonlinear iterations
    int numNonLinSolves_;

    //! Error norm for that step
    double solnErrorNorm_;

    //! Error norm for the source term for that step
    double wdotErrorNorm_;

    //! Volts
    /*!
     *  Volts assumed during the step (t_final)
     */
    double volts_;

    //! SrcTermElectrons
    /*!
     *   source term for the electrons ( kmol s-1)
     */
    double srcTermStepElectrons_;

    //! Value of the current on this particular step
    /*!
     *   This will be equal to the srcTerm divided by the step length
     */
    double currentStep_;
};
//===================================================================================================================================
//! Complete time step history for a subtime step integration 
/*!
 *  SubIntegrations are multiple time steps taken by the Electrode object per global time step
 *  This object stores all of the information, which consists of multiple TimeStepHistory objects along with some header information.
 *
 *  We also use this class to repete time step histories for calculating more accurate jacobians.
 *  We keep a counter iCounter, which is the current counter. Then, we ask for the next time step step 
 *  using getNextRegularTimeStep() so that we can repeat the time step history during a numerical jacobian calculation.
 */
class SubIntegrationHistory
{
public:
   //! Subintration history default constructor
    SubIntegrationHistory();

    //! Copy constructor
    /*!
     *  @param[in]           right              Object to be copied
     */
    SubIntegrationHistory(const SubIntegrationHistory& right);

    //! Default destructor
    ~SubIntegrationHistory();


    //! Assignment operator
    /*!
     *  @param[in]           right              Object to be copied
     *  @return                                 Returns a reference to the current object
     */
    SubIntegrationHistory& operator=(const SubIntegrationHistory& right);

    //! Wipe out the history and all of the steps
    void clear();

    //! Add a step to the time history of the subintegration
    /*!
     *    @param[in]         t_init             Initial time of the step
     *    @param[in]         t_final             Final time of the calc called for
     *    @param[in]         t_final_calc        Final time of calc actually computed
     *    @param[in]         timeTypeSoln        type of the solution (did we run into a discontinuous point or not)
     *    @param[in]         numNonLinSolves     Number of nonlinear solves
     *    @param[in]         solnErrorNorm       Error norm of the solution during integration
     *    @param[in]         wdotErrorNorm       Error norm on the integrated source terms that are the final product
     *                                           of this electrode object
     *    @param[in]         volts               Returns the volts at the end of the step (volts)
     *    @param[in]         srcElectronsStep    Source term for electrons over the step (kmol)
     *    @param[in]         currentStep         Current during the step    (amps)
     */ 
    void addTimeStep(double t_init, double t_final, double t_final_calc, int timeTypeSoln, 
		     int numNonLinSolves, double solnErrorNorm, double wdotErrorNorm, double volts,
	             double srcElectronsStep, double currentStep);

    //!  zero out the time step counter.
    void zeroTimeStepCounter();

    //! Advance iCounter if it's below nTimeSteps_
    void advanceTimeStepCounter();

    // Return the next regular time that a time step ended given the value of iCounter
    /*!
     *  If we ran out of times, use the stored value of time_step_next to create a good deltaT.
     *  We skip over special time steps that had discontinuities, because their times were calculated. 
     *  And, the current calculation will just calculate a slightly different time.
     *
     *  (this is a messy routine with lots of logic)
     *
     *  @param[in]           currentTime         Current time of the new calculation
     *  @return                                  Returns the next time to integrate to.
     */
    double getNextRegularTime(double currentTime) const;


    //! This routine does some error checking on the times
    /*!
     *  @param[in]           gtinit              value of the global initial time
     *  @param[in]           gtfinal             value of the global final time
     *
     *  @return                                  returns the final counter value
     */
    int assureTimeInterval(double gtinit, double gtfinal);

    //! Supply a history for the integrator to follow
    /*!
     *  This creates a constant time step history 
     *
     *   @param tinit Starting global time step
     *   @param tfinal     Final global time step
     *   @param nsteps     number of local regular time steps
     */
    void setConstantStepSizeHistory(double tinit, double tfinal, int nsteps);

    //! Returns the global start time
    /*!
     *  @return                                  Returns the global start time
     */
    double globalStartTime() const;

    //! Returns the global end time
    /*!
     *  @return                                  Returns the global end time
     */
    double globalEndTime() const;

    //! Prints out the contents of the object to stdout
    /*!
     *  @param[in]           lvl                 Print lvl
     *                                           0  - one line
     *                                           1  - table of values
     *                                           2  - extensive printing
     */
    void print(int lvl) const;

    //! Overload of the "==" operator compares two histories to see if they are equivalent
    /*!
     *  @param[in]           other               Reference to the other SubIntegrationHistory object to compare against
     *
     *  @return                                  Returns true if equivalent
     */
    bool operator==(const SubIntegrationHistory& other) const;

    //! Overload of the "!=" operator compares two histories to see if they are different
    /*!
     *  @param[in]           other               Reference to the other SubIntegrationHistory object to compare against
     *
     *  @return                                  Returns true if they are different
     */
    bool operator!=(const SubIntegrationHistory& other) const;
	
    //! Number of regular time steps, where there are no special events 
    int nTimeStepsRegular_;

    //! Number of time steps, no matter what the kind.
    int nTimeSteps_;

    //! Vector of stored TimeStepHistory objects for the global time step
    /*!
     *  Length:  nTimeSteps_
     */
    std::vector<TimeStepHistory> TimeStepList_;

    //! Global solution error norm value for the global step
    double GsolnErrorNorm_;

    //! Global solution source term error norm value for the global step
    double GwdotErrorNorm_;

    //! Counter for replaying the integration history
    mutable int iCounter;

    //! Internal storage for the next time step amount to be taken after the current step
    double time_step_next;
};
//==================================================================================================================================

//! This class is a derived class used to model phase - change electrodes
/*!
 *  It inherits from Electrode and from ResidJacEval.
 *
 *  Therefore, it can be called by any DAE integrator as a residual function evaluator.
 */
class Electrode_Integrator : public Electrode , public ZZCantera::ResidJacEval
{
public:

    // ---------------------------------------------------------------------------------------------
    // ----------------------- BASIC SETUP ROUTINES  -----------------------------------------------
    // ---------------------------------------------------------------------------------------------

    //! Constructor
    Electrode_Integrator();

    //! Destructor
    virtual  ~Electrode_Integrator();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    Electrode_Integrator(const Electrode_Integrator& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     *
     *  @return                                  Returns a reference to the current object
     */
    Electrode_Integrator& operator=(const Electrode_Integrator& right);

    //!  Create the electrode model from the key input
    /*!
     *   (virtual from Electrode)
     *   @param[in]          ei                  Pointer to the ELECTRODE_KEY_INPUT containing the parsed input for the model
     *
     *   @return                                 Returns 0 if successful and -1 if not successful
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
     *  @param               ei                  ELECTRODE_KEY_INPUT pointer object
     *
     *  @return                                  Returns zero if successful, and -1 if not successful.
     */
    virtual int setInitialConditions(ELECTRODE_KEY_INPUT* ei) override;

    //! Calculate the number of equations that will be solved during the nonlinear solver step.
    /*!
     *  (virtual from Electrode_Integrator)
     *  All classes which inherit from this routine must have a class that determines this value.
     *  The parent class, this function, will throw an error.
     *
     *  @return                                  Returns the number of unknowns in the nonlinear problem and time-stepping problem.
     */
    virtual size_t nEquations_calc() const;

    //! Create and malloc the solvers used to solve the nonlinear problem
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *   As a prerequisite before calling, this function must have a valid count of the number of unknowns in the nonlinear problem
     *   when we call this. If we don't, then the solver arrays won't be malloced with the right amount of space.
     *   This means that nEquations() must return a good value.
     *   This means that the function must be called after or at least at the end of electrode_model_create().
     *
     *  @return                                  Returns 1 if ok
     */
    virtual int create_solvers();

    // ---------------------------------------------------------------------------------------------
    // ----------------------- SPECIFY AND OBTAIN PROBLEM PARAMETERS -------------------------------
    // ---------------------------------------------------------------------------------------------

    //! Setup the vectors for handling the error control on source terms
    /*!
     *  This routine is called at setup time.
     *
     *  The default behavior is to designate the entire vector, spMoleIntegratedSourceTerm_,
     *  as the source term whose error will be controlled.
     *
     *  This routine is responsible for calculating numIntegratedSrc_.
     *
     *  @return Returns the number of source terms whose error will be controlled.
     */
    virtual int setupIntegratedSourceTermErrorControl();

    // -----------------------------------------------------------------------------------------------------------------------------
    // ----------------------------- CARRY OUT INTEGRATION OF EQUATIONS-------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------------------

    //! Reset the Internal state of the electrode at the start of a new global integration step or a redo of the current global step
    /*!
     *  (virtual from Electrode)
     *  This function advances the initial state to the final state that was calculated in the last integration step.
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
    virtual void resetStartingCondition(double Tinitial, bool doResetAlways = false) override;

    //!  Calculate the change in the state of the system when integrating from T_initial_initial_  to t_final_final_
    /*!
     *  All information is kept internal within this routine. This may be done continuously and the solution is not updated.
     *
     *  Note the tolerance parameters refere to the nonlinear solves within the calculation
     *  They do not refer to time step parameters.
     *
     *  @param deltaT        DeltaT for the integration step.
     *  @param GlobalRtolSrcTerm Relative tolerance for the source term vector calcualted from
     *                       this routine.
     *                       Defaults to 1.0E-3
     *  @param fieldInterpolationType Type of interpolation of field variables defaults to T_FINAL_CONST_FIS,
     *  @param subIntegrationType     Type of subintegration. Defaults to BASE_TIMEINTEGRATION_SIR.
     *                                In this integration, the program determines its own strategy
     *                                for the time step.
     *
     *  @return Returns the number of subcycle steps it took to complete the full step.
     *          Failures to complete the integration due to time truncation error issues return a -1.
     *          Failures due to invalid function calculation attempts return a -2.
     *          Failures due to invalid arguments return a -3.
     */
    virtual int integrate(double deltaT, double  GlobalRtolSrcTerm = 1.0E-3,
                          Electrode_Exterior_Field_Interpolation_Scheme_Enum fieldInterpolationType = T_FINAL_CONST_FIS,
                          Subgrid_Integration_RunType_Enum subIntegrationType = BASE_TIMEINTEGRATION_SIR) override;

    //! Calculate the largest mole fraction in each of the phases
    /*!
     *  We use this to determine the equation system
     */
    virtual void determineBigMoleFractions();

    //! Check the consistency of the initial state
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *   This is a checker routine. Therefore, it can be taken out in nondebug mode
     */
    virtual void check_init_consistency() const;

    //!   Prepare the problem for predicting the solution
    /*!
     *      -> determine the largest mole fraction in a phase
     *      -> determine if a surface phase can have reactions turned on
     *      -> do stuff that affects the equation system that will be solved at the
     *         current step.
     */
    virtual void prepareProblemStatement();

    //! Predict the solution
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *  Ok at this point we have a time step and initial conditions consisting of phaseMoles_init_ and spMF_init_.
     *  We now calculate predicted solution components from these conditions.
     *  Additionally, we evaluate whether any multispecies phases are going to pop into existence,
     *  adding in potential seed values, or die, creating a different equation set and time step equation.
     *  We set the phaseExistence flags in the kinetics solver to reflect phase pops.
     *
     *  @return                                  Returns the success of the operation
     *                                           -  1  A predicted solution is achieved
     *                                           -  2  A predicted solution with a multispecies phase pop is acheived
     *                                           -  0  A predicted solution is not achieved, but go ahead anyway
     *                                           - -1  The predictor suggests that the time step be reduced and a retry occur.
     */
    virtual int predictSoln();

    //! Predict the solution at t_final_ from the derivative of the solution at the current time t_init_
    /*!
     * (virtual from Electrode_Integrator)
     * Predicts the solution at the final time from the current derivative of the solution at the initial time.
     *
     * @return                                   Returns the success of the operation
     *                                           -  1  A predicted solution is achieved
     *                                           -  2  A predicted solution with a multispecies phase pop is acheived
     *                                           -  0  A predicted solution is not achieved, but go ahead anyway
     *                                           - -1  The predictor suggests that the time step be reduced and a retry occur.
     */
    virtual int predictSolnDot();

    //! Extract the ROP of the multiple reaction fronts from Cantera within this routine
    /*!
     *  (virtual function from Electrode_Integrator)
     *
     *  In this routine we calculate the rates of progress of reactions and species on all active reacting surfaces.
     *  The function loops over active ReactingSurDomain objects. It calculates ROP_ for each reaction by calling                   
     *                         getNetRatesOfProgress()
     * 
     *  It gets spNetProdPerArea by calling:
     *            getNetSurfaceProductionRates(isk, spNetProdPerArea);
     *   
     *  It may also fill in 
     *            justBornPhase_[]
     */
    virtual void extractInfo();

    //! Collect mole change information
    /*!
     *   (inherited from Electrode_Integrator)
     *   We take the ROP_inner_[] and ROP_outer_[] rates of progress, combine them with the surface area calculation,
     *   and the stoichiometric coefficients to calculate the DspMoles_final_[], which is production rate for
     *   all species in the electrode due to surface reactions
     */
    virtual void updateSpeciesMoleChangeFinal();

    //! Pack the nonlinear solver proplem
    /*!
     *  Formulate the nonlinear solver problem to be solved.
     *     Fields to be filled in
     *             yvalNLS_
     *             ylowNLS_
     *             yhighNLS_
     *             deltaBoundsMagnitudesNLS_
     */
    virtual void initialPackSolver_nonlinFunction();

    //! Formulate and/or check the value of yvalNLS_init_
    /*!
     *  Create the initial value of the nonlinear solution vector at the init time, yvalNLS_init_[].
     *
     *  @param[in]           doOthers            If true copy this solution vector into the the times yvalNLS_init_init_[],
     *                                           yvalNLS_[]_ and yvalNLS_final_final_[]
     */
    virtual void check_yvalNLS_init(bool doOthers);

    //! Check the nonlinear residual equations for completeness and the ability to be solved
    /*!
     *  @return                                   0 Everything is good
     *                                           -1 residual isn't good. We need to cut the time step and retry again.
     */
    virtual int check_nonlinResidConditions();


    // ----------------------------------------  ERROR ANALYSIS OF THE INTEGRATION ------------------------------------

    //! Report the number of state variables and their relative integration errors during the
    //! current global integration step
    /*!
     *   Note rtol doesn't factor into this immediately. Therefore, a value or 1E-3
     *                                  would mean the error in the value is 1 part in 1000.
     *
     *  @param[out] numSV               Returns the number of state variables
     *  @param[out] errorVector         Returns a vector of errors in the state variables for the global step
     *                                  Note rtol doesn't factor into this immediately. Therefore, a value or 1E-3
     *                                  would mean the error in the value is 1 part in 1000.
     *  @return     Returns the large value of the errors in the errorVector.
     */
    double reportStateVariableIntegrationError(int& numSV, double* const errorVector) const;


    //---------------------------------------------------------------------------------------------
    // ---------------------------- INTEGRATED SOURCE TERM CALCULATIONS ------- ----------------------------------------------
    // ---------------------------------------------------------------------------------------------

    //!  Returns the number of integrated source terms whose errors are controlled by this integrator
    /*!
     *   @return                                 Returns the number of integrated source terms
     */
    size_t numIntegratedSrcTerms() const;

    //! Gather the predicted solution values and the predicted integrated source terms
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *  Both the predicted solution values and the predicted integrated source terms are used
     *  in the time step control
     */
    virtual void gatherIntegratedSrcPrediction();


    //! Calculate the integrated source terms and do other items now that we have a completed time step
    /*!
     *  Calculate source terms on completion of a step. At this point we have solved the nonlinear problem
     *  for the current step, and we are calculating post-processed quantities like source terms.
     */
    virtual void calcSrcTermsOnCompletedStep();

    //! Accumulate src terms and other results from the local step into the global holding bins.
    /*!
     *  Accumulate source terms on completion of a step. At this point we have solved the nonlinear problem
     *  for the current step and we have satisfied all accuracy requirements.
     *  The step is good. We now accumulate the results before going on to a new local step.
     *  For example, for the species source terms, we add the value of spMoleIntegratedSourceTermLast_[i], which is the source term
     *  from the last sub-global interval time step to the accumulation vector for the global time step,
     *  spMoleIntegratedSourceTerm_[i] in order to keep track of the global time step terms.
     *
     *  @param[in]           remove              If true we remove the contributions from the last sub-time step
     *                                           from the sum of the previous sub-time steps. default is false
     */
    virtual void accumulateSrcTermsOnCompletedStep(bool remove = false);

    //! Set the base tolerances for the nonlinear solver within the integrator
    /*!
     *  The tolerances are based on controlling the integrated global source term (e.g., electron source term)
     *  for the electrode over the integration interval.  The integrated source term has units of kmol.
     *
     *  Because the electron is only one molar quantity within a bunch of molar quantities,
     *  this requirement will entail that we control the source terms of all species within the
     *  electrode to the tolerance requirements of the electron source term.
     *
     *  @param[in]           rtolResid           Relative tolerance allowed for the global source term over the interval.
     *                                           This is a unitless quantity
     */
    virtual void setNLSGlobalSrcTermTolerances(double rtolResid);

    //! Zero vectors that are accumulated over local time step to represent global time step quantities
    /*!
     *  Zero the integrated source term contributions such as spMoleIntegratedSourceTerm_. We do this when
     *  we are starting a new global time step or reintegrating the current global time step.
     *  Estimated global time step errors.
     */
    virtual void zeroGlobalStepAccumulationTerms();

    // ---------------------------------------------------------------------------------------------
    // ---------------------------- SOLUTION OF NONLINEAR TIME DEPENDENT SYSTEM  ----------------------------------------------
    // ---------------------------------------------------------------------------------------------

protected:
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
    virtual void setResidAtolNLS();

public:
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
     * @return                                   Returns a value of 1 if everything went well
     *                                           Returns negative numbers to indicate types of failures
     */
    virtual int evalResidNJ(const double t, const double delta_t, const double* const y,
                            const double* const ydot, double* const resid,
                            const ResidEval_Type evalType = ResidEval_Type::Base_ResidEval,
                            const int id_x = -1, const double delta_x = 0.0) override;

    //! Calculate the residual
    /*!
     *  (virtual function from Electrode_Integrator)
     *
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
    virtual int calcResid(double* const resid, const ResidEval_Type evalType);

    //! Evaluate the residual function for the "Globally Fully Coupled Electrode Object" formulation
    /*!
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
     * @return                                   Returns a value of 1 if everything went well
     *                                           Returns negative numbers to indicate types of failures
     */
    virtual int GFCEO_evalResidNJ(const double t, const double delta_t,
                            const double* const y,
                            const double* const ydot,
                            double* const resid,
                            const ResidEval_Type evalType = ResidEval_Type::Base_ResidEval,
                            const int id_x = -1,
                            const double delta_x = 0.0);

    //! Calculate the residual for the Electrode object for the Global problem
    /*!
     *  @param[out]          resid               Value of the residual
     *  @param[in]           evalType            ResidEval_Type_Enum type
     *
     *  @return                                  1 Means a good calculation that produces a valid result
     *                                           0 Bad calculation that means that the current nonlinear iteration should be terminated
     */
    virtual int GFCEO_calcResid(double* const resid, const ResidEval_Type evalType);

    //!  Return a vector of delta y's for calculation of the numerical Jacobian
    /*!
     * (virtual from ZZCantera::ResidJacEval)
     *
     *   There is a default algorithm provided.
     *
     *        delta_y[i] = atol[i] + 1.0E-6 ysoln[i]
     *        delta_y[i] = atol[i] + MAX(1.0E-6 ysoln[i] * 0.01 * solnWeights[i])
     *
     * @param[in]            t                   Time                    (input)
     * @param[in]            ySoln               Solution vector (input, do not modify)
     * @param[in]            ySolnDot            Rate of change of solution vector. (input, do not modify)
     * @param[in]            deltaYSoln          Value of the delta to be used in calculating the numerical jacobian
     * @param[in]            solnWeights         Value of the solution weights that are used in determining convergence (default = 0)
     *
     * @return                                   Returns a flag to indicate that operation is successful.
     *                                           - 1  Means a successful operation
     *                                           - 0  Means an unsuccessful operation
     */
    virtual int calcDeltaSolnVariables(const double t, const double* const ySoln,
                                       const double* const ySolnDot, double* const deltaYSoln,
                                       const double* const solnWeights);

    //! Unpack the solution vector on return from the time stepper
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *  This function unpacks the solution vector into  phaseMoles_final_,  spMoles_final_, and spMf_final_[]
     *  If there is a bounds error, this routine will return an unsuccessful operation.
     *
     *  @param[in]           ySoln               Solution vector as returned from the time stepper.
     *                                           What the solution actually refers to depends on the individual Electrode objects.
     *                                           y[0] is always the current  deltaTsubcycleCalc_ value
     *
     *  @return                                   Returns a flag to indicate that operation is successful.
     *                                           - 1  Means a successful operation
     *                                           - 0  Means an unsuccessful operation
     */
    virtual int unpackNonlinSolnVector(const double* const ySoln);

    //! Check to see that the preceding step is a successful one
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *   We check to see if the preceding step is a successful one.
     *
     *  @return Returns a bool true if the step is acceptable, and false if it is unacceptable.
     */
    virtual bool  checkSubIntegrationStepAcceptable() const;

    //! Calculate the value of the solution derivatives, solnDot_final_[], at the end of the time step
    /*!
     *  We set the solnDot_final_[0] value to zero, as this is the calculated  deltaTsubcycleCalc_ value.
     *  If we don't have a good idea of the initial solnDot at the start of the subcycle, we copy solnDot_final_[] into solnDot_init_[]
     *  Right now, we only signal a bad initial solnDot for the entrance to this routine.
     *
     *  This routine must agree with calc_ydot() within the Zuzax' nonlinear solver, NonlinearSolver, for the methods to be consistent.
     *  However we zero the derivative of deltaT, because it doesn't make sense to carry that value around.
     */
    virtual void calc_solnDot_final();

     //! Calculate the current value of the solution derivatives ydotNLS_[] that are consistent with yvalNLS_[]
     /*!
      *  (virtual from Electrode_Integrator)
      *
      *  We calculate the consistent values of ydotNLS that are consistent with the current value of yvalNLS_[] during the nonlinear iterations
      *  It returns the derivative of ydot wrt y, which for the backwards Euler method is equal to 1/deltaTsubcycleCalc;
      *  This value is needed to compute analytical values of the Jacobian for time stepping problems.
      *
      *  This routine must agree with calc_ydot() within the Zuzax' nonlinear solver, NonlinearSolver
      *  ydotNLS_[] is updated within the nonlinear solver as yvalNLS_ is updated. This routine serves to initialize the values on
      *  going into the nonlinear solver.
      *  
      *  @return                                 Returns the derivative of ydot wrt y.
      */
     virtual double calc_ydotNLS_final();

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
    virtual double predictorCorrectorWeightedSolnNorm(const std::vector<double>& yvalNLS);

    //! On a completed local step, accumulate local errors into the global error vectors for the global time step
    /*!
     *
     */
    void accumulateLocalErrorToGlobalErrorOnCompletedStep();

    //! When the norm is high this routine decides what to do with the step
    /*!
     *  @param[in]    pnorm                            Norm of the step
     *  @param[in]    num_newt_its                     Number number of newton iterations
     *  @param[in]    iterSubCycle                     Number of iterations
     *
     *  @return   Returns whether the step is accepted or not
     */
    virtual bool decide_normHighLogic(double pnorm, int num_newt_its, int iterSubCycle);

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

    //! Create the L0 relative norm for the vector against a vector or predicted values
    /*!
     * This creates a relative norm that is scaled by the current value of rtol, and that can always be compared to 1
     * for satisfaction of the error criteria.
     *
     *                L0_i    =     || Actual_i - Pred_i || / ( rtol * MAX(||Actual_i||, || Pred_i ||)
     *
     *           or if ||Abstol_i|| > rtol * MAX(||Actual_i||, || Pred_i || 
     * 
     *                L0_i    =     || Actual_i - Pred_i || / ( MAX(||Abstol_i|| )
     *
     *  Also stores:
     *
     *                errorLocalNLS[i] = L0[i] * rtol
     *
     *  so that errorLocalNLS has the magnitude of the normalized error associated with it, mitigated by Abstol[i]
     *
     *  @param[in]           v1                  First vector to compare against. Length >= num
     *  @param[in]           v2                  Second vector to compare against. Length >= num
     *  @param[in]           num                 Length of the comparision. (usually equal to v1.size()
     *  @param[in]           atolVec             Absolute tolerance vector. Length >= num
     *  @param[in]           rtol                Relative tolerance
     * 
     *   @return                                 returns the maximum value of Lo_i for all i.
     */
    double l0norm_PC_NLS(const std::vector<double>& v1, const std::vector<double>& v2, size_t num,
                         const std::vector<double>& atolVec, const double rtol);

    //------------------------------------------------------------------------------------------------------------------
    // -------------------------------  SetState Functions -------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------

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

    // ----------------------------------------------------------------------------------------------
    // ----------------------------- GET CONDITIONS OUT ---------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------
    //
    //       (unless specified this is always at the final conditions and time
    //
    // ---------------------------------------------------------------------------------------------
    // ----------------------------- GET INSTANTANEOUS SOURCE TERMS ---------------------------------------------------------
    // ---------------------------------------------------------------------------------------------


    // ---------------------------------------------------------------------------------------------
    // ---------------------------- GET INTEGRATED SOURCE TERMS -------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------

    // ---------------------------------------------------------------------------------------------
    // --------------------------- GET MOLE NUMBERS ------------------------------------------------
    // ---------------------------------------------------------------------------------------------

    // ---------------------------------------------------------------------------------------------
    // --------------------------- GET THERMO AND VOLTAGE  -----------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------

    // ---------------------------------------------------------------------------------------------
    // ------------------------- GET VOLUMES -------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------

    // ---------------------------------------------------------------------------------------------
    // ---------------------- GET SURFACE AREAS ----------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------



    // ---------------------------------------------------------------------------------------------
    // --------------------------- INTERNAL UPDATE FUNCTIONS ---------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------

    // ---------------------------------------------------------------------------------------------
    // -------------------------------  SET STATE FUNCTIONS ----------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------

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
     *  (non-virtual function)  -> function should onionize in-first.
     *
     *  Set the final state from the init state. This is commonly called during a failed time step
     *
     */
    virtual void setFinalStateFromInit_Oin();

    //! Set the internal final intermediate and from the internal init state
    /*!
     *  (virtual function from Electrode)
     *
     *  Set the final state from the init state. This is commonly called during a failed time step
     *
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

    //! Set the internal final global state from the internal final intermediate state
    /*!
     *  (virtual function from Electrode)
     *
     *  Set the final_final state from the final state. This is commonly called at the end of successful base integration
     */
    virtual void setFinalFinalStateFromFinal();

    // ---------------------------------------------------------------------------------------------
    // ------------------------------ PRINT ROUTINES -----------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------

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

    //! Print a header for the residual printout
    /*!
     *  (virtual from Eelctrode_Integrator)
     */
    virtual void printResid_TimeHeader();

    //! Check for problems with the residual
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *  Checks here will cause the current nonlinear solve to fail
     *
     *    @return Return a negative value if there is a fatal problem.
     */
    virtual int residEval_BaseChecks();

    //!  Print details about the satisfaction of the residual
    /*!
     *  (virtual from Electrode_Integrator)
     */
    virtual void printResid_ResidSatisfaction();

    // Sets the base timeHistory record to the current timeHistory record
    /*!
     *  @return                                  Returns the number of time steps in the current timeHistory record
     */
    int setTimeHistoryBaseFromCurrent();

    //! Sets the timeHistory_base_ member to the argument object
    /*!
     *  During calculation of jacobians, the base time step history is followed.
     *
     *  @param[in]           timeHistory         SubIntegrationHistory object to set the base to
     *
     *  @return                                  Returns the number of steps in the timeHistory record
     */
    int setTimeHistoryBase(const SubIntegrationHistory& timeHistory);

    //! Returns an object containing the subcycle time integration history
    /*!
     *  Return the history record of the current substep time integration.
     *
     *  @param[in]           returnBase          If true this returns the timeHistory_base_ member.
     *                                           If false this returns the current timeHistory_current_ member.
     *
     *  @return                                  Returns a changeable SubIntegrationHistory object
     */
    SubIntegrationHistory& timeHistory(bool returnBase = false);

    //! Set the Maximum number of subcycles of integration before throwing an error 
    /*!
     *   This is the maximum number of subcycles no matter what the type of subsycles
     *
     *  @param[in]           maxN                Number of time step subcyles. Default value is 50000
     */
    void setMaxNumberSubCycles(int maxN);

    // ----------------------------------------------------------------------------------------------------------------
    // -------------------------------- DATA  -------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------

protected:

    //! Value of the computed time constant
    /*!
     *  This is a component of the solution vector using a trivial dirichlet equation. For special sub steps, it is solved for
     *  quantity indicating the deltaT needed for a phase to disappear from the solution domain.
     */
    double deltaTsubcycleCalc_;

    //! Absolute tolerance for nonlinear residual associated with time integration
    /*!
     *  length: neq_  
     *  Units:  Units of the individual components of the residual equations.
     *  Value depends on the units of the residual equations. This is determined from the Jacobian of the problem by default
     *  multiplied by acceptable delta_yvals.
     */
    std::vector<double> atolResidNLS_;

    //! Relative tolerance for nonlinear residual
    /*!
     *  Value is dimensionless
     */
    double rtolResidNLS_;

    //! Vector of absolute tolerances for the nonlinear solution values for the problem
    /*!
     *  length: ResidJacEval::neq_  
     *  Units:  Value of the solution components
     *  These have units particular to the value of the solution compoent
     */
    std::vector<double> atolNLS_;

    //! Relative tolerance for the solution variables for the nonlinear solution of the time stepping problem
    /*!
     *  Value is dimensionless
     */
    double rtolNLS_;

    //! Vector of low bounds on the solution variables for the nonlinear solver
    /*!
     *  Length: ResidJacEval::neq_  
     *  Units:  Value of the solution components
     *  values at all global and local times
     */
    std::vector<double> ylowNLS_;

    //! Vector of high bounds on the solution variables for the nonlinear solver
    /*!
     *  Length: ResidJacEval::neq_  
     *  Units:  Value of the solution components
     *  values at all global and local times
     */
    std::vector<double> yhighNLS_;

    //! Vector of solution values for the nonlinear solver for the current final local time step
    /*!
     *  Length: ResidJacEval::neq_  
     *  Units:  Value of the solution components
     *  values at t_n
     */
    std::vector<double> yvalNLS_;

    //! Vector of initial solution values for the nonlinear solver for the current  local time step
    /*!
     *  Length: ResidJacEval::neq_  
     *  Units:  Value of the solution components
     *  values at t_n-1
     */
    std::vector<double> yvalNLS_init_;

    //! Vector of initial solution values for the nonlinear solver for the current global time step
    /*!
     *  Length: ResidJacEval::neq_  
     *  Units:  Value of the solution components
     *  values at tGlob_m-1
     */
    std::vector<double> yvalNLS_init_init_;

    //! Vector of final solution values for the nonlinear solver for the current global time step
    /*!
     *  Length: ResidJacEval::neq_  
     *  Units:  Value of the solution components
     *  values at tGlob_m
     */
    std::vector<double> yvalNLS_final_final_;

  public:

    //! Vector of solution dot values for the nonlinear solver for the current final local time step
    /*!
     *  This is the vector that is kept up to date as the nonlinear solver converges on the solution by changing the
     *  value of yvalNLS_[]. ydotNLS_[] is kept consistent with yvalNLS_[] by calling calc_ydot() within the NonlinearSolver
     *  routine. After a successful converged step, ydotNLS_[] is transfered to solnDot_final_[]
     *
     *  Length: ResidJacEval::neq_  
     *  Units   Value of the solution components divided by deltaTsubgridcalc
     *  values at t_n
     */
    std::vector<double> ydotNLS_;

    //! Vector of normalized error values at the current local time step
    //! multiplied by rtolNLS_. Thus a 10^-3 value would indicate the error is 1 part in 1000.
    /*!
     *  Length: ResidJacEval::neq_  
     *  values at t_n
     */
    std::vector<double> errorLocalNLS_;

    //! Vector of normalized error values at the current global time step
    //! multiplied by rtolNLS_. Thus a 10^-3 value would indicate the error is 1 part in 1000.
    /*!
     *  This is the accumulated error for an unknown over multiple sub time steps for the global m time step
     *  Length: ResidJacEval::neq_  
     *  values at tGlob_m
     */
    std::vector<double> errorGlobalNLS_;

    //! Vector of normalized values used to scale the nonlinear solver damping strategy.
    /*!
     *  Step sizes that are smaller than this value are not controlled by numerical damping.
     *  Usually,  this is set to (1000 atol[i]). 
     *
     *  Length: ResidJacEval::neq_ 
     *  Units:  These are the solution unknown units for the time stepping problem
     *  values at all time steps.
     */
    std::vector<double> deltaBoundsMagnitudesNLS_;

    //! Boolean vector indicating a phase just died on this subgrid time integration step
    /*!
     *  Length:   m_NumTotPhases
     *  Units:    boolean 
     *  Values are reset to 0 at the start of every local subgrid time step
     */
    std::vector<int> phaseJustDied_;

    //! Boolean vector indicating a phase just was born on this subgrid time integration step
    /*!
     *  Length = m_NumTotPhases
     *  Units = boolean 
     *  Values are reset to 0 at the start of every local subgrid time step
     */
    std::vector<int> phaseJustBorn_;

    //! Vector of predicted solutions for the t_n given values at t_n-1 and previous
    /*!
     *  We use the equation system to come up with predictions. This approximation doesn't use solnDot to produce a guess.
     *  We augment this vector with the onBoundaryRegion value prediction at the end of the vector
     *
     *   soln_predict_[neq_] = onRegionBoundary_final_predicted;
     *
     *  Length: ResidJacEval::neq_ + 2
     *  Units:  These are the solution unknown units for the time stepping problem
     *  Values: These are predictions at t_n given previous conditions
     */
    std::vector<double> soln_predict_;

    //! Boolean indicating that solnDot_init_ is filled with an acceptable approximation of the solution derivatives.
    /*!
     *  This value starts out as false. Then, it turns positive after the first step. 
     *  When the value is false, we put the derivative of the current step (i.e., failed step) into solnDot_init_
     */
    bool haveGood_solnDot_init_;

    //! Vector of initial solution dot values for the nonlinear solver for the current local time step
    /*!
     *  Length: ResidJacEval::neq_  
     *  Units:  Solution components values divided by delta_t 
     *  values at t_n-1
     */
    std::vector<double> solnDot_init_;

    //! Vector of final solution dot values for the nonlinear solver for the current local time step
    /*!
     *  This is used to store the solnDot values at the end of a converged solution. During the solution
     *  phase, ydotNLS_[] is used.
     *
     *  Length: ResidJacEval::neq_  
     *  Units:  Solution components values divided by delta_t 
     *  values at t_n
     */
    std::vector<double> solnDot_final_;

    //! Vector of final solution dot values for the nonlinear solver for the current global time step
    /*!
     *  Length: ResidJacEval::neq_  
     *  Units:  Solution components values divided by delta_t
     *  values at tGlob_m
     */
    std::vector<double> solnDot_final_final_;

    //! Vector of initial solution dot values for the nonlinear solver for the current global time step
    /*!
     *  Length: ResidJacEval::neq_  
     *  Units:  Solution components values divided by delta_t 
     *  values at tGlob_m
     */
    std::vector<double> solnDot_init_init_;

    //! Vector of the predicted solution values at t_n using soln_init_[] and solnDot_init_[]
    /*!
     *  We calculate two predicted solution values, one from the equations and one from the solved solution
     *  and its derivative at the previous time step. This is the later guess. The best predictor is used to
     *  control the time step. 
     *
     *  Length: ResidJacEval::neq_  
     *  Units:  Solution components values divided by delta_t 
     *  values at t_n
     */
    std::vector<double> soln_predict_fromDot_;

    //! Current boolean flag for whether the solnDot predictor is better than the explicit predictSoln predictor.
    /*!
     *  We use this to figure out whether to use the predictSoln() or predictSolnDot() initial guess.
     *  At the start, we set this to false, because solnDot is zero. However, whenever we can, we use the last
     *  result from predictorCorrectorWeightedSolnNorm() to set this flag.
     */
    bool predictDotBetter_;

    //! Pointer to the nonlinear solver
    /*!
     *  We use Zuzax' nonlinear solver to relax the equations
     */
#ifdef useZuzaxNamespace
    Zuzax::NonlinearSolver* pSolve_;
#else
    Cantera::NonlinearSolver* pSolve_;
#endif

    //! Jacobian matrix
    SquareMatrix* jacPtr_;

    //! Number of degrees of freedom in the integrated source terms that constitute the output
    size_t numIntegratedSrc_;

    //! Predicted integrated source term vector for the local step
    std::vector<double> IntegratedSrc_Predicted;

    //! final integrated source term vector for the local step
    std::vector<double> IntegratedSrc_final_;

    //! Evaluation of the current error in the integrated source term for the current local step
    std::vector<double> IntegratedSrc_Errors_local_;

    //! Evaluation of the current error in the integrated source term for the current global step
    std::vector<double> IntegratedSrc_Errors_globalStep_;

    //! Relative tolerance for the integrated global src term vectors
    double rtol_IntegratedSrc_global_;

    //! Absolute tolerance for the integrated global src term vectors
    std::vector<double> atol_IntegratedSrc_global_;

    //! Maximum number of subcycles of integration before yielding an error message
    /*!
     *   This is the maximum number of subcycles no matter what the type of subsycles
     *
     *   Default value is 50000
     */
    int maxNumberSubCycles_;

    //! Maximum number of subGlobal time step iterations
    int maxNumberSubGlobalTimeSteps_;

    //! Norm of the error in the source terms for the current subgrid time step
    /*!
     *  This is the L2 norm of the estimated errors in the source term vector produced by this object
     */
    double IntegratedSrc_normError_local_;

    //! Norm of the sum of the error in the source terms for the current global time step
    /*!
     *  We add the errors for the local grid time steps
     */
    double IntegratedSrc_normError_global_;

    //! This is the base time step history for use with numerical jacobian calculations
    /*!
     *  A numerical jacobian calculation consists of a base calculation and a delta calculation, one where the external variable
     *  value is delta'd and the calculation is repeted, using the same time step history. This variable stores the base time
     *  step history so that it can be repeted again.
     */
    SubIntegrationHistory timeHistory_base_;

    //! This variable stores the current time step history for the current subgrid time step calculation
    /*!
     *  Stores the time step history of the current calculation
     */
    SubIntegrationHistory timeHistory_current_;

public:
    //! relative minimum time step ratio
    double relativeLocalToGlobalTimeStepMinimum_;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

