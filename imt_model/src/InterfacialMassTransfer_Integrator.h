/*
 * $Id: InterfacialMassTransfer_Integrator.h 192 2012-06-04 16:03:01Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _INTERFACIALMASSTRANSFER_INTEGRATOR_H
#define _INTERFACIALMASSTRANSFER_INTEGRATOR_H

#include "InterfacialMassTransfer.h"
#include "zuzax/integrators.h"
#include "zuzax/numerics/ResidJacEval.h"

namespace Zuzax
{

  class NonlinearSolver_JAC;
  class SquareMatrix;



  //! This class is a derived class used to model phase - change electrodes
  /*!
   * Complete problem statement
   *
   */
  class InterfacialMassTransfer_Integrator : public InterfacialMassTransfer , public Zuzax::ResidJacEval {
  public:
  
    // ---------------------------------------------------------------------------------------------
    // ----------------------- BASIC SETUP ROUTINES  -----------------------------------------------
    // ---------------------------------------------------------------------------------------------

    //! Constructor
    InterfacialMassTransfer_Integrator();

    //! Destructor
    virtual  ~InterfacialMassTransfer_Integrator();

    //! Copy Constructor
    /*!1
     * @param right Object to be copied
     */
    InterfacialMassTransfer_Integrator(const InterfacialMassTransfer_Integrator &right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    InterfacialMassTransfer_Integrator & operator=(const InterfacialMassTransfer_Integrator &right);
				
    //! Create the electrode model
    int model_create(IMT_KEY_INPUT *ei);


    //!  Set the initial conditions from the input file.
    /*!   
     *   (virtual from InterfacialMassTransfer)
     *   (This is a serial virtual function or an overload function)
     *
     *    This is one of the most important routines. It sets up the initial conditions of the interface
     *    from the input file. The interface object itself has been set up from a call to model_create().
     *    After the call to this routine, the interface should be internally ready to be integrated and reacted. 
     *    It takes its input from an IMT_KEY_INPUT object which specifies the setup of the interface
     *    object and the initial state of that object.
     *    
     *    The routine works like an onion initialization. The parent object is initialized before the 
     *    child. This means the child object first calls the parent, before it does its own initializations.
     * 
     * @param ei    IMT_KEY_INPUT pointer object
     *  
     *  @return  Returns zero if successful, and -1 if not successful.
     */
    virtual int setInitialConditions(IMT_KEY_INPUT *ei);

    // ---------------------------------------------------------------------------------------------
    // ----------------------- SPECIFY AND OBTAIN PROBLEM PARAMETERS -------------------------------
    // ---------------------------------------------------------------------------------------------

    // ------------------------------ OBTAIN STATIC PROBLEM INFORMATION ----------------------------


    virtual int create_solvers();

    //! Setup the vectors for handling the error control on source terms
    /*!
     *  This routine formally identifies what source terms will be error-controlled.
     *
     *  The default is to assume that m_NumTotSpecies are created and the vector is called
     *  spMoleIntegratedSourceTerm_
     * 
     *  @return Returns the number of degrees of freedom actually created.
     */
    virtual int setupIntegratedSourceTermErrorControl();

    // ------------------------------ SPECIFY BASIC THERMO CONDITIONS  ------------------------------------
 

    // ------------------------------ SPECIFY PROBLEM PARAMETERS ------------------------------------

 
 
    // ---------------------------------------------------------------------------------------------
    // ----------------------------- CARRY OUT INTEGRATIONS -----------------------------------------
    // ---------------------------------------------------------------------------------------------

    //! The internal state of the electrode must be kept for the initial and 
    //! final times of an integration step.
    /*
     *  This function advances the initial state to the final state that was calculated
     *  in the last integration step. 
     *
     * @param Tinitial   This is the New initial time. This time is compared against the "old"
     *                   final time, to see if there is any problem.
     */
    virtual void  resetStartingCondition(double Tinitial);

    //!  Calculate the change in the state of the system when integrating from Tinitial to Tfinal
    /*!
     *  All information is kept internal within this routine. This may be done continuously
     *  and the solution is not updated. 
     *
     *  Note the tolerance parameters refere to the nonlinear solves within the calculation
     *  They do not refer to time step parameters.
     *
     *  @param deltaT              DeltaT for the integration step.
     *  @param GlobalRtolSrcTerm   Relative tolerance for nonlinear solves within the calculation
     *                             Defaults to 1.0E-3
     *  @param GlobalAtolSrcTerm   Absolute tolerance for nonlinear solves within the calculation
     *                             Defaults to 1.0E-12
     *
     *  @param fieldInterpolationType
     *
     *  @param subIntegrationType
     *  @param sih
     *
     *  @return Returns the number of subcycle steps it took to complete the full step.
     */
    virtual int integrate(double deltaT, double GlobalRtolSrcTerm = 1.0E-3, double GlobalAtolSrcTerm = 1.0E-12,
                          int fieldInterpolationType = 0, int subIntegrationType = 0, SubIntegrationHistory * sih = 0);

    //! Calculate the largest mole fraction in each of the phases
    /*!
     *  We use this to determine the equation system
     */
    virtual void determineBigMoleFractions();

    //! Check the consistency of the initial state
    /*!
     *  (virtual from InterfacialMassTransfer_Integrator) 
     *
     *   This is a checker routine. Therefore, it can be taken out in nondebug mode
     */
    virtual void check_init_consistency();

    //! Set the base tolerances for the nonlinear solver within the integrator for the current global step
    /*!
     *   The tolerances are based on controlling the integrated electron source term 
     *   for the electrode over the integration interval.  The integrated source term
     *   has units of kmol.
     *
     *   Because the electron is only one molar quantity within a bunch of molar quantities,
     *   this requirement will entail that we control the source terms of all species within the
     *   electrode to the tolerance requirements of the electron source term. 
     *
     *   @param rtolResid  Relative tolerance allowed for the source term over the interval.
     *                     This is a unitless quantity
     *   @param atolResid  Absolute tolerance cutoff for the source term over the interval.
     *                     Below this value we do not care about the results.
     *                     atol has units of kmol. 
     */
    virtual void setNLSGlobalSrcTermTolerances(double rtolResid, double atolResid);

    //! Returns the number of integrated source terms whose errors are controlled by this integrator
    int numIntegratedSrcTerms() const;

    //! Set the Residual absolute error tolerances
    /*!
     *  (virtual from InterfacialMassTransfer_Integrator)
     *
     *   Set the absolute error tolerances fror the nonlinear solvers. This is called at the top
     *   of the integrator() routine.
     *
     *   Calculates residAtolNLS_[]
     */
    virtual void setResidAtolNLS(double GlobalRtolSrcTerm);

    //! Predict the solution
    /*!
     *  (virtual from InterfacialMassTransfer_Integrator)
     *
     * Ok at this point we have a time step 
     * and initial conditions consisting of phaseMoles_init_ and spMF_init_.
     * We now calculate predicted solution components from these conditions.
     * Additionally, we evaluate whether any multispecies phases are going to pop into existence,
     * adding in potential seed values, or die, creating a different equation set and time step equation.
     * We set the phaseExistence flags in the kinetics solver to reflect phase pops.
     *
     * @return   Returns the success of the operation
     *                 1  A predicted solution is achieved 
     *                 2  A predicted solution with a multispecies phase pop is acheived
     *                 0  A predicted solution is not achieved, but go ahead anyway
     *                -1  The predictor suggests that the time step be reduced and a retry occur.
     */
    virtual int predictSoln();

    //! Extract information from cantera
    /*!
     *  (virtual fucntion from InterfacialMassTransfer_Integrator)
     */
    virtual void extractInfo();

    //! Collect mole change information
    /*!
     *  (virtual fucntion from InterfacialMassTransfer_Integrator)
     */
    virtual void updateSpeciesMoleChangeFinal();

    //! calculate the residual
    /*!
     *   
     *  (virtual fucntion from InterfacialMassTransfer_Integrator)
     *
     *  @return Returns 1 if everything is ok.
     *          Returns 0 if the current conditions can not be calculated.
     */
    virtual int calcResid(doublevalue * const resid, const ResidEval_Type evalType);

    //!  Gather the predicted solution values and the predicted integrated source terms
    /*!
     *  (virtual from InterfacialMassTransfer_Integrator)
     *
     *  Both the predicted solution values and the predicted integrated source terms are used
     *  in the time step control
     */
    virtual void gatherIntegratedSrcPrediction();

      
    //!  Calculate the integrated source terms and do other items now that we have a completed time step
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
     */
    virtual void accumulateSrcTermsOnCompletedStep();

 

    //! Pack the nonlinear solver problem
    /*!
     *  formulate the nonlinear solver problem to be solved.
     *     Fields to be filled in
     *             yvalNLS_
     *             ylowNLS_
     *             yhighNLS_
     *             atolNLS_
     *             deltaBoundsMagnitudesNLS_     
     */
    virtual void initialPackSolver_nonlinFunction();

    //! Check the nonlinear residual equations for completeness and the ability to be solved
    /*!
     *
     */
    virtual void check_nonlinResidConditions();
 
 
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
    virtual int evalResidNJ(const doublevalue t, const doublevalue delta_t,
			    const doublevalue * const y,
			    const doublevalue * const ydot,
			    doublevalue * const resid,
			    const ResidEval_Type evalType = ResidEval_Type::Base_ResidEval,
			    const int id_x = -1, 
			    const doublevalue delta_x = 0.0);

    //! Fill in the initial conditions
    /*! 
     * (virtual from NonlinearSolver)
     *
     * Values for both the solution and the value of ydot may be provided.
     *
     * @param t0            Time                    (input) 
     * @param y             Solution vector (output)
     * @param ydot          Rate of change of solution vector. (output)
     */
    virtual int getInitialConditions(const doublevalue t0, doublevalue * const y, doublevalue * const ydot);
      
      
    //! Return the number of equations in the equation system
    /*!
    * (virtual from NonlinearSolver)
     *
     */
    virtual int nEquations() const;

    //!  Return a vector of delta y's for calculation of the numerical Jacobian 
    /*!  
     * (virtual from Zuzax::ResidJacEval)
     *
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
    virtual int calcDeltaSolnVariables(const doublevalue t, const doublevalue * const ySoln,
				       const doublevalue * const ySolnDot, doublevalue * const deltaYSoln,
				       const doublevalue *const solnWeights);

  
    //! Unpack the soln vector
    /*!
     *  (virtual from InterfacialMassTransfer_Integrator)
     *
     *  This function unpacks the solution vector into  phaseMoles_final_,  spMoles_final_, and spMf_final_[]
     */
    virtual void unpackNonlinSolnVector(const double * const y);

    //! Check to see that the preceding step is a successful one
    /*!
     *  (virtual from InterfacialMassTransfer_Integrator)
     *
     *   We check to see if the preceding step is a successful one.
     *
     *  @return Returns a bool true if the step is acceptable, and false if it is unacceptable.
     */
    virtual bool checkSubIntegrationStepAcceptable() const;

    //!  Calculate the norm of the difference between the predicted answer and the final converged answer 
    //!  for the current time step
    /*!
     *  (virtual from InterfacialMassTransfer_Integrator)
     *
     *   The norm calculated by this routine is used to determine whether the time step is accurate enough.
     *
     *  @return    Returns the norm of the difference. Normally this is the L2 norm of the difference
     */  
    virtual double predictorCorrectorWeightedSolnNorm(const std::vector<double> & yvalNLS);

    //! Calculate the vector of predicted errors in the source terms that this integrator is responsible for
    /*!
     *  (virtual from InterfacialMassTransfer_Integrator)
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
     *  (virtual from InterfacialMassTransfer_Integrator)
     *
     *   This routine make use of the source term error vector along with rtols and atols for the 
     *   individual source terms to calculated a normalized error measure. This is the single number
     *   that the integration routine will try to control as it calculates a time stepping strategy.
     *    
     *   @return  Returns a single nondimensional number representing the normalized error
     *            for the calculation of the source term
     */
    virtual double predictorCorrectorGlobalSrcTermErrorNorm();

  


    //! Possibly change the solution due to phase births and deaths. 
    /*!
     *   (virtual from InterfacialMassTransfer_Integrator)
     *
     *  @return  Returns true if the solution step is still good. It returns false if there is a problem. 
     */
    virtual bool changeSolnForBirthDeaths();

    //! Possibly change the solution due to phase births and deaths after phase has been accepted.
    /*!
     *   (virtual from InterfacialMassTransfer_Integrator)
     *
     *  This routine is carried out after the step is deemed a success. Massaging of the solution
     *  must be carried out within strict tolerances.
     */
    virtual void manageBirthDeathSuccessfulStep();

    //! Error check on the routine step
    /*!
     *    (virtual from InterfacialMassTransfer_Integrator)
     *
     *   Error checks go here. All errors are fatal exits.
     */
    virtual void check_final_state();

    //! Print a header for the residual printout
    /*!
     *  (virtual from Eelctrode_Integrator)
     */
    virtual void printResid_TimeHeader();

    //! Check for problems with the residual
    /*!
     *  (virtual from InterfacialMassTransfer_Integrator)
     *
     *  Checks here will cause the current nonlinear solve to fail
     *
     *    @return Return a negative value if there is a fatal problem. 
     */
    virtual int residEval_BaseChecks();

    // ----------------------------------------------------------------------------------------------
    // ----------------------------- GET CONDITIONS OUT --------------------------------------------
    // ----------------------------------------------------------------------------------------------
    //
    //       (unless specified this is always at the final conditions and time
    //

    // ----------------------------- GET INSTANTANEOUS SOURCE TERMS --------------------------------


    //! Get the net production rates of all species in the object
    //! at the current final conditions from all surface kinetics objects
    /*!
     * @param net   Species net production rates [kmol/s]. Return the species
     */
    void getNetProductionRates(doublevalue* const net) const;

    //! Get the net production rates of all species in the object
    //! at the current final conditions from one surface kinetics object
    /*!
     * @param isk   Surface index to get the net production rates from
     * @param net   Species net production rates [kmol/s]. Return the species
     */
    void getNetProductionRatesRSD(const size_t isk, doublevalue* const net) const;

    //!  Returns the current and the net production rates of the phases in kg/m2/s from a single surface
    //!  at the current final conditions and t_final
    /*!
     *  Returns the net production rates of all phases from reactions on a single surface
     *  
     *  @param isk Surface ID to get the fluxes from.      
     *  @param phaseMassFlux  Returns the mass fluxes of the phases
     */
    void getPhaseMassFlux(doublevalue* const phaseMassFlux) const;

    //!  Returns the current and the net production rates of the phases in kmol/m2/s from a single surface
    //!  at the current final conditions and t_final
    /*!
     *  Returns the net production rates of all phases from reactions on a single surface
     *  
     *  @param isk Surface ID to get the fluxes from.      
     *  @param phaseMassFlux  Returns the mass fluxes of the phases
     */
    void getPhaseMoleFlux(const size_t isk, doublevalue* const phaseMoleFlux) const;

    //! Returns the computed Stefan velocity for phase A in the object
    //! at the current final conditions. 
    /*!
     *  Multiply by the surface area to get the stefan volume production rate at the
     *  current final conditions.
     * 
     *    @return returns stefan velociy created in phase A (m s-1)
     */
    virtual double StefanVelocityPhaseA() const;

    //! Returns the integrated volume creation rate for phase B in the object
    //! over the current global time step
    /*!
     *  This can be turned into a velocity by dividing by the surface area within the
     *  object.
     *    @return returns volume created in phase B (m3 s-1)
     */
    virtual double StefanVelocityPhaseB() const;


    // ---------------------------- GET INTEGRATED SOURCE TERMS -------------------------------------

    // --------------------------- GET MOLE NUMBERS ------------------------------------------------

    // --------------------------- GET THERMO AND VOLTAGE  ------------------------------------------------


    // ------------------------- GET VOLUMES -----------------------------------------------------------

  
    // ---------------------- GET SURFACE AREAS -------------------------------------------------------


    // --------------------------- INTERNAL UPDATE FUNCTIONS --------------------------------------

    //! Update the velocities within the object at the final conditions.
    /*!
     *  (virtual from InterfacialMassTransfer)
     * 
     *  In order to do this, we need the instantaenous source terms and the current value of the
     *  surface area over the local time step. Therefore, this routine must be called 
     *  after DspMoles_final_ has been calculated. 
     */
    virtual void updateVelocities();

    // ---------------------------------------------------------------------------------------------
    // -------------------------------  SetState Functions -------------------------------------------------------
    // ---------------------------------------------------------------------------------------------

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
     *  (virtual function from InterfacialMassTransfer) 
     *
     *  Set the final state from the init state. This is commonly called during a failed time step
     *
     */
    virtual void setFinalStateFromInit();

    // ---------------------------------------------------------------------------------------------
    // ------------------------------ PRINT ROUTINES -----------------------------------------------
    // ---------------------------------------------------------------------------------------------
    //! Print table representing prediction vs. corrector information
    /*!
     *  @param yval           Vector of corrector values
     *  @param pnormSrc       Norm of the predictor-corrector comparison for the source vector.
     *  @param pnormSoln      Norm of the predictor-corrector comparison for the solution vector.
     */
    virtual void predictorCorrectorPrint(const std::vector<double> &yval, 
					 double pnormSrc, double pnormSoln);

    //!  Print details about the satisfaction of the residual
    /*!
     *  (virtual from InterfacialMassTransfer_Integrator)
     */
    virtual void printResid_ResidSatisfaction();


    //! Print conditions of the electrode for the current integration step to stdout
    /*!
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values 
     */
    virtual void printInterfacialMassTransfer(int pSrc = 1, bool subTimeStep = false);

    //! Print condition of a phase in the electrode
    /*!
     *  @param iPhase        Print the phase
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values 
     */
    virtual void printInterfacialMassTransferPhase(size_t iPhase, int pSrc = 1, bool subTimeStep = false);
 
    // ---------------------------------------------------------------------------------------------
    // --------------------------------- MEMBER DATA -----------------------------------------------
    // ---------------------------------------------------------------------------------------------
  protected:

    //! value of the comptued time constant
    doublevalue deltaTsubcycleCalc_;
 

    //! Absolute tolerance for nonlinear residual
    std::vector<doublevalue> residAtolNLS_;

    //! Relative tolerance for nonlinear residual
    doublevalue rtolResidNLS_;

    //! Vector of absolute tolerances for the nonlinear solution values for the problem
    std::vector<doublevalue> atolNLS_;

    //!  Vector of relative tolerances for the solution variables for the nonlinear solution of the
    //!  time stepping problem
    doublevalue rtolNLS_;

    //! Vector of low bounds on the solution variables for the nonlinear solver
    std::vector<doublevalue> ylowNLS_;

    //! Vector of high bounds on the solution variables for the nonlinear solver
    std::vector<doublevalue> yhighNLS_;

    std::vector<doublevalue> yvalNLS_;
    std::vector<doublevalue> ydotNLS_;

  

    std::vector<doublevalue> deltaBoundsMagnitudesNLS_;

    //! Boolean vector indicating a phase just died on this subgrid integration step
    /*!
     *  Length = m_NumTotPhases
     */
    std::vector<size_t> phaseJustDied_;

    std::vector<size_t> phaseJustBorn_;

    //! Vector of predicted solutions
    /*!
     *  we augment this vector with the onBoundaryRegion value prediction.
     */
    std::vector<double> soln_predict_;

    //! Pointer to the nonlinear solver
    NonlinearSolver_JAC * pSolve_;

    //! Jacobian matrix
    SquareMatrix * jacPtr_;

    //! Number of degrees of freedom in the integrated source terms that are controlled
    //! in terms of accuracy
    int numIntegratedSrc_;

    //! Predicted integrated source term vector for the local step
    std::vector<double> IntegratedSrc_Predicted;

    //! final integrated source term vector for the local step
    std::vector<double> IntegratedSrc_final_;

    //! Evaluation of the current error in the integrated source term for the current local step
 
    std::vector<double> IntegratedSrc_Errors_local_;

    //! Evaluation of the current error in the integrated source term for the current global step
    std::vector<double> IntegratedSrc_Errors_globalStep_;

    //! Relative tolerance for the integrated global src term vectors
    doublevalue rtol_IntegratedSrc_global_;

    //! Absolute tolerance for the integrated global src term vectors
    std::vector<double> atol_IntegratedSrc_global_;

    
    doublevalue IntegratedSrc_normError_local_;
    doublevalue IntegratedSrc_normError_global_;

    char buf_[1024];

    //! Vector of rates of progress from all Reacting domains in the object
    /*!
     *  This is filled in 
     */
    std::vector< std::vector<double> > ROP_RSD_List_;

    //! flag to indicate we have zero reaction rate
    int goNowhere_;

  public:
    //! relative minimum time step ratio
    double relativeLocalToGlobalTimeStepMinimum_;

  };

}


#endif 
/*****************************************************************************/
