/*
 * $Id: InterfacialMassTransfer_Integrator.cpp 192 2012-06-04 16:03:01Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"



#include "cantera/base/mdp_allo.h"

#include "InterfacialMassTransfer_Integrator.h"
#include "cantera/integrators.h"
#include "cantera/numerics/NonlinearSolver.h"

using namespace Cantera;
using namespace std;
using namespace BEInput;
using namespace TKInput;
using namespace mdpUtil;

#ifndef SAFE_DELETE
#define SAFE_DELETE(x)  if (x) { delete x;  x = 0;}
#endif

#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif           

namespace Cantera {
  //======================================================================================================================
  SubIntegrationHistory::SubIntegrationHistory() :
    nTimeSteps_(0),
    TimeStepList_(0),
    GsolnNorm_(0.0),
    GwdotNorm_(0.0)
  {

  }
  //======================================================================================================================
  SubIntegrationHistory::SubIntegrationHistory(const SubIntegrationHistory &right) :
    nTimeSteps_(right.nTimeSteps_),
    TimeStepList_(TimeStepList_),
    GsolnNorm_(right.GsolnNorm_),
    GwdotNorm_(right.GwdotNorm_)
  {

  }
  //======================================================================================================================
  SubIntegrationHistory::~SubIntegrationHistory() 
  {

  }
  //======================================================================================================================
  SubIntegrationHistory &SubIntegrationHistory::operator=(const SubIntegrationHistory &right)
  {
    if (this == &right) {
      return *this;
    }
    nTimeSteps_ = right.nTimeSteps_;
     TimeStepList_ = right.TimeStepList_;
     GsolnNorm_ = right.GsolnNorm_;
     GwdotNorm_ = right.GwdotNorm_;
     return *this;
  }
  //======================================================================================================================
  /*
   *  ELECTRODE_INPUT: constructor
   *
   *  We initialize the arrays in the structure to the appropriate sizes.
   *  And, we initialize all of the elements of the arrays to defaults.  
   */
  InterfacialMassTransfer_Integrator::InterfacialMassTransfer_Integrator() :
    InterfacialMassTransfer(),
    ResidJacEval(),
    deltaTsubcycleCalc_(0),
    residAtolNLS_(0),
    rtolResidNLS_(1.0E-6),
    atolNLS_(1.0E-6),
    ylowNLS_(0),
    yhighNLS_(0),
    yvalNLS_(0),
    ydotNLS_(0),
    deltaBoundsMagnitudesNLS_(0),
    phaseJustDied_(0),
    phaseJustBorn_(0),
    pSolve_(0),
    jacPtr_(0),
    numIntegratedSrc_(0),
   
    IntegratedSrc_Predicted(0),
    IntegratedSrc_final_(0),
    IntegratedSrc_Errors_local_(0),
    IntegratedSrc_Errors_globalStep_(0),
    rtol_IntegratedSrc_global_(1.0E-4),
    atol_IntegratedSrc_global_(0),
    IntegratedSrc_normError_local_(0.0),
    IntegratedSrc_normError_global_(0.0),
    ROP_RSD_List_(0),
    goNowhere_(0),
    relativeLocalToGlobalTimeStepMinimum_(0.0001)
  {
    /*
     * Set the number of equations in ResidJacEval. This will be used in the nonlinear solver
     */
    neq_ = 1;
  }
  //======================================================================================================================
  // Copy Constructor
  /*
   * @param right Object to be copied
   */
  InterfacialMassTransfer_Integrator::InterfacialMassTransfer_Integrator(const InterfacialMassTransfer_Integrator &right) :
    InterfacialMassTransfer(), 
    ResidJacEval(),
    deltaTsubcycleCalc_(0),
    residAtolNLS_(0),
    rtolResidNLS_(1.0E-6),
    atolNLS_(1.0E-6),
    ylowNLS_(0),
    yhighNLS_(0),
    yvalNLS_(0),
    ydotNLS_(0),
    deltaBoundsMagnitudesNLS_(0),
    phaseJustDied_(0),
    phaseJustBorn_(0),
    pSolve_(0),
    jacPtr_(0),
    numIntegratedSrc_(0),
    IntegratedSrc_Predicted(0),
    IntegratedSrc_final_(0),
    IntegratedSrc_Errors_local_(0),
    IntegratedSrc_Errors_globalStep_(0),
    rtol_IntegratedSrc_global_(1.0E-4),
    atol_IntegratedSrc_global_(0),
    IntegratedSrc_normError_local_(0.0),
    IntegratedSrc_normError_global_(0.0),
    ROP_RSD_List_(0),
    goNowhere_(0),
    relativeLocalToGlobalTimeStepMinimum_(0.0001)
  { 
    /*
     *  set the number of equations
     */
    neq_ = 1;
    /*
     * Call the assignment operator.
     */
    *this = operator=(right);
  }
  //======================================================================================================================
  // Assignment operator
  /*
   *  @param right object to be copied
   */
  InterfacialMassTransfer_Integrator &
  InterfacialMassTransfer_Integrator::operator=(const InterfacialMassTransfer_Integrator &right)
  {
    /*
     * Check for self assignment.
     */
    if (this == &right) return *this;
 
    InterfacialMassTransfer::operator=(right);
    Cantera::ResidJacEval::operator=(right);

    neq_                                = right.neq_;
    deltaTsubcycleCalc_                 = right.deltaTsubcycleCalc_;

    residAtolNLS_                       = right.residAtolNLS_;
    rtolResidNLS_                       = right.rtolResidNLS_;
    atolNLS_                            = right.atolNLS_;
    rtolNLS_                            = right.rtolNLS_;
    ylowNLS_                            = right.ylowNLS_;
    yhighNLS_                           = right.yhighNLS_;
    yvalNLS_                            = right.yvalNLS_;
    ydotNLS_                            = right.ydotNLS_;
    deltaBoundsMagnitudesNLS_           = right.deltaBoundsMagnitudesNLS_;

    phaseJustDied_                      = right.phaseJustDied_;
    phaseJustBorn_                      = right.phaseJustBorn_;
   
    SAFE_DELETE(pSolve_);
    pSolve_                             = new NonlinearSolver(this);
    SAFE_DELETE(jacPtr_);
    jacPtr_                             = new SquareMatrix(*right.jacPtr_);

    numIntegratedSrc_                   = right.numIntegratedSrc_;
    IntegratedSrc_Errors_local_         = right.IntegratedSrc_Errors_local_;
    IntegratedSrc_Predicted             = right.IntegratedSrc_Predicted;
    IntegratedSrc_final_                = right.IntegratedSrc_final_;
    IntegratedSrc_Errors_globalStep_    = right.IntegratedSrc_Errors_globalStep_;
    rtol_IntegratedSrc_global_          = right.rtol_IntegratedSrc_global_;
    atol_IntegratedSrc_global_          = right.atol_IntegratedSrc_global_;
    IntegratedSrc_normError_local_      = right.IntegratedSrc_normError_local_;
    IntegratedSrc_normError_global_     = right.IntegratedSrc_normError_global_;
    ROP_RSD_List_                       = right.ROP_RSD_List_;
    goNowhere_                          = right.goNowhere_;
    relativeLocalToGlobalTimeStepMinimum_ = right.relativeLocalToGlobalTimeStepMinimum_;
    /*
     * Return the reference to the current object
     */
    return *this;
  }
  //======================================================================================================================
  /*
   *
   *  ELECTRODE_INPUT:destructor
   *
   * We need to manually free all of the arrays.
   */
  InterfacialMassTransfer_Integrator::~InterfacialMassTransfer_Integrator() 
  {
    SAFE_DELETE(jacPtr_);
    SAFE_DELETE(pSolve_);
  }
  //======================================================================================================================
  int 
  InterfacialMassTransfer_Integrator::model_create(IMT_KEY_INPUT *ei) {
    
    InterfacialMassTransfer::model_create(ei);

    
    phaseJustDied_.resize(m_NumTotPhases, 0);
    phaseJustBorn_.resize(m_NumTotPhases, 0);
 
    setupIntegratedSourceTermErrorControl();

    /*
     *  Create the ROP list that will be used later to gather results in extractInfo()
     */
    ROP_RSD_List_.clear();
    for (size_t ii = 0; ii < RSD_List_.size(); ii++) {
      int nr = RSD_List_[ii]->nReactions();
      std::vector<double> netROP(nr, 0.0);
      ROP_RSD_List_.push_back(netROP);
    }

    return 0;
  }
  //===================================================================================================================
  //  Set the initial conditions from the input file.
  /*   
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
  int  InterfacialMassTransfer_Integrator::setInitialConditions(IMT_KEY_INPUT *ei)
  {
    return InterfacialMassTransfer::setInitialConditions(ei);
  }
  //===================================================================================================================
  int 
  InterfacialMassTransfer_Integrator::create_solvers() {
    
    int neqNL = nEquations();
    int neq = neqNL;
    if (neqNL == 0) {
      neq = 1;
    }
      
    /*
     * Resize dimensions of other workspaces needed to solve the nonlinear solver.
     */
    atolNLS_.resize(neq, 1.0E-16);
    residAtolNLS_.resize(neq, 1.0E-16);

    yhighNLS_.resize(neq, 0.0);
    ylowNLS_.resize(neq, 0.0);

    yvalNLS_.resize(neq, 0.0);
    ydotNLS_.resize(neq, 0.0);
    deltaBoundsMagnitudesNLS_.resize(neq, 0.0);
    soln_predict_.resize(neq, 0.0);
    jacPtr_ = new SquareMatrix(neq, neq);

    pSolve_ = new NonlinearSolver(this);
    return neqNL;
  }
  //===================================================================================================================
  //! Setup the vectors for handling the error control on source terms
  /*!
   *  This routine formally identifies what source terms will be error-controlled.
   *
   *  The default is to assume that m_NumTotSpecies are created and the vector is called
   *  spMoleIntegratedSourceTerm_
   * 
   *  @return Returns the number of degrees of freedom actually created.
   */
  int InterfacialMassTransfer_Integrator::setupIntegratedSourceTermErrorControl()
  {
  
    numIntegratedSrc_ = m_NumTotSpecies;

    IntegratedSrc_Errors_local_.resize(numIntegratedSrc_, 0.0);
    IntegratedSrc_Predicted.resize(numIntegratedSrc_, 0.0);
    IntegratedSrc_final_.resize(numIntegratedSrc_, 0.0);
    IntegratedSrc_Errors_globalStep_.resize(numIntegratedSrc_, 0.0);

    rtol_IntegratedSrc_global_ = 1.0E-5;
    atol_IntegratedSrc_global_.resize(numIntegratedSrc_, 0.0);
    /*
     *   The default tolerances is to set the kmol tolerances to 1E-14 times the magnitude of the moles in the solid
     *   However, this may not be a good estimate. 
     *   The absolute tolerance is actually set by the numerical roundoff created by adding forward and reverse
     *   rates of progress together on reaction rates that are near equilibrium. Basically, this means taking 1.0E-14 
     *   times the production rate of these reactions as the absolute tolerance for the calculation. 
     */
    double sm = totalMoles();
    for (int i = 0; i <  numIntegratedSrc_; i++) {
      atol_IntegratedSrc_global_[i] = 1.0E-14 * sm;
    }

    return  numIntegratedSrc_;
  }
  //===================================================================================================================
  //    The internal state of the electrode must be kept for the initial and 
  //    final times of an integration step.
  /*
   *  This function advances the initial state to the final state that was calculated
   *  in the last integration step. 
   *
   * @param Tinitial   This is the New initial time. This time is compared against the "old"
   *                   final time, to see if there is any problem.
   */
  void  InterfacialMassTransfer_Integrator::resetStartingCondition(doublereal Tinitial) 
  {
    /*
     * If the initial time is input, then the code doesn't advance
     */
    double tbase = MAX(t_init_init_, 1.0E-50);
    if (fabs(Tinitial - t_init_init_) < (1.0E-9 * tbase)) {
      return;
    }
    InterfacialMassTransfer::resetStartingCondition(Tinitial);

  } 
  //====================================================================================================================
  //  Calculate the change in the state of the system when integrating from Tinitial to Tfinal
  /*
   *  All information is kept internal within this routine. This may be done continuously
   *  and the solution is not updated. 
   *
   *  Note the tolerance parameters refere to the nonlinear solves within the calculation
   *  They do not refer to time step parameters.
   *
   *  @param deltaT        DeltaT for the integration step.
   *  @param GlobalRtolResid     Relative tolerance allowed for the electron source term over the interval.
   *                       This is a unitless quantity 
   *                       Defaults to 1.0E-3
   *  @param GlobalAtolResid     Absolute tolerance cutoff for the electron source term over the interval.
   *                       Below this value we do not care about the results.
   *                       atolResid has units of kmol.
   *                       Defaults to 1.0E-12
   *
   *  @return Returns the number of subcycle steps it took to complete the full step.
   */
  int InterfacialMassTransfer_Integrator::
  integrate(double deltaT, double  GlobalRtolSrcTerm, double GlobalAtolSrcTerm,
	    int fieldInterpolationType, int subIntegrationType, SubIntegrationHistory * sih)
  {

    counterNumberIntegrations_++;
    counterNumberLastSubIntegrations_ = 0;
    //double atolVal = 1.0E-20;	
    
    int nsteps_est;
    double pnorm = 0.0;

    /*
     *  Guard against a zero deltaT
     */
    if (deltaT <= 0.0) {
      return 0;
    }

    /*
     *  When we call this routine successfully we have an integration for the current step pending
     *  Tempory data now exist for t_final_
     */
    pendingIntegratedStep_ = 1;
#ifdef DEBUG_HKM_FAILURES
    // checking out capability to recover from nonlin failures
 
#endif
    /*
     *  Create matrices, vectors and nonlinear solver for the integration procedure, if we haven't
     *  done so already
     */
    if (pSolve_ == 0) {
      create_solvers();
    }

    /*
     * We do a complete subcycle system here, even though the default is to use the same step as 
     *  the main calculation. There will be times when we break up the step, and there may be
     *  problems where subcycling is warranted.
     */
    setNLSGlobalSrcTermTolerances(GlobalRtolSrcTerm, GlobalAtolSrcTerm); 

    /*
     *  Set the absolute tolerance vector for the nonlinear solver
     */
    setResidAtolNLS(GlobalRtolSrcTerm);

    /*
     *  Set the counter for the number of subcycles
     */
    int iterSubCycle = 0; 
    /*
     *  Set the counter for the number of consequative failures of the time stepping method that happens for any reason
     */
    int conseqFailures = 0;
    // Set the flag for overall time stepping to the goal of t_final_final_
    bool notDone = true;

    /*
     *  Set the inital time steps and final_final time step based on the overall step size, deltaT
     */
    t_init_ = t_init_init_;
    t_final_ = t_init_;
    t_final_final_ = t_init_init_ + deltaT;

    /*
     * Choose the starting value of the subcycle time step
     */
    deltaTsubcycleNext_ = MIN(deltaT, deltaTsubcycleMax_);
    deltaTsubcycleNext_ = MIN(deltaTsubcycleNext_, deltaTsubcycle_init_init_);
    deltaTsubcycle_  = deltaTsubcycleNext_;

    //     Put a max of two to the growth of deltaT for the next global time step
    if (choiceDeltaTsubcycle_init_ == 0) {
      deltaTsubcycle_init_next_ =  deltaTsubcycle_;
    } else {
      deltaTsubcycle_init_next_ = 2.0*deltaT;
    }
    /*
     *  We calculate the number of estimated steps to take in the subcycling. This will be used to adjust the
     *  tolerance requirements needed for each of the steps
     */
    nsteps_est = deltaT * 0.99 / deltaTsubcycleNext_ + 1;
  
    /*
     *  Set the init and final state from the init_init state
     */
    setInitStateFromInitInit(true);
    /*
     *  Update the state
     */
    updateState();


  
    //
    // *****************************************************************************************************************
    //
    // Loop over any possible subcycling of the time step
    do {
 
      /*
       *  Increment the local subcycle counter
       */
      iterSubCycle++;
      /*
       *  Increment the grand counter, which can be used to isolate one bad integration amongst the thousands
       */
      counterNumberSubIntegrations_++;
      counterNumberLastSubIntegrations_++;
      /*
       *  Calculate the time interval from the previous step or from the initial value coming intot the routine
       */
      deltaTsubcycle_ = deltaTsubcycleNext_;
      /*
       *  Set the initial time interval
       */
      t_init_ = t_final_;
      /*
       * set the final time of the interval
       */
      t_final_ += deltaTsubcycleNext_;
      /*
       *  Bound the time interval by the global time step
       *      If we are close to the final step, put the final time at that step. 
       *      This avoids a case that occurred in practice where the final deltat after the current
       *      step was extremely small.
       */
      if (t_final_ > t_final_final_ - 0.1 * deltaTsubcycleNext_) {
	t_final_ = t_final_final_;
	deltaTsubcycle_ = t_final_ - t_init_;
	nsteps_est =  iterSubCycle;
      }
 
      /*
       *   Update initial values of the state vector with the final values from the last subcycle
       *   iteration, if we are beyond the first subscycle
       */
      if (iterSubCycle > 1) {
	setInitStateFromFinal(false);
      }

      /*
       *  Determine the big mole fractions of multispecies phases
       */
      determineBigMoleFractions();

      /*
       *  Check internal consistency of init solution
       */
      check_init_consistency();


      int nonlinConverged = 0;
      int attemptSubCycle = 0;
      while (nonlinConverged <= 0) {
      topConvergence:
	// increment the subcycle attempts
	attemptSubCycle++;
        if (iterSubCycle == 1) {
	if (attemptSubCycle > 40) {
	  throw CanteraError("InterfacialMassTransfer_Integrator::integrate()",
			     "FAILURE: Number attempted subcycles at first it is greater than 40. Something is wrong");
	}

        } else {
	if (attemptSubCycle > 10) {
	  throw CanteraError("InterfacialMassTransfer_Integrator::integrate()",
			     "FAILURE: Number attempted subcycles is greater than 10. Something is wrong");
	}
        }

	// Zero needed counters
	mdp::mdp_init_int_1(DATA_PTR(phaseJustDied_), 0, m_NumTotPhases);
	mdp::mdp_init_int_1(DATA_PTR(phaseJustBorn_), 0, m_NumTotPhases);

	/*
	 * Ok at this point we have a time step deltaTsubcycle_
	 * and initial conditions consisting of phaseMoles_init_ and spMF_init_.
	 * We now calculate predicted solution components from these conditions.
	 * Additionally, we evaluate whether any multispecies phases are going to pop into existence,
	 * adding in potential seed values and we set the phaseExistence flags in the kinetics solver
	 */
	int info = predictSoln();

	if (info != 1) {
	  conseqFailures++ ;
	  nonlinConverged = 0;
	  writeCSVData(-2);
	  t_final_ = t_init_;
	  deltaTsubcycle_ = deltaTsubcycle_ * 0.25;
	  deltaTsubcycleNext_ = deltaTsubcycle_;
	  deltaTsubcycleCalc_ = deltaTsubcycle_;
	  t_final_ += deltaTsubcycle_;
	  /*
	   * Revert to the old solution -> copy _init_ to _final_
	   */
	  setFinalStateFromInit();
	  goto topConvergence;
	}
	/*
	 *  Gather the predicted integrated src prediction
	 */
	gatherIntegratedSrcPrediction();
	
	nsteps_est =  iterSubCycle + (t_final_final_ - t_init_) / deltaTsubcycle_;
	if (nsteps_est > 100) {
	  nsteps_est = 100;
	}
	if (nsteps_est < 1) {
	  nsteps_est = 1;
	}

	double localRtolResidAllowed =  GlobalRtolSrcTerm/ nsteps_est;

	/*
	 *  We have a successful prediction: go and prep the nonlinear problem
	 */
	int num_newt_its = 0;
	int num_linear_solves = 0;
	int numBacktracks = 0;
	int loglevelInput = MAX(0, printLvl_ - 2);

	/*
	 *  formulate the nonlinear solver problem to be solved.
	 *     Fields to be filled in
	 *             yval
	 *             ylow
	 *             yhigh
	 *             yatol
	 *             ydeltaBoundsMagnitudes
	 *             ysType
	 *             
	 */
	initialPackSolver_nonlinFunction();
	
	/*
	 *
	 */
	check_nonlinResidConditions();
	int solnType = NSOLN_TYPE_STEADY_STATE;

	/*
	 *  
	 */
	pSolve_->setAtol(DATA_PTR(atolNLS_));
	pSolve_->setRtol(localRtolResidAllowed);
	/*
	 *  Set the tolerances on the residual evaluation of the nonlinear solver
	 *    We set the relative tolerance of the residual to GlobalRtolElectronSrcTerm/ nsteps_est. This is the
	 *    same value as the relative tolerance on the electrode stat variables.
	 *
	 *    Set the absolute tolerances for the residual calculation to a vector residAtolNLS_[] which is calculated
	 *    for the problem and is associated with the degrees of freedom in the problem. residAtolNLS_[] has
	 *    units of the residual. 
	 *
	 *    The residual tolerance is given by the minimum of the row weight sums and the user tolerances given below.
	 */
	pSolve_->setResidualTols(localRtolResidAllowed,  &residAtolNLS_[0]);

	pSolve_->setBoundsConstraints(&ylowNLS_[0], &yhighNLS_[0]);
	pSolve_->setDeltaBoundsMagnitudes(DATA_PTR(deltaBoundsMagnitudesNLS_));
	num_newt_its = 0;
	double CJ = 0.0;
	double time_curr = t_final_;
#ifdef DEBUG_MODE
	loglevelInput = 5;
	printLvl_ = 5;
	//	enableExtraPrinting_ = 3;
        //detailedResidPrintFlag_ = 10;
	//loglevelInput = 15;
        pSolve_->m_min_newt_its = 2;
#endif
	int nonlinearFlag = pSolve_->solve_nonlinear_problem(solnType, &yvalNLS_[0], &ydotNLS_[0], CJ,
							     time_curr, *jacPtr_,  num_newt_its, num_linear_solves,  numBacktracks, loglevelInput);
	if (nonlinearFlag < 0) {
	  if (printLvl_ > 2) {
	    printf("InterfacialMassTransfer_MultiPlateau_NoDiff::integrate(): Unsuccessful Nonlinear Solve flag = %d\n", nonlinearFlag);
	  }
	}
      

	if (nonlinearFlag < 0) {
	  conseqFailures++ ;
	  nonlinConverged = 0;
	  writeCSVData(-1);
	  t_final_ = t_init_;
	  deltaTsubcycle_ = deltaTsubcycle_ * 0.25;
	  deltaTsubcycleNext_ = deltaTsubcycle_;
	  t_final_ += deltaTsubcycle_;
	  /*
	   * Revert to the old solution -> copy _init_ to _final_
	   */
	  setFinalStateFromInit();
	} else {
	  unpackNonlinSolnVector(DATA_PTR(yvalNLS_));
	  /*
	   *  Correct the final time
	   */
	  if (deltaTsubcycleCalc_ < deltaTsubcycle_ * (1.0 - 1.0E-10)) {
	    if (printLvl_ > 1) {
	      printf("deltaT reduced to %g from %g for reason %d\n",
		     deltaTsubcycleCalc_,  deltaTsubcycle_, 0);
	    }
	    // Pick the next delta T to be equal to the current delta T
	    deltaTsubcycleNext_ = deltaTsubcycle_;
	    deltaTsubcycle_ = deltaTsubcycleCalc_;
	    t_final_ = t_init_ + deltaTsubcycle_;
	  } else  if (deltaTsubcycleCalc_ > deltaTsubcycle_ * (1.0 + 1.0E-10)) {
	    if (printLvl_ > 1) {
	      printf("deltaT increased to %g from %g for reason %d\n",
		     deltaTsubcycleCalc_,  deltaTsubcycle_, 1);
	    }
	    // Pick the next delta T to be equal to the current delta T
	    deltaTsubcycleNext_ = deltaTsubcycle_;
	    deltaTsubcycle_ = deltaTsubcycleCalc_;
	    t_final_ = t_init_ + deltaTsubcycle_;
	  } else {
	    /*
	     *  If we didn't fail this step, and if the nonlinear iterations are less than ~4,
	     *  increase the time step by a factor of 2. If we took more than 7 newton iterations
	     *  decrease the next time step
	     */
	    if (conseqFailures == 0) {
	      if (num_newt_its < 4) {
		deltaTsubcycleNext_ = MAX(deltaTsubcycleNext_, 2.0 * deltaTsubcycle_);
		deltaTsubcycleNext_ = MIN(deltaTsubcycleNext_, deltaTsubcycleMax_);
	      } else if (num_newt_its > 7) {
		if (deltaTsubcycle_ > relativeLocalToGlobalTimeStepMinimum_ * deltaT) {
		  deltaTsubcycleNext_ = MIN(deltaTsubcycleNext_, 0.5 * deltaTsubcycle_);
		}
	      }
	    }
	  }
	  conseqFailures--;
	  conseqFailures = MAX(0, conseqFailures);
	  nonlinConverged = 1;
	}
      }  // End of convergence of the nonlinear iteration for the current subcycle
      //
      //  ACTIVITIES THAT ARE NEEDED TO DETERMINE IF THE SUBCYCLE IS SUCCESSFUL EVEN IF NONLINEAR SOLVER CONVERGES
      //
      bool stepAcceptable = true;
      stepAcceptable = checkSubIntegrationStepAcceptable();
      /*
       *  Massage the results due to birth and death, before gauging the accuracy of the time stepping.
       */
      bool badCases = changeSolnForBirthDeaths();
      if (badCases) {
	stepAcceptable = false;
      }

      /*
       * Print out intermediate results if the print level is high enough
       */
      if (printLvl_ > 4) {
	printInterfacialMassTransfer(true, true);
      }

      
	
      //
      // ACTIVITIES THAT ARE NEEDED TO DETERMINE IF THE SUBCYCLE IS ACCURATE ENOUGH
      //
      if (stepAcceptable) {
      
	/*
	 *  Calculate the integrated source terms and do other items now that we have a completed time step
	 */
	calcSrcTermsOnCompletedStep();

	/*
	 *  Estimate the errors associated with the vector of source terms that are used by users of this object
	 */
	predictorCorrectorGlobalSrcTermErrorVector();

	/*
	 *  Now create a norm out of the error vector, making one number that will feed back into the time stepping
	 *  algorithm here.
	 */
	double pnormSrc = predictorCorrectorGlobalSrcTermErrorNorm();

	double pnormSoln = predictorCorrectorWeightedSolnNorm(yvalNLS_);

	pnorm = std::max(pnormSrc, pnormSoln);

	/*
	 * Adjust the next time step according to a predictor-corrector critera
	 */ 

	/*
	 * If we are way over the allowed local subcycle tolerance requirement, reject the step
	 */
	if (pnorm > 3.0) {
	  if (deltaTsubcycle_ > relativeLocalToGlobalTimeStepMinimum_ * deltaT) {
	    stepAcceptable = false;
	  }
	}

	if (printLvl_ > 4) {
	  predictorCorrectorPrint(yvalNLS_, pnormSrc, pnormSoln);
	}
      }

      /*
       *  HANDLE THE CASE WHERE WE HAVE NOT FOUND THE CURRENT STEP TO BE ACCEPTABLE
       *   - we decrease the time step and then redo the calculation
       */
      if (! stepAcceptable) {
#ifdef DEBUG_HKM_NOT
	printLvl_ = 9;
	enableExtraPrinting_ = 9;
	detailedResidPrintFlag_ = 9;
#endif
	conseqFailures++ ;
	nonlinConverged = 0;
	writeCSVData(-1);
	t_final_ = t_init_;
	deltaTsubcycle_ = deltaTsubcycle_ * 0.25;
	deltaTsubcycleNext_ = deltaTsubcycle_;
	t_final_ += deltaTsubcycle_;
	/*
	 * Revert to the old solution -> copy _init_ to _final_
	 */
	setFinalStateFromInit();


	goto topConvergence;
      } 
      /*
       *  Massage the results now that we have a successful step. This will mean zeroing out small phases
       */
      manageBirthDeathSuccessfulStep();

      /*
       * Accumulate local results into global vectors
       */ 
      accumulateSrcTermsOnCompletedStep();

      /*
       * Adjust the next time step according to a predictor-corrector critera
       */
 
      if (pnorm < 0.05) {
	if (deltaTsubcycleNext_ <= deltaTsubcycle_) {
	  // unusual case here
	  if (t_final_ < t_final_final_ - 1.0E-8) {
	    deltaTsubcycleNext_ = 1.5 * deltaTsubcycle_;
	  }
	}
      }
      if (deltaTsubcycle_ > relativeLocalToGlobalTimeStepMinimum_ * deltaT) {
	if (pnorm > 0.10) {
	  if (deltaTsubcycleNext_ > deltaTsubcycle_) {
	    deltaTsubcycleNext_ = deltaTsubcycle_;
	  }
	}
	if (pnorm > 0.20) {
	  if (deltaTsubcycleNext_ > 0.5 * deltaTsubcycle_) {
	    deltaTsubcycleNext_ = 0.5 * deltaTsubcycle_;
	    if (deltaTsubcycleNext_ < relativeLocalToGlobalTimeStepMinimum_ * deltaT) {
	      deltaTsubcycleNext_ = relativeLocalToGlobalTimeStepMinimum_ * deltaT;
	    }
	  }
	}
      } else {
	deltaTsubcycleNext_ = relativeLocalToGlobalTimeStepMinimum_ * deltaT;
      }


      /*
       *  If this is the first subcycle time step, calcualte the first time step for the next global iteration
       */
      if ((choiceDeltaTsubcycle_init_ == 1) && (iterSubCycle == 1)) {
	deltaTsubcycle_init_next_ = MIN(deltaTsubcycle_init_next_, deltaTsubcycleNext_);
      }


      if (t_final_ >= t_final_final_) {
	notDone = false;
      }

      if (printLvl_ > 1) {
	if (notDone) {
	  writeCSVData(1);
	} else {
	  writeCSVData(2);
	}
      }

      /*
       * Check final results
       */
      check_final_state();

      /*
       *  Calculate the first time step for the next global iteration
       */
      if (choiceDeltaTsubcycle_init_ == 0) {
	if (deltaTsubcycleNext_ < deltaTsubcycle_init_next_ ) {
	  // printf(" we are here %g, %g\n", deltaTsubcycleNext_ , deltaTsubcycle_init_next_);

	}
	deltaTsubcycle_init_next_ = deltaTsubcycleNext_;
	//  deltaTsubcycle_init_next_ = MIN(deltaTsubcycle_init_next_, deltaTsubcycleNext_);
      }


      /*
       *  Save the state of the object
       */


      //-------------------------------- End of Subcycle ---------------------------------------------------------------



    } while (notDone);


    /*
     *  Calcualte the first time step for the next global iteration
     */
    if (choiceDeltaTsubcycle_init_ == 0) {
      // deltaTsubcycle_init_next_ = MIN(deltaTsubcycle_init_next_, deltaTsubcycleNext_);
      deltaTsubcycle_init_next_ = MIN(deltaTsubcycle_init_next_, 2.0 * deltaT);
    }

    /*
     *  Copy the results into the final holding pens.
     */
    setFinalFinalStateFromFinal();



    return iterSubCycle;

  }



  //==================================================================================================================
  //   Calculate the largest mole fraction in each of the phases
  /* 
   *  We use this to determine the equation system
   */
  void InterfacialMassTransfer_Integrator::determineBigMoleFractions()
  {
  }
  //==================================================================================================================
  // Check the consistency of the initial state
  /*
   *  (virtual from InterfacialMassTransfer_Integrator) 
   *
   *   This is a checker routine. Therefore, it can be taken out in nondebug mode
   */
  void  InterfacialMassTransfer_Integrator::check_init_consistency()
  {

  }
  //==================================================================================================================
  // Set the base tolerances for the nonlinear solver within the integrator
  /*
   *   The tolerances are based on controlling the integrated source term 
   *   for molar production over the integration interval.  The integrated source term
   *   has units of kmol.
   *
   *   Because the electron is only one molar quantity within a bunch of molar quantities,
   *   this requirement will entail that we control the source terms of all species within the
   *   electrode to the tolerance requirements of the electron source term. 
   *
   *   @param rtolResid  Relative tolerance allowed for the electron source term over the interval.
   *                     This is a unitless quantity
   *   @param atolResid  Absolute tolerance cutoff for the electron source term over the interval.
   *                     Below this value we do not care about the results.
   *                     atol has units of kmol. 
   */
  void  InterfacialMassTransfer_Integrator::setNLSGlobalSrcTermTolerances(double rtolResid, double atolResid)
  {


    //! Relative tolerance for the integrated global src term vectors
    rtol_IntegratedSrc_global_ = rtolResid;

    //! Absolute tolerance for the integrated global src term vectors
    double sm = totalMoles();
    for (int i = 0; i <  numIntegratedSrc_; i++) {
      atol_IntegratedSrc_global_[i] = 1.0E-14 * sm;
    }


  } 
  //==================================================================================================================
  // Returns the number of integrated source terms whose errors are controlled by this integrator
  int  InterfacialMassTransfer_Integrator::numIntegratedSrcTerms() const
  {
    return numIntegratedSrc_;
  }
 //==================================================================================================================
  //   Set the Residual absolute error tolerances and the solution absolute error tolerances
  /* 
   *  (virtual from InterfacialMassTransfer_Integrator)
   *
   *   Set the absolute error tolerances fror the nonlinear solvers. This is called at the top
   *   of the integrator() routine.  We control both the solution error tolerances and also the
   *   residual error tolerances.
   *
   *                vector     residAtolNLS_(0),
   *                double     rtolResidNLS_(1.0E-6),
   *                vector     atolNLS_
   *                double     rtolNLS_
   */
  void  InterfacialMassTransfer_Integrator::setResidAtolNLS(double GlobalRtolSrcTerm)
  {
    residAtolNLS_[0] = 1.0E-14;
    atolNLS_[0] = 1.0E-14;
    //! Absolute tolerance for the integrated global src term vectors
    double sm = totalMoles();
    if (sm <= 1.0E-200) {
      throw CanteraError(" InterfacialMassTransfer_Integrator::setResidAtolNLS()",
			 "total moles is small or negative");
    }
    for (int i = 0; i <  neq_; i++) {
     residAtolNLS_[i] = 1.0E-14 * sm;
     atolNLS_[i] = 1.0E-14 * sm;
    }
    rtolResidNLS_ = MIN(1.0E-5, GlobalRtolSrcTerm);
    rtolNLS_ = MIN(1.0E-3, GlobalRtolSrcTerm);
  }
  //==================================================================================================================
  //  Predict the solution
  /* 
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
  int  InterfacialMassTransfer_Integrator::predictSoln()
  {
    throw CanteraError(" InterfacialMassTransfer_Integrator::predictSoln()","unimplemented");
  }
  //==================================================================================================================
  // Extract information from cantera's interface kinetics objects
  /* 
   *  (virtual fucntion from InterfacialMassTransfer_Integrator)
   *  In this routine we calculate the rates of progress of reactions and species on all active reacting surfaces.
   *
   *     ROP_RSD_List_[]
   */
  void InterfacialMassTransfer_Integrator::extractInfo() {
    /*
     * Loop over surface phases, filling in the phase existence fields within the
     * kinetics operator
     */
    for (int isk = 0; isk < numSurfaces_; isk++) {
      /*
       *   Only loop over surfaces that have kinetics objects associated with them
       */
      std::vector<double> & netROP =  ROP_RSD_List_[isk];
      size_t nr = netROP.size();
      if (nr > 0) {
	mdp::mdp_zero_dbl_1(DATA_PTR(netROP), nr);
      }
      if (ActiveKineticsSurf_[isk]) {   
	/*
         *  For each Reacting surface
         *
         *  Get the species production rates for the reacting surface
         */
        RSD_List_[isk]->getNetRatesOfProgress(DATA_PTR(netROP));
	/*
         *  Sometimes the ROP for goNowhere is nonzero, when it should be zero.
         */
        if (goNowhere_) {
          for (size_t i = 0; i < nr; i++) {
            netROP[i] = 0.0;
          }
        } 
      }
    }
  }
  //==================================================================================================================
  // Collect mole change information
  /*
   *  (virtual fucntion from InterfacialMassTransfer_Integrator)
   *
   *   We take the ROP_inner_[] and ROP_outer_[] rates of progress, combine them with the surface area calculation,
   *   and the stoichiometric coefficients to calculate the DspMoles_final_[], which is production rate for
   *   all species.
   */
  void InterfacialMassTransfer_Integrator::updateSpeciesMoleChangeFinal()
  {
 
    mdp::mdp_zero_dbl_1(DATA_PTR(DspMoles_final_), m_NumTotSpecies);

    /*
     * Loop over surface phases, filling in the phase existence fields within the
     * kinetics operator
     */
    for (int isk = 0; isk < numSurfaces_; isk++) {
      /*
       *   Only loop over surfaces that have kinetics objects associated with them
       */
      std::vector<double> & netROP =  ROP_RSD_List_[isk];
      size_t nr = netROP.size();
      if (nr > 0) {

	if (ActiveKineticsSurf_[isk]) {
	  double mult =  0.5* (surfaceAreaRS_init_[isk] + surfaceAreaRS_final_[isk]); 
	  for (int i = 0; i < m_NumTotSpecies; i++){
	    for (int j = 0; j < numRxns_[isk]; j++){
	      DspMoles_final_[i] +=  mult * productStoichCoeff(isk,i,j) * netROP[j];
	      DspMoles_final_[i] -=  mult * reactantStoichCoeff(isk,i,j) * netROP[j];
	    }
	  }
	}
      }
    }
  
  }
  //====================================================================================================================
  // Update the velocities within the object at the final conditions.
  /*
   *  (virtual from InterfacialMassTransfer)
   * 
   *  In order to do this, we need the instantaenous source terms and the current value of the
   *  surface area over the local time step. Therefore, this routine must be called 
   *  after DspMoles_final_ has been calculated. 
   */
  void InterfacialMassTransfer_Integrator::updateVelocities()
  {
    double massfluxA = getPhaseAMassSourceTerm();
    ThermoPhase &tpA = thermo(solnAPhase_);
    double densA = tpA.density();
    double area = 0.5 * (surfaceAreaRS_init_[0] + surfaceAreaRS_final_[0]);
    double sVelocA = massfluxA  / area / densA;

    double massfluxB = getPhaseBMassSourceTerm();
     ThermoPhase &tpB = thermo(solnBPhase_);
    double densB = tpB.density();
    double sVelocB = massfluxB  / area / densB;

    Velo_A_Interface_final_ = Velo_ReferenceFrame_;
    Velo_S_Interface_final_ = sVelocA + Velo_A_Interface_final_;
    Velo_B_Interface_final_ = sVelocB + Velo_S_Interface_final_;
  }
  
  //==================================================================================================================
  // calculate the residual
  /*
   *   
   *  (virtual fucntion from InterfacialMassTransfer_Integrator)
   *
   */
  int InterfacialMassTransfer_Integrator::calcResid(doublereal * const resid, const ResidEval_Type_Enum evalType)
  {
    throw CanteraError(" InterfacialMassTransfer_Integrator::calcResid()",  "unimplemented");
  }
  //==================================================================================================================
  //  Gather the predicted solution values and the predicted integrated source terms
  /*
   *  (virtual from InterfacialMassTransfer_Integrator)
   *
   *  Both the predicted soln and the predicted integrated source terms are used
   *  in the time step control.
   *
   *  Here we derive a value for the predicted source term for the local step, 
   *  IntegratedSrc_Predicted[]
   */
  void InterfacialMassTransfer_Integrator::gatherIntegratedSrcPrediction()
  {
    /*
     *  We have a new predicted solution vector at the new time. Update the state with this new solution vector
     */
    updateState();

    /*
     *  Extract rates of progress from Cantera's interfacial kinetics objects
     */
    extractInfo();
 
    /*
     *  Update the values for  DspMoles_final_
     */
    updateSpeciesMoleChangeFinal();

    /*
     *  Update the Stefan Velocity internal variables for the current final conditions
     */
    updateVelocities();

    /*
     *  Calculate the integrated source term
     */
    for (int i = 0; i < m_NumTotSpecies; i++) {
      IntegratedSrc_Predicted[i] = DspMoles_final_[i] * deltaTsubcycle_;
    }

  }
  //==================================================================================================================
  //   Calculate the integrated source terms and do other items now that we have a completed time step
  /*
   *  Calculate source terms on completion of a step. At this point we have solved the nonlinear problem
   *  for the current step, and we are calculating post-processed quantities like source terms.
   */
 void InterfacialMassTransfer_Integrator::calcSrcTermsOnCompletedStep()
  {
    /*
     *  We have a new predicted solution vector at the new time. Update the state with this new solution vector
     */
    updateState();

    /*
     *  Extract rates of progress from Cantera's interfacial kinetics objects
     */
    extractInfo();
 
    /*
     *  Update the values for  DspMoles_final_
     */
    updateSpeciesMoleChangeFinal();

    /*
     *  Update the Stefan Velocity internal variables for the current final conditions
     */
    updateVelocities();

    /*
     *  Calculate the integrated source term
     */
    for (int i = 0; i < m_NumTotSpecies; i++) {
      spMoleIntegratedSourceTermLast_[i] = DspMoles_final_[i] * deltaTsubcycle_;
    }
  }
  //==================================================================================================================     
  //   Accumulate src terms and other results from the local step into the global holding bins.
  /*
   *  Accumulate source terms on completion of a step. At this point we have solved the nonlinear problem
   *  for the current step and we have satisfied all accuracy requirements.
   *  The step is good. We now accumulate the results before going on to a new local step.
   */
  void InterfacialMassTransfer_Integrator::accumulateSrcTermsOnCompletedStep()
    {
      for (int i = 0; i < m_NumTotSpecies; i++) {
	spMoleIntegratedSourceTerm_[i] += spMoleIntegratedSourceTermLast_[i]; 
      }
    }
  //==================================================================================================================
  // Set the internal final intermediate and from the internal init state
  /*
   *  (non-virtual function)  -> function should onionize in-first.
   *
   *  Set the final state from the init state. This is commonly called during a failed time step
   *
   */
  void InterfacialMassTransfer_Integrator::setFinalStateFromInit_Oin()
  {
    InterfacialMassTransfer::setFinalStateFromInit_Oin();
  } 
  //==================================================================================================================
  // Set the internal final intermediate and from the internal init state
  /*
   *  (virtual function from InterfacialMassTransfer) 
   *
   *  Set the final state from the init state. This is commonly called during a failed time step
   *
   */
  void InterfacialMassTransfer_Integrator::setFinalStateFromInit()
  {
    InterfacialMassTransfer_Integrator::setFinalStateFromInit_Oin();
  }
  //==================================================================================================================
  // Return the number of equations in the equation system that is used to solve the ODE integration
  int InterfacialMassTransfer_Integrator::nEquations() const
  {
    return neq_;
  }
  //==================================================================================================================
  //  Return a vector of delta y's for calculation of the numerical Jacobian 
  /*  
   * (virtual from ResidJacEval)
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
  int InterfacialMassTransfer_Integrator::calcDeltaSolnVariables(const doublereal t, const doublereal * const ySoln,
						   const doublereal * const ySolnDot, doublereal * const deltaYSoln,
						   const doublereal * const solnWeights)
  {
    int retn = ResidJacEval::calcDeltaSolnVariables(t, ySoln, ySolnDot,deltaYSoln,solnWeights);
    return retn;
  }
  //==================================================================================================================
  //   Unpack the soln vector
  /* 
   *  (virtual from InterfacialMassTransfer_Integrator)
   *
   *  This function unpacks the solution vector into  phaseMoles_final_,  spMoles_final_, and spMf_final_[]
   */
  void InterfacialMassTransfer_Integrator::unpackNonlinSolnVector(const double * const y)
  {  
    throw CanteraError(" InterfacialMassTransfer_Integrator::unpackNonlinSolnVector()", "unimplemented");
  }
  //====================================================================================================================
  // Check to see that the preceding step is a successful one
  /*
   *   We check to see if the preceding step is a successful one.
   *
   *  @return Returns a bool true if the step is acceptable, and false if it is unacceptable.
   */
  bool InterfacialMassTransfer_Integrator::checkSubIntegrationStepAcceptable() const
  {
    return true;
  } 
  //=================================================================================================================
  double l0norm2M(const std::vector<double> &v1,const std::vector<double> &v2, int num, double rtol,
		const std::vector<double> &atolVec)
  {
    double max0 = 0.0;
    double denom, diff, ee;
    
    for (int k = 0; k < num; k++) {
      diff =  fabs(v1[k] - v2[k]);
      denom = 0.5 * (fabs(v1[k]) + fabs(v2[k])) * rtol + atolVec[k];
      ee = fabs( diff/ denom);
      if (ee > max0) {
        max0 = ee;
      }
    }
    return max0;
  }
  //====================================================================================================================
  //  Calculate the norm of the difference between the predicted answer and the final converged answer 
  //  for the current time step
  /*
   *   The norm calculated by this routine is used to determine whether the time step is accurate enough.
   *
   *  @return    Returns the norm of the difference. Normally this is the L2 norm of the difference
   */  
  double InterfacialMassTransfer_Integrator::predictorCorrectorWeightedSolnNorm(const std::vector<double> & yvalNLS)
  {  
    double nn = l0norm2M(soln_predict_, yvalNLS_, neq_, rtolNLS_, atolNLS_);
    return nn;
  }
  //====================================================================================================================
  // Calculate the vector of predicted errors in the source terms that this integrator is responsible for
  /*
   *  (virtual from InterfacialMassTransfer_Integrator)
   *
   *    In the base implementation we assume that the there are just one source term, the electron
   *    source term.
   *    However, this will be wrong in almost all cases. 
   *    The number of source terms is unrelated to the number of unknowns in the nonlinear problem.
   *    Source terms will have units associated with them. 
   *    For example the integrated source term for electrons will have units of kmol 
   */
  void InterfacialMassTransfer_Integrator::predictorCorrectorGlobalSrcTermErrorVector()
  {
    for (int k = 0; k < m_NumTotSpecies; k++) {
      IntegratedSrc_Errors_local_[k] = IntegratedSrc_Predicted[k] - spMoleIntegratedSourceTermLast_[k];
    }
  }
 
  //====================================================================================================================
  //  Calculate the normalized norm of the errors in the global source terms 
  /*
   *  (virtual from InterfacialMassTransfer_Integrator)
   *
   *   This routine make use of the source term error vector along with rtols and atols for the 
   *   individual source terms to calculated a normalized error measure. This is the single number
   *   that the integration routine will try to control as it calculates a time stepping strategy.
   *    
   *   @return  Returns a single nondimensional number representing the normalized error
   *            for the calculation of the source term
   */
  double InterfacialMassTransfer_Integrator::predictorCorrectorGlobalSrcTermErrorNorm()
  {
    double nn = l0norm2M(IntegratedSrc_Predicted, spMoleIntegratedSourceTermLast_, m_NumTotSpecies,
			 rtol_IntegratedSrc_global_, atol_IntegratedSrc_global_);
    return nn;
  }
  //====================================================================================================================
  // Print table representing prediction vs. corrector information
  /*
   *  @param yval           Vector of corrector values
   *  @param pnormSrc       Norm of the predictor-corrector comparison for the source vector.
   *  @param pnormSoln      Norm of the predictor-corrector comparison for the solution vector.
   */
  void  InterfacialMassTransfer_Integrator::predictorCorrectorPrint(const std::vector<double> &yval, 
						      double pnormSrc, double pnormSoln)
  {

    double denom;
    double tmp;
    printf(" -------------------------------------------------------------------------------------------------------------------\n");
    printf(" PREDICTOR_CORRECTOR  SubIntegrationCounter = %7d       t_init = %12.5E,       t_final = %12.5E\n",
	   counterNumberSubIntegrations_, t_init_, t_final_);
    printf("                         IntegrationCounter = %7d  t_init_init = %12.5E, t_final_final = %12.5E\n",
	   counterNumberIntegrations_, t_init_init_, t_final_final_);
    printf(" -------------------------------------------------------------------------------------------------------------------\n");
    printf("                            Initial        Prediction     Actual           Difference       Tol          Contrib   |\n");
  
    denom = MAX(fabs(yval[0]), fabs(soln_predict_[0]));
    denom = MAX(denom, atolNLS_[0]);
    tmp = fabs((yval[0] - soln_predict_[0])/ denom);
    printf(" DeltaT                  | %14.7E %14.7E %14.7E | %14.7E | %10.3E | %10.3E |\n",
	   deltaTsubcycle_, soln_predict_[0],  yval[0], yval[0] - soln_predict_[0], atolNLS_[0], tmp);
    for (int i = 1; i < (int) yval.size(); i++) {
      denom = 0.5 * (fabs(yval[i]) + fabs(soln_predict_[i])) * rtolNLS_ + atolNLS_[i];
      tmp = fabs((yval[i] - soln_predict_[i])/ denom);
      printf(" soln %3d                |                %14.7E %14.7E | %14.7E | %10.3E | %10.3E | \n",
	     i, soln_predict_[i],  yval[i], yval[i] - soln_predict_[i], atolNLS_[i], tmp);
    }
    printf(" -------------------------------------------------------------------------------------------------------------------\n");
    printf("                                                                                                        %10.3E\n",
	   pnormSoln);
  }
  //====================================================================================================================
  // Possibly change the solution due to phase births and deaths. 
  /*
   *   (virtual from InterfacialMassTransfer_Integrator)
   *
   */
  bool InterfacialMassTransfer_Integrator::changeSolnForBirthDeaths()
  {
    return false;
  } 
  //====================================================================================================================
  // Possibly change the solution due to phase births and deaths. 
  /*
   *   (virtual from InterfacialMassTransfer_Integrator)
   *   
   *  This routine is carried out after the step is deemed a success. Massaging of the solution
   *  must be carried out within strict tolerances.
   */
  void InterfacialMassTransfer_Integrator::manageBirthDeathSuccessfulStep()
  {
  } 
  //====================================================================================================================
  //   Error check on the routine step
  /*
   *    (virtual from InterfacialMassTransfer_Integrator)
   *
   *   Error checks go here. All errors are fatal exits.
   */
  void InterfacialMassTransfer_Integrator::check_final_state()
  {
  } 
  //====================================================================================================================
  // Print a header for the residual printout
  /*
   *  (virtual from Eelctrode_Integrator)
   */
  void InterfacialMassTransfer_Integrator::printResid_TimeHeader()
  {
  } 
  //====================================================================================================================
  //   Check for problems with the residual
  /*
   *  (virtual from InterfacialMassTransfer_Integrator)
   *
   *  Checks here will cause the current nonlinear solve to fail
   *
   *    @return Return a negative value if there is a fatal problem. 
   */
  int  InterfacialMassTransfer_Integrator::residEval_BaseChecks()
  {
    return 1;
  } 



  // ----------------------------------------------------------------------------------------------
  // ----------------------------- GET CONDITIONS OUT --------------------------------------------
  // ----------------------------------------------------------------------------------------------
  //
  //       (unless specified this is always at the final conditions and time
  //
  
  // ----------------------------- GET INSTANTANEOUS SOURCE TERMS --------------------------------

  // Get the net production rates of all species in the electrode object
  // at the current conditions
  /*
   *
   *  This routine assumes that the underlying objects have been updated.
   *  It uses the default Backwards Euler integration rule, which doesn't take into account of issues with surfaces going away
   */
  void InterfacialMassTransfer_Integrator::getNetProductionRates(doublereal* const net) const
  {
    mdp::mdp_zero_dbl_1(net, m_NumTotSpecies);
    /*
     *  This routine basically translates between species lists for the reacting surface
     *  domain and the InterfacialMassTransfer.
     */
    /*
     *  For each Reacting surface
     */ 
    for (int isk = 0; isk < numSurfaces_; isk++) {
      if (ActiveKineticsSurf_[isk]) {
	/*
	 *  Just assume surface area is equal to final value
	 */
	double area = 0.5 * (surfaceAreaRS_init_[isk] + surfaceAreaRS_final_[isk]);
	/*
	 *  Get the species production rates for the reacting surface
	 */
	const vector<double> &rsSpeciesProductionRates = RSD_List_[isk]->calcNetProductionRates();
	/*
	 * Number of phases in the reacting species.
	 */
	int nphRS = RSD_List_[isk]->nPhases();

	int kIndexKin = 0;
	/*
	 *  Loop over the phases defined in the Reacting Surface Domain object
	 */
	for (int kph = 0; kph < nphRS; kph++) {
	  /*
	   *  Find the phase Id of the current phase in the InterfacialMassTransfer object
	   */
	  int jph = RSD_List_[isk]->kinOrder[kph];

	  /*
	   *  Find the starting species index within the InterfacialMassTransfer object
	   */
	  int istart = m_PhaseSpeciesStartIndex[jph];
	  int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
	  for (int k = 0; k < nsp; k++) {
	    net[istart + k] += rsSpeciesProductionRates[kIndexKin] * area;
	    kIndexKin++;
	  }
	}
      }
    }
  } 
 //====================================================================================================================
  // Get the net production rates of all species in the electrode object
  // at the current conditions
  /*
   *
   *  This routine assumes that the underlying objects have been updated.
   *  It uses the default Backwards Euler integration rule, which doesn't take into account of issues with surfaces going away
   */
  void InterfacialMassTransfer_Integrator::getNetProductionRatesRSD(const int isk, doublereal* const net) const
  {
    mdp::mdp_zero_dbl_1(net, m_NumTotSpecies);
    /*
     *  This routine basically translates between species lists for the reacting surface
     *  domain and the InterfacialMassTransfer.
     */
    /*
     *  For each Reacting surface
     */
    if (ActiveKineticsSurf_[isk]) {
      /*
       *  Just assume surface area is equal to final value
       */
      double area = 0.5 * (surfaceAreaRS_init_[isk] + surfaceAreaRS_final_[isk]);
      /*
       *  Get the species production rates for the reacting surface
       */
      const vector<double> &rsSpeciesProductionRates = RSD_List_[isk]->calcNetProductionRates();
      /*
       * Number of phases in the reacting species.
       */
      int nphRS = RSD_List_[isk]->nPhases();

      int kIndexKin = 0;
      /*
       *  Loop over the phases defined in the Reacting Surface Domain object
       */
      for (int kph = 0; kph < nphRS; kph++) {
	/*
	 *  Find the phase Id of the current phase in the InterfacialMassTransfer object
	 */
	int jph = RSD_List_[isk]->kinOrder[kph];

	/*
	 *  Find the starting species index within the InterfacialMassTransfer object
	 */
	int istart = m_PhaseSpeciesStartIndex[jph];
	int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
	for (int k = 0; k < nsp; k++) {
	  net[istart + k] += rsSpeciesProductionRates[kIndexKin] * area;
	  kIndexKin++;
	}
      }
    }
  }
  //====================================================================================================================
  //  Returns the current and the net production rates of the phases in kg/s from a single surface
  /*
   *  Returns the net production rates of all phases from reactions on a single surface
   *  
   *  @param isk Surface ID to get the fluxes from.      
   *  @param phaseMassFlux  Returns the mass fluxes of the phases
   */
  void InterfacialMassTransfer_Integrator::getPhaseMassFlux(doublereal* const phaseMassFlux) const
  {
    mdp::mdp_zero_dbl_1(phaseMassFlux, m_NumTotPhases);

    for (int isk = 0; isk < numSurfaces_; isk++) {
      if (ActiveKineticsSurf_[isk]) { 
	const vector<double> &rsSpeciesProductionRates = RSD_List_[isk]->calcNetProductionRates();
	double area =  0.5 * (surfaceAreaRS_init_[isk] + surfaceAreaRS_final_[isk]);
	int nphRS = RSD_List_[isk]->nPhases();
	int kIndexKin = 0;
	for (int kph = 0; kph < nphRS; kph++) {
	  int jph = RSD_List_[isk]->kinOrder[kph];
	  ThermoPhase &tp = thermo(jph);
	  int istart = m_PhaseSpeciesStartIndex[jph];
	  int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
	  for (int k = 0; k < nsp; k++) {
	    double net = rsSpeciesProductionRates[kIndexKin] * area;
	    double mw = tp.molecularWeight(k);
	    phaseMassFlux[jph] += net * mw;
	    kIndexKin++;
	  }
	}
      }
    }
  }
  //================================================================================================
  //  Returns the current and the net production rates of the phases in kmol/s from a single surface
  /*
   *  Returns the net production rates of all phases from reactions on a single surface
   *  
   *  @param isk Surface ID to get the fluxes from.      
   *  @param phaseMassFlux  Returns the mass fluxes of the phases
   */
  void InterfacialMassTransfer_Integrator::getPhaseMoleFlux(const int isk, doublereal* const phaseMoleFlux) const
  { 
    mdp::mdp_zero_dbl_1(phaseMoleFlux, m_NumTotPhases);
    for (int isk = 0; isk < numSurfaces_; isk++) {
      if (ActiveKineticsSurf_[isk]) {
	const vector<double> &rsSpeciesProductionRates = RSD_List_[isk]->calcNetProductionRates();
	double area = 0.5 * (surfaceAreaRS_init_[isk] + surfaceAreaRS_final_[isk]);
	int nphRS = RSD_List_[isk]->nPhases();
	int PLph, RSph;
	int kIndexKin = 0;
	for (RSph = 0; RSph < nphRS; RSph++) {
	  PLph = RSD_List_[isk]->kinOrder[RSph];
	  int istart = m_PhaseSpeciesStartIndex[PLph];
	  int nsp = m_PhaseSpeciesStartIndex[PLph+1] - istart;
	  for (int k = 0; k < nsp; k++) {
	    double net = rsSpeciesProductionRates[kIndexKin]* area;
	    phaseMoleFlux[PLph] += net;
	    kIndexKin++;
	  }
	}
      }
    }
  }
  //====================================================================================================================
  // Returns the computed Stefan velocity for phase A in the object
  // at the current final conditions. 
  /*
   *  Multiply by the surface area to get the stefan volume production rate at the
   *  current final conditions.
   * 
   *    @return returns stefan velocity created in phase A (m s-1)
   */
  double InterfacialMassTransfer_Integrator::StefanVelocityPhaseA() const
  {
    double massflux = getPhaseAMassSourceTerm();
    ThermoPhase &tpA = thermo(solnAPhase_);
    double dens = tpA.density();
    double area = 0.5 * (surfaceAreaRS_init_[0] + surfaceAreaRS_final_[0]);
    return massflux / area / dens;
  }
  //====================================================================================================================
  //  Returns the computed Stefan velocity for phase B in the object
  //  at the current final conditions.
  /* 
   *  Multiply by the surface area to get the stefan volume production rate at the
   *  current final conditions.
   * 
   *    @return returns stefan velocity created in phase B (m s-1)
   */
  double InterfacialMassTransfer_Integrator::StefanVelocityPhaseB() const
  {
    double massflux = getPhaseBMassSourceTerm();
    ThermoPhase &tpB = thermo(solnBPhase_);
    double dens = tpB.density(); 
    double area = 0.5 * (surfaceAreaRS_init_[0] + surfaceAreaRS_final_[0]);
    return massflux / area / dens;
  }
  //====================================================================================================================
  //     Print details about the satisfaction of the residual
  /* 
   *  (virtual from InterfacialMassTransfer_Integrator)
   */
  void  InterfacialMassTransfer_Integrator::printResid_ResidSatisfaction()
  {
  }
  //====================================================================================================================
  //  Pack the nonlinear solver proplem
  /*
   *  formulate the nonlinear solver problem to be solved.
   *     Fields to be filled in
   *     Fields to be filled in
   *             yvalNLS_
   *             ylowNLS_
   *             yhighNLS_
   *             atolNLS_
   *             deltaBoundsMagnitudesNLS_     
   */
  void InterfacialMassTransfer_Integrator::initialPackSolver_nonlinFunction()
  {
    throw CanteraError("InterfacialMassTransfer_Integrator::initialPackSolver_nonlinFunction()",
		       "unimplemented");
  }
  //====================================================================================================================
  // Check the nonlinear residual equations for completeness and the ability to be solved
  /*
   *
   */
  void InterfacialMassTransfer_Integrator::check_nonlinResidConditions()
  { 
  }

  //====================================================================================================================
  //  Residual calculation for the solution of the Nonlinear integration problem
  /*
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
   */
  int InterfacialMassTransfer_Integrator::evalResidNJ(const doublereal t, const doublereal delta_t,
					const doublereal * const y, const doublereal * const ySolnDot,
					doublereal * const resid,
					const ResidEval_Type_Enum evalType, const int id_x,
					const doublereal delta_x)
  {
    int retn = 1;
    if ((evalType == Base_ShowSolution) || (enableExtraPrinting_ && detailedResidPrintFlag_ > 1)) {
      printf("\t\t===============================================================================================================================\n");
      printf("\t\t  EXTRA PRINTING FROM NONLINEAR RESIDUAL: ");
      if (evalType ==  Base_ResidEval) {
        printf(" BASE RESIDUAL");
      } else if (evalType == JacBase_ResidEval) {
        printf(" BASE JAC RESIDUAL");
      } else  if  (evalType == JacDelta_ResidEval) {
        printf(" DELTA JAC RESIDUAL");
        printf(" var = %d delta_x = %12.4e Y_del = %12.4e Y_base = %12.4e", id_x, delta_x, y[id_x], y[id_x] - delta_x);
      } else  if  (evalType == Base_ShowSolution) {
        printf(" BASE RESIDUAL - SHOW SOLUTION");
      }
      printf(" DomainNumber = %2d , CellNumber = %2d , IntegrationCounter = %d\n",
             DomainNumber_, CellNumber_, counterNumberIntegrations_);
    }
    if ((evalType != JacDelta_ResidEval)) {
      mdp::mdp_init_int_1(DATA_PTR(phaseJustDied_), 0, m_NumTotPhases);
    }
    /*
     *  UNPACK THE SOLUTION VECTOR
     */
    unpackNonlinSolnVector(y);

    if (evalType != JacDelta_ResidEval && (evalType != Base_LaggedSolutionComponents)) {
      //    mdp::mdp_copy_dbl_1(DATA_PTR(phaseMoles_final_lagged_),(const double *)DATA_PTR(phaseMoles_final_), m_NumTotPhases);
    }


    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
      printResid_TimeHeader();
 
    }
 
    /*
     *  Take the unpacked solution and calculate the consistent and full final state
     */
    updateState();

    /*
     *  Query Cantera for all of the rate information at the final state (and the initial state if we are doing higher order)
     */
    extractInfo();

    /*
     *   We take the ROP_inner_[] and ROP_outer_[] rates of progress, combine them with the surface area calculation,
     *   and the stoichiometric coefficients to calculate the DspMoles_final_[], which is production rate for
     *   all species in the electrode.
     */
    updateSpeciesMoleChangeFinal();

    /*
     *  Update the Stefan Velocity internal variables for the current final conditions
     */
    updateVelocities();

    /*
     * Calculate the residual
     */
    calcResid(resid, evalType);

    /*
     * Change the problem specification for the nonlinear solve in certain circumstances when the solver
     *  is calculating the base residual of a jacobian 
     */
    if (evalType == JacBase_ResidEval) {
      /*
       *  Perform some checks that would lead to the return flag indicating an error condition
       */
      retn = residEval_BaseChecks();
    }

 
    if ((evalType == Base_ShowSolution) || (enableExtraPrinting_ && detailedResidPrintFlag_ > 1)) {
      printResid_ResidSatisfaction();

    }
    int rr = 1;
    if (retn < 0) {
      rr = retn;
    }
    return rr;
  }
  //====================================================================================================================
  // Fill in the initial conditions
  /* 
   * (virtual from NonlinearSolver)
   *
   * Values for both the solution and the value of ydot may be provided.
   *
   * @param t0            Time                    (input) 
   * @param y             Solution vector (output)
   * @param ydot          Rate of change of solution vector. (output)
   */
  int InterfacialMassTransfer_Integrator::getInitialConditions(const doublereal t0, doublereal * const y, doublereal * const ydot)
  {
    return 0;
  }
  //====================================================================================================================
  void InterfacialMassTransfer_Integrator::printInterfacialMassTransfer(int pSrc, bool subTimeStep) {
    int iph;
    double *netROP = new double[m_NumTotSpecies];
    double egv = TotalVol();
    printf("   ===============================================================\n");
    if (subTimeStep) {
      printf("      InterfacialMassTransfer at intermediate-step time final = %g\n", t_final_);
      printf("                   intermediate-step time init  = %g\n", t_init_);
    } else {
      printf("      InterfacialMassTransfer at time final = %g\n", t_final_final_);
      printf("                   time init  = %g\n", t_init_init_);
    }
    printf("   ===============================================================\n");
    printf("          Number of surfaces = %d\n", numSurfaces_);
    printf("          Total Volume = %10.3E\n", egv);
    printf("          Temperature = %g\n", Temp_);
    printf("          Pressure A = %g\n", Pres_A_Interface_final_);
    printf("          Pressure B = %g\n", Pres_B_Interface_final_);


    for (iph = 0; iph < m_NumTotPhases; iph++) {
      printInterfacialMassTransferPhase(iph, pSrc);
      printf("     ===============================================================\n");
    }
    delete [] netROP;
  }
  //===================================================================================================================
 
  void InterfacialMassTransfer_Integrator::printInterfacialMassTransferPhase(int iph, int pSrc, bool subTimeStep) {
    int isph;
    double *netROP = new double[m_NumTotSpecies];
    ThermoPhase &tp = thermo(iph);
    string pname = tp.id();
    int istart = m_PhaseSpeciesStartIndex[iph];
    int nsp = tp.nSpecies();
    printf("     ===============================================================\n");
    printf("          Phase %d %s \n", iph,pname.c_str() );
    printf("                Total moles = %g\n", phaseMoles_final_[iph]);
 
 
    printf("                  Voltage = %g\n", tp.electricPotential());
    
    if (iph >= NumVolPhases_) {
      isph = iph - NumVolPhases_;
      printf("                surface area (final) = %g\n",  surfaceAreaRS_final_[isph]);
      printf("                surface area (init)  = %g\n",  surfaceAreaRS_init_[isph]);
    }
    printf("\n");
    printf("                Name               MoleFrac_final  kMoles_final kMoles_init SrcTermLastStep(kMoles)\n");
    for (int k = 0; k < nsp; k++) {
      string sname = tp.speciesName(k);
      if (pSrc) {
	if (subTimeStep) { 
	  printf("                %-22s %10.3E %10.3E   %10.3E  %10.3E\n", sname.c_str(), spMf_final_[istart + k],
		 spMoles_final_[istart + k], spMoles_init_[istart + k],
		 spMoleIntegratedSourceTermLast_[istart + k]);
	} else {
	  printf("                %-22s %10.3E %10.3E   %10.3E  %10.3E\n", sname.c_str(), spMf_final_[istart + k],
	       spMoles_final_[istart + k], spMoles_init_init_[istart + k],
	       spMoleIntegratedSourceTerm_[istart + k]);
	}
      } else {
	if (subTimeStep) { 
	  printf("                %-22s %10.3E %10.3E   %10.3E\n", sname.c_str(), spMf_final_[istart + k], 
		 spMoles_final_[istart + k],   spMoles_init_[istart + k]);
	} else {
	  printf("                %-22s %10.3E %10.3E   %10.3E\n", sname.c_str(), spMf_final_[istart + k], 
		 spMoles_final_[istart + k],   spMoles_init_init_[istart + k]);
	}
      }
    }
    if (iph >= NumVolPhases_) {
      const vector<double> &rsSpeciesProductionRates = RSD_List_[isph]->calcNetProductionRates();
      RSD_List_[isph]->getNetRatesOfProgress(netROP);
      
      doublereal * spNetProdPerArea = (doublereal *) spNetProdPerArea_List_.ptrColumn(isph);
      mdp::mdp_zero_dbl_1(spNetProdPerArea, m_NumTotSpecies);
      int nphRS = RSD_List_[isph]->nPhases();
      int kIndexKin = 0;
      for (int kph = 0; kph < nphRS; kph++) {
	int jph = RSD_List_[isph]->kinOrder[kph];
	int istart = m_PhaseSpeciesStartIndex[jph];
	int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
	for (int k = 0; k < nsp; k++) {
	  spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
	  kIndexKin++;
	}
      }
      printf("\n");
      printf("                           spName                  Source (kmol/m2/s) \n");
      for (int k = 0; k <  m_NumTotSpecies; k++) {
	string ss = speciesName(k);
	printf("                           %-22s %10.3E\n", ss.c_str(), spNetProdPerArea[k]);
      }
    }
    delete [] netROP;

  }


} // End of namespace Cantera
//======================================================================================================================
