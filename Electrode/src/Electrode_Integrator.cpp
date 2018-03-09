/**
 *  @file Electrode_Integrator.cpp
 *     Definitions of the Electrode_Integrator class, used to perform subtimestep integrations on top of the Electrode class
 *     (see \ref electrode_mgr and class \link Zuzax::Electrode_Integrator Electrode_Integrator\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "Electrode_Integrator.h"
#include "cantera/numerics/NonlinearSolver_JAC.h"

#ifndef SAFE_DELETE
//! Define a delete and set to null operation
#define SAFE_DELETE(x)  if (x) { delete x;  x = nullptr;}
#endif

#ifndef MAX
//! define a quick max op
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif

#ifndef MIN
//! define a quick min op
#define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif

// #define DEBUG_CHECK_XML

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
SubIntegrationHistory::SubIntegrationHistory() :
    nTimeStepsRegular_(0),
    nTimeSteps_(0),
    TimeStepList_(0),
    GsolnErrorNorm_(0.0),
    GwdotErrorNorm_(0.0),
    iCounter(0),
    time_step_next(0.0)
{
}
//==================================================================================================================================
SubIntegrationHistory::SubIntegrationHistory(const SubIntegrationHistory& right) :
    nTimeStepsRegular_(right.nTimeStepsRegular_),
    nTimeSteps_(right.nTimeSteps_),
    TimeStepList_(right.TimeStepList_),
    GsolnErrorNorm_(right.GsolnErrorNorm_),
    GwdotErrorNorm_(right.GwdotErrorNorm_),
    iCounter(right.iCounter),
    time_step_next(right.time_step_next)
{
}
//==================================================================================================================================
SubIntegrationHistory::~SubIntegrationHistory()
{
}
//==================================================================================================================================
SubIntegrationHistory& SubIntegrationHistory::operator=(const SubIntegrationHistory& right)
{
    if (this == &right) {
        return *this;
    }
    nTimeStepsRegular_ = right.nTimeStepsRegular_;
    nTimeSteps_ = right.nTimeSteps_;
    TimeStepList_ = right.TimeStepList_;
    GsolnErrorNorm_ = right.GsolnErrorNorm_;
    GwdotErrorNorm_ = right.GwdotErrorNorm_;
    iCounter = right.iCounter;
    time_step_next = right.time_step_next;
    return *this;
}
//==================================================================================================================================
void SubIntegrationHistory::clear()
{
    nTimeStepsRegular_ = 0;
    nTimeSteps_ = 0;
    TimeStepList_.clear();
    GsolnErrorNorm_ = 0.0;
    GwdotErrorNorm_ = 0.0;
    iCounter = 0;
    time_step_next = 0.0;
}
//==================================================================================================================================
void  
SubIntegrationHistory::addTimeStep(double t_init, double t_final, double t_final_calc, int timeTypeSoln,
                                   int numNonLinSolves, double solnErrorNorm, double wdotErrorNorm, double volts,
				   double srcElectronsStep, double currentStep)
{
    /*
     * check for consistency
     */
    if (iCounter > 0) {
	TimeStepHistory& tshPrevious =  TimeStepList_[iCounter-1];
	double t_final_calc_Previous =  tshPrevious.t_final_calc_;
	if (fabs(t_init - t_final_calc_Previous) > 1.0E-13) {
	    throw ZZCantera::Electrode_Error("Shouldn't be here","");
	}
    }
    if (timeTypeSoln != 2) {
	if (fabs(t_final - t_final_calc) > 1.0E-13) {
	    throw Electrode_Error("Shouldn't be here", "");
	}
    } else {
	if (fabs(t_final - t_final_calc) < 1.0E-16) {
	    throw Electrode_Error("warning t_final == t_final_calc", "");
	}
    }
    // store time step (note inefficient op - may replace)
    TimeStepList_.push_back(
	TimeStepHistory(t_init, t_final, t_final_calc , timeTypeSoln, numNonLinSolves, solnErrorNorm, wdotErrorNorm)
	);
    TimeStepHistory &tsh = TimeStepList_[nTimeSteps_];
    tsh.volts_ = volts;
    tsh.srcTermStepElectrons_ = srcElectronsStep;
    tsh.currentStep_ = currentStep;
    iCounter++;

    time_step_next = t_final - t_init;
    GsolnErrorNorm_ += solnErrorNorm;
    GwdotErrorNorm_ += wdotErrorNorm;
    nTimeSteps_++;
    if (timeTypeSoln != 2) {
	nTimeStepsRegular_++;
    }
}
//==================================================================================================================================
void SubIntegrationHistory::zeroTimeStepCounter()
{
    iCounter = 0;
}
//==================================================================================================================================
void SubIntegrationHistory::advanceTimeStepCounter()
{
    if (iCounter + 1 < nTimeSteps_) {
	iCounter++;
    }
}
//==================================================================================================================================
double SubIntegrationHistory::getNextRegularTime(double currentTime) const
{
    // If we are beyond the stored time step history, then we used the stored value of time_step_next to make
    // up a good next time.
    if (iCounter >= nTimeSteps_) {
	const TimeStepHistory& tshCurrent =  TimeStepList_[iCounter];
	double tfinal = tshCurrent.t_final_;
	tfinal += time_step_next;
	return tfinal;
    }

    const TimeStepHistory* tshCurrent_ptr = & TimeStepList_[iCounter];
    double t_final_calc_First = tshCurrent_ptr->t_final_calc_;
    double t_final_First = tshCurrent_ptr->t_final_;
    double t_final = t_final_First;
    double t_final_calc = t_final_calc_First;

    // Step forward in the stored time steps until we find a t_final_calc beyond the input time, currentTime
    while (t_final_calc <= currentTime) {
	if (iCounter +1 < nTimeSteps_) {
	    iCounter++;
	    tshCurrent_ptr = & TimeStepList_[iCounter];
	    t_final_calc = tshCurrent_ptr->t_final_calc_;
	    t_final = tshCurrent_ptr->t_final_;
	} else {
	    t_final_calc = currentTime + time_step_next;
	    return t_final_calc;
	}
    }

    // step past special time steps that had discontinuities in them, until you find a regular step
    while (tshCurrent_ptr->timeTypeSoln_ == 2) {
	if (t_final_calc > t_final) {
	    return t_final;
	}
	if (iCounter +1 < nTimeSteps_) {
	    iCounter++;
	} else {
	    
	}
	tshCurrent_ptr = & TimeStepList_[iCounter];
	t_final = tshCurrent_ptr->t_final_;
	t_final_calc = tshCurrent_ptr->t_final_;
    }
 
    // Return a default next time based on the default deltaT if we have run out of times
    if (t_final_calc <= currentTime) {
	t_final_calc = currentTime + time_step_next;
    }
 
    // Normal return of the next stored time
    return t_final_calc;
}
//==================================================================================================================================
int SubIntegrationHistory::assureTimeInterval(double gtinit, double gtfinal)
{
    int iC = 0;
    TimeStepHistory* tshCurrent_ptr = & TimeStepList_[iC];
    double tinit =  tshCurrent_ptr->t_init_;
    if (fabs(gtinit - tinit) > 1.0E-13) {
	throw Electrode_Error("SubIntegrationHistory::assureTimeInterval", "The global time, " + fp2str(gtinit) +
                           ", doesn't agree with the storred initial time for the first step");
    }
    tshCurrent_ptr = &TimeStepList_[nTimeSteps_-1];
    double tfinal = tshCurrent_ptr->t_final_;
    if (tshCurrent_ptr->timeTypeSoln_ == 2) {
	tshCurrent_ptr->timeTypeSoln_ = 0;
    }
    if (fabs(gtfinal - tfinal) > 1.0E-13) {
	throw Electrode_Error("SubIntegrationHistory::assureTimeInterval", "The global final time, " + fp2str(gtfinal) +
                           ", doesn't agree with the storred final time for the last step");
    }
    int iCsave = iCounter;
    for (int i = 0; i < nTimeStepsRegular_; i++) {
	tfinal = getNextRegularTime(tinit);
	if (tfinal <= tinit) {
	   throw Electrode_Error("SubIntegrationHistory::assureTimeInterval", 
                              "final time is less than the initial time");
	}
	advanceTimeStepCounter();
	tinit = tfinal;
    }
    if (fabs(tfinal - gtfinal) > 1.0E-13) {
        throw Electrode_Error("SubIntegrationHistory::assureTimeInterval", "The global final time, " + fp2str(gtfinal) +
                           ", doesn't agree with the storred final time for the last step");
    }
    iCounter = iCsave;
    return iCounter;
}
//==================================================================================================================================
double SubIntegrationHistory::globalStartTime() const
{
    if  (TimeStepList_.size() > 0) {
	const TimeStepHistory& tshCurrent =  TimeStepList_[0];
	return tshCurrent.t_init_;
    }
    return 0.0;
}
//==================================================================================================================================
double SubIntegrationHistory::globalEndTime() const
{
    if  (nTimeSteps_ > 0) {
	const TimeStepHistory& tshCurrent =  TimeStepList_[nTimeSteps_-1];
	return tshCurrent.t_final_;
    }
    return 0.0;
}
//==================================================================================================================================
void SubIntegrationHistory::print(int lvl) const
{
    if (lvl > 1) {
	printf("     ===============================================================================================================================\n");
    }
    if (lvl > 0) {
	printf("     SubStepIntegrationHistory t_init = %12.5g to t_final = %12.5g nstepsReg = %d numSpecialSteps = %d\n",
	       globalStartTime(), globalEndTime(), nTimeStepsRegular_, nTimeSteps_- nTimeStepsRegular_);
    } 
    if (lvl == 2) {
	printf("            t_init     t_final      delta_t  solnType    t_requ numNonlin   solnErrDot wdotError\n");
	for (int i = 0; i < nTimeSteps_; i++) {
	    const TimeStepHistory& tshCurrent =  TimeStepList_[i];
	    double delta_t = tshCurrent.t_final_ - tshCurrent.t_init_;
	    printf("     %12.5g %12.5g %12.5g %5d %12.5g %5d %12.5g %12.5g\n", tshCurrent.t_init_, tshCurrent.t_final_calc_, 
		    delta_t, tshCurrent.timeTypeSoln_,  tshCurrent.t_final_calc_,
		   tshCurrent.numNonLinSolves_, tshCurrent.solnErrorNorm_, tshCurrent.wdotErrorNorm_);
	}
	printf("     ----------------------------------------------------------------------------------------------------------------------------\n");
	printf("                                                                     %12.5g %12.5g \n",    GsolnErrorNorm_, GwdotErrorNorm_);

	printf("     ==================================================================================================================================\n");
    }

    if (lvl > 2) {
	printf("            t_init     t_final      delta_t  solnType    t_requ numNonlin   solnErrDot wdotError "
	       " volts       srcElect       Current\n");
	double sSum = 0.0;
	double gtinit = globalStartTime();
	double gtfinal = globalEndTime();
	for (int i = 0; i < nTimeSteps_; i++) {
	    const TimeStepHistory& tshCurrent =  TimeStepList_[i];
	    double delta_t = tshCurrent.t_final_ - tshCurrent.t_init_;
	    sSum +=  tshCurrent.srcTermStepElectrons_;
	    printf("     %12.6g %12.6g %12.5g %5d %12.5g %5d %12.5g %12.5g %12.8f %12.5g %12.5g\n", 
		   tshCurrent.t_init_, tshCurrent.t_final_calc_, 
		   delta_t, tshCurrent.timeTypeSoln_,  tshCurrent.t_final_calc_,
		   tshCurrent.numNonLinSolves_, tshCurrent.solnErrorNorm_, tshCurrent.wdotErrorNorm_,
		   tshCurrent.volts_, tshCurrent.srcTermStepElectrons_, tshCurrent.currentStep_);
	}
	printf("     --------------------------------------------------------------------------------------------------------------------------------\n");
	double gCurr = 0.0;
	if (gtfinal > gtinit) {
	    gCurr = sSum / (gtfinal - gtinit) * Faraday;
	}
	printf("                                                                     %12.5g %12.5g              %12.5g %12.5g \n", 
	       GsolnErrorNorm_, GwdotErrorNorm_,
	       sSum, gCurr);

	printf("     =================================================================================================================================\n");
    }
}
//==================================================================================================================================
void SubIntegrationHistory:: setConstantStepSizeHistory(double gtinit, double gtfinal, int nsteps)
{
    zeroTimeStepCounter();
    double delta = (gtfinal - gtinit) / nsteps;
    for (int i = 0; i < nsteps; i++) {
	double tinit = gtinit + i * delta;
	double tfinal = gtinit + (i+1) * delta;
	addTimeStep(tinit, tfinal, tfinal, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
}
//==================================================================================================================================
bool SubIntegrationHistory::operator==(const SubIntegrationHistory& other) const
{
    if (nTimeStepsRegular_ != other.nTimeStepsRegular_) {
	return false;
    }
    int i_l = 0;
    int i_r = 0;
    do {
	const TimeStepHistory* tshCurrent =  &TimeStepList_[i_l];
	const TimeStepHistory* tshCurrent_other =  & other.TimeStepList_[i_r];
	if (! (*tshCurrent == *tshCurrent_other)) {
	    if (tshCurrent->timeTypeSoln_ == 2) {
		i_l++;
		break;
	    }
	    if (tshCurrent_other->timeTypeSoln_ == 2) {
		i_r++;
		break;
	    }
	    return false;
	}
	i_l++;
	i_r++;
    } while (i_l < nTimeSteps_ && (i_r < other.nTimeSteps_));
   
    return true;
}
//===================================================================================================================================
bool SubIntegrationHistory::operator!=(const SubIntegrationHistory& other) const
{
    return ! (*this == other);
}
//==================================================================================================================================
//==================================================================================================================================
//==================================================================================================================================
Electrode_Integrator::Electrode_Integrator() :
    Electrode(),
    ZZCantera::ResidJacEval(),
    deltaTsubcycleCalc_(0.0),
    rtolResidNLS_(1.0E-6),
    haveGood_solnDot_init_(false),
    predictDotBetter_(false),
    pSolveJAC_(nullptr),
    jacPtr_(nullptr),
    jacMng_(nullptr),
    numIntegratedSrc_(0),
    rtol_IntegratedSrc_global_(1.0E-4),
    maxNumberSubCycles_(50000),
    maxNumberSubGlobalTimeSteps_(1000),
    IntegratedSrc_normError_local_(0.0),
    IntegratedSrc_normError_global_(0.0),
    relativeLocalToGlobalTimeStepMinimum_(1.0E-3),
    extraRtolNonlinearSolver_(1.0)
{
}
//==================================================================================================================================
Electrode_Integrator::Electrode_Integrator(const Electrode_Integrator& right) :
    Electrode_Integrator()
{
    operator=(right);
}
//==================================================================================================================================
Electrode_Integrator& Electrode_Integrator::operator=(const Electrode_Integrator& right)
{
    if (this == &right) {
        return *this;
    }

    Electrode::operator=(right);
    ZZCantera::ResidJacEval::operator=(right);

    //neq_                                = right.neq_;
    deltaTsubcycleCalc_                 = right.deltaTsubcycleCalc_;

    atolResidNLS_                       = right.atolResidNLS_;
    rtolResidNLS_                       = right.rtolResidNLS_;
    atolNLS_                            = right.atolNLS_;
    rtolNLS_                            = right.rtolNLS_;
    ylowNLS_                            = right.ylowNLS_;
    yhighNLS_                           = right.yhighNLS_;
    yvalNLS_                            = right.yvalNLS_;
    yvalNLS_init_                       = right.yvalNLS_init_;
    yvalNLS_init_init_                  = right.yvalNLS_init_init_;
    yvalNLS_final_final_                = right.yvalNLS_final_final_;
    ydotNLS_                            = right.ydotNLS_;
    errorLocalNLS_                      = right.errorLocalNLS_;
    errorGlobalNLS_                     = right.errorGlobalNLS_;

    deltaBoundsMagnitudesNLS_           = right.deltaBoundsMagnitudesNLS_;

    soln_predict_                       = right.soln_predict_;
    haveGood_solnDot_init_              = right.haveGood_solnDot_init_;
    solnDot_init_                       = right.solnDot_init_;
    solnDot_final_                      = right.solnDot_final_;
    solnDot_final_final_                = right.solnDot_final_final_;
    solnDot_init_init_                  = right.solnDot_init_init_;
    soln_predict_fromDot_               = right.soln_predict_fromDot_;
    predictDotBetter_                   = right.predictDotBetter_;
   // pSolve_                             = new NonlinearSolver(this);
    SAFE_DELETE(pSolveJAC_);
    pSolveJAC_                          = new NonlinearSolver_JAC(this);

    SAFE_DELETE(jacPtr_);
    if (right.jacPtr_) {
        jacPtr_                             = new SquareMatrix(*right.jacPtr_);
    }
    SAFE_DELETE(jacMng_);
    if (right.jacMng_) {
        jacMng_                             = new JacobianManager(*right.jacMng_);
        jacMng_->reinstall_ptrs(this, jacPtr_, false);
    }

    numIntegratedSrc_                   = right.numIntegratedSrc_;
    IntegratedSrc_Predicted             = right.IntegratedSrc_Predicted;
    IntegratedSrc_final_                = right.IntegratedSrc_final_;
    IntegratedSrc_Errors_local_         = right.IntegratedSrc_Errors_local_;
    IntegratedSrc_Errors_globalStep_    = right.IntegratedSrc_Errors_globalStep_;
    rtol_IntegratedSrc_global_          = right.rtol_IntegratedSrc_global_;
    atol_IntegratedSrc_global_          = right.atol_IntegratedSrc_global_;
    maxNumberSubCycles_                 = right.maxNumberSubCycles_;
    maxNumberSubGlobalTimeSteps_        = right.maxNumberSubGlobalTimeSteps_;
    IntegratedSrc_normError_local_      = right.IntegratedSrc_normError_local_;
    IntegratedSrc_normError_global_     = right.IntegratedSrc_normError_global_;

    timeHistory_base_                   = right.timeHistory_base_;
    timeHistory_current_                = right.timeHistory_current_;

    relativeLocalToGlobalTimeStepMinimum_ = right.relativeLocalToGlobalTimeStepMinimum_;
    extraRtolNonlinearSolver_           = right.extraRtolNonlinearSolver_;
    return *this;
}
//==================================================================================================================================
Electrode_Integrator::~Electrode_Integrator()
{
    SAFE_DELETE(jacPtr_);
    SAFE_DELETE(jacMng_);
    SAFE_DELETE(pSolveJAC_);
}
//==================================================================================================================================
Electrode* Electrode_Integrator::duplMyselfAsElectrode() const 
{
    Electrode_Integrator* dd = new Electrode_Integrator(*this);
    return dd;
}
//==================================================================================================================================
int Electrode_Integrator::electrode_model_create(ELECTRODE_KEY_INPUT* ei)
{
    Electrode::electrode_model_create(ei);

    setupIntegratedSourceTermErrorControl();
    //
    // Gather the input parameters for time stepping that were located in the Electrode input deck
    //
    relativeLocalToGlobalTimeStepMinimum_ = ei->relativeLocalToGlobalTimeStepMinimum;
    extraRtolNonlinearSolver_             = ei->extraRtolNonlinearSolver;

    maxNumberSubGlobalTimeSteps_ = ei->maxNumberSubGlobalTimeSteps;

    return 0;
}
//==================================================================================================================================
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
int  Electrode_Integrator::setInitialConditions(ELECTRODE_KEY_INPUT* ei)
{
    Electrode::setInitialConditions(ei);
    nEquations_calc();
    create_solvers();
    return 0;
}
//==================================================================================================================================
size_t Electrode_Integrator::nEquations_calc() const
{
    throw Electrode_Error("Electrode_Integrator::nEquations_calc()", "Parent member function called");
    return 0;
}
//==================================================================================================================================
int Electrode_Integrator::create_solvers()
{
    //
    //  Note: We must have a valid count of the number of unknowns in the nonlinear problem
    //        when we call this. If we don't, then the solver arrays won't be malloced with the right amount of space.
    size_t neqNLS = nEquations();
    neq_ = neqNLS;

    if (pSolveJAC_) {
        int m = nEquations();
        if (m != (int) yvalNLS_.size()) {
            printf("shouldn't be here\n");
            exit(-1);
        }
        return nEquations();
    }

    yvalNLS_.resize(neqNLS, 0.0);
    yvalNLS_final_final_.resize(neqNLS, 0.0);
    yvalNLS_init_init_.resize(neqNLS, 0.0);
    yvalNLS_init_.resize(neqNLS, 0.0);

    ydotNLS_.resize(neqNLS, 0.0);
    ylowNLS_.resize(neqNLS, 0.0);
    yhighNLS_.resize(neqNLS, 0.0);
    errorLocalNLS_.resize(neqNLS, 0.0);
    errorGlobalNLS_.resize(neqNLS, 0.0);

    // Set all absolute tolerances to a default of 1.0E-12. We will do better when we know what each unknown is
    atolNLS_.resize(neqNLS, 1.0E-12);
    // Set all absolute tolerances to a default of 1.0E-12. We will do better when we know what each unknown is
    atolResidNLS_.resize(neqNLS, 1.0E-12);
    deltaBoundsMagnitudesNLS_.resize(neqNLS, 1.0E300);

    // Add a couple of extra doubles, to the predictor, because some objects store extra info in those slots
    soln_predict_.resize(neqNLS + 2, 0.0);
    soln_predict_fromDot_.resize(neqNLS + 2, 0.0);

    solnDot_init_.resize(neqNLS, 0.0);
    solnDot_final_.resize(neqNLS, 0.0);
    solnDot_init_init_.resize(neqNLS, 0.0);
    solnDot_final_final_.resize(neqNLS, 0.0);

    IntegratedSrc_Predicted.resize(numIntegratedSrc_, 0.0);
    IntegratedSrc_final_.resize(numIntegratedSrc_, 0.0);
    IntegratedSrc_Errors_local_.resize(numIntegratedSrc_, 0.0);
    IntegratedSrc_Errors_globalStep_.resize(numIntegratedSrc_, 0.0);
    atol_IntegratedSrc_global_.resize(numIntegratedSrc_, 0.0);


    jacPtr_ = new SquareMatrix(neqNLS, 0.0);

    jacMng_ = new JacobianManager(this, jacPtr_);
    pSolveJAC_ = new NonlinearSolver_JAC(this);

    return neqNLS;
}
//==================================================================================================================================
int Electrode_Integrator::setupIntegratedSourceTermErrorControl()
{
    int numDofs = 0;
    if (metalPhase_ >= 0) {
        numDofs++;
    }
    if (solnPhase_ >= 0) {
        ThermoPhase& tp = thermo(solnPhase_);
        int nsp = tp.nSpecies();
        numDofs += nsp;
    }
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
    double sm = SolidTotalMoles();
    for (size_t i = 0; i <  numIntegratedSrc_; i++) {
        atol_IntegratedSrc_global_[i] = 1.0E-14 * sm;
    }
    return numDofs;
}
//==================================================================================================================================
bool Electrode_Integrator::resetStartingCondition(double Tinitial, bool doAdvancementAlways)
{
    bool resetToInitInit = false;
    /*
     * If the initial time is input, then the code doesn't advance. We clear stuff and redo the current global time step
     */
    double tbase = std::max(t_init_init_, 1.0E-50);
    if (fabs(Tinitial - t_init_init_) < (1.0E-13 * tbase) && !doAdvancementAlways) {
        resetToInitInit = true;
    }

    /*
     *  Clear the time step histories for the base and current timeHistories.
     */
    timeHistory_base_.clear();
    timeHistory_current_.clear();
    //  Call the base class function
    bool rr = Electrode::resetStartingCondition(Tinitial, doAdvancementAlways);
    if (rr != resetToInitInit) {
        throw Electrode_Error("Electrode_Integrator::resetStartingCondition()",
                              "Inconsistent resetToInitInit value ");
    }
    
    /*
     *  Zero the global error vectors
     */
    std::fill(errorGlobalNLS_.begin(), errorGlobalNLS_.end(), 0.0);
    int neqNLS = nEquations();
    /*
     *  Advance the solution to the state of the final_final values from the previous global step
     */
    if (!resetToInitInit) {
        for (int i = 0; i < neqNLS; ++i) {
            solnDot_init_init_[i] = solnDot_final_final_[i];
            solnDot_init_[i]      = solnDot_final_final_[i];
            solnDot_final_[i]     = solnDot_final_final_[i];
            yvalNLS_init_init_[i] = yvalNLS_final_final_[i];
            yvalNLS_init_[i]      = yvalNLS_final_final_[i];
            yvalNLS_[i]           = yvalNLS_final_final_[i];
        }
    }
    return resetToInitInit;
}
//==================================================================================================================================
//  Calculate the change in the state of the system when integrating from T_initial_initial_
//  to t_final_final_
/*
 *  All information is kept internal within this routine. This may be done continuously
 *  and the solution is not updated.
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
 *          Failure due to maxsubcycles exceeded returns a -4. 
 *
 *  When the max subcycles are exceeded, can check the subIntegrationHistory for the time that the
 *  cycle did get to. Also, the value of t_final_final_ is set to the final time that the integration got to.
 *
 */
int  Electrode_Integrator::integrate(double deltaT, double  GlobalRtolSrcTerm,
                                     Electrode_Exterior_Field_Interpolation_Scheme_Enum fieldInterpolationType,
                                     Subgrid_Integration_RunType_Enum subIntegrationType)
{
    double tfinal_start;
    /*
     *  Need to turn off following electrolyte moles. There is no solution if this is not true.
     */
    turnOffFollowElectrolyteMoles();
    /*
     *  Clear the time history for the current step
     */
    timeHistory_current_.clear();
   

    double pnorm = 0.0;
    int num_newt_its = 0; // Number of newton iterations in the nonlinear solve


    // Create the solvers and other vectors if we haven't already
    //  -> Note we haven't malloced the memory initially because we need to know the
    //     number of equations.
    if (pSolveJAC_ == 0) {
        size_t neqNLS = nEquations_calc();
        create_solvers();
        if (neqNLS != neq_) {
            exit(-1);
        }
    }

    /*
     *  Guard against a zero deltaT
     */
    if (deltaT <= 0.0) {
        return -3;
    }

    /*
     *  When we call this routine successfully we have an integration for the current step pending
     *  Tempory data now exist for t_final_
     */
    pendingIntegratedStep_ = 1;

    /*
     *  Set the local counter for the number of subcycles
     */
    int iterSubCycle = 0;

    /*
     *  Increment the global counter for the number of calls to this routine
     */
    counterNumberIntegrations_++;

    /*
     *  Set the counter for the number of consequative failures of the time stepping method that happens for any reason
     */
    int conseqFailures = 0;

    // Set the flag for overall time stepping to the goal of t_final_final_
    bool notDone = true;

    /*
     *  Set the inital time steps and final_final time step based on the overall step size, deltaT
     */
    tinit_ = t_init_init_;
    tfinal_ = tinit_;
    t_final_final_ = t_init_init_ + deltaT;

    /*
     * Choose the starting value of the subcycle time step
     *  We won't enforce a beginning min time step due to relativeLocalToGlobalTimeStepMinimum or max number of time
     *  steps.  We may need to start really small for some test problems. And the deltaT coming in will
     *  reflect these trade offs.
     */
    deltaTsubcycleNext_ = MIN(deltaT, deltaTsubcycleMax_);
    deltaTsubcycleNext_ = MIN(deltaTsubcycleNext_, deltaTsubcycle_init_init_);
    if (choiceDeltaTsubcycle_init_ == 2) {
        deltaTsubcycleNext_ = deltaT / 10.;
    }

    //     Put a max of two to the growth of deltaT for the next global time step
    if (choiceDeltaTsubcycle_init_ == 0) {
        deltaTsubcycle_init_next_ = deltaTsubcycleNext_;
    } else {
        deltaTsubcycle_init_next_ = 2.0*deltaT;
    }

    if (subIntegrationType == FVDELTA_TIMEINTEGRATION_SIR || subIntegrationType == FIXEDSUBCYCLENUMBER_TIMEINTEGRATION_SIR) {
	// Set counter to zero
	timeHistory_base_.zeroTimeStepCounter();
	timeHistory_base_.assureTimeInterval(t_init_init_, t_final_final_);
	double tt = timeHistory_base_.getNextRegularTime(tinit_);
	deltaTsubcycleNext_ = tt - tinit_;
    }

    /*
     *   We calculate the number of estimated steps to take in the subcycling. This will be used to adjust the
     *   tolerance requirements needed for each of the steps
     */
    int  nsteps_est = deltaT * 0.99 / deltaTsubcycleNext_ + 1;

    /*
     *   Set the init and final state from the init_init state
     */
    setInitStateFromInitInit(true);

    /*
     *   Figure out what phases exist and tell the kinetics operator
     */
    setPhaseExistenceForReactingSurfaces(false);

    /*
     *   Set up the global tolerances for the source term vector. Here we firm up what the source term vector is
     */
    setNLSGlobalSrcTermTolerances(GlobalRtolSrcTerm);

    /*
     *   Set the relative tolerance on the nonlinear solver
     */
    rtolNLS_      = GlobalRtolSrcTerm;
    rtolResidNLS_ = GlobalRtolSrcTerm * 1.0E-2;

    //  Set the absolute tolerance vector for the nonlinear solver, both the residual and the solution tolerances.
    //    -  atolNLS_[]
    //    -  atolResidNLS_[]
    //    We control both because we have had trouble making sure that the equations are solved to a sufficient accuracy.
    setResidAtolNLS();
    
    //  Zero vectors that are accumulated over local time step to represent global time step quantities,
    //       the integrated source term. This will be accumulated over the subcycling estimated global time step errors
    zeroGlobalStepAccumulationTerms();

    //  Update the state
    updateState();
    
    //  Determine the big mole fractions of multispecies phases
    determineBigMoleFractions();

    //  Save the Electrode state into an XML state object
    if (eState_save_) {
#ifdef DEBUG_CHECK_XML
        if (xmlStateData_final_) {
            bool retn = check_XML_valid(xmlStateData_final_);
            if (!retn) {
                throw Electrode_Error("Electrode_Integrate() 1", "xmlStateData_final_ corrupted");
            }
        }
        if (xmlStateData_init_) {
            bool retn = check_XML_valid(xmlStateData_init_);
            if (!retn) {
                throw Electrode_Error("Electrode_Integrate() 1", "xmlStateData_init_ corrupted");
            }
        }
#endif
        eState_save_->copyElectrode_intoState(this, false);
        SAFE_DELETE(xmlStateData_init_);
        xmlStateData_init_ = eState_save_->write_electrodeState_ToXML();
        startXML_TI_final();
#ifdef DEBUG_CHECK_XML
        if (xmlStateData_final_) {
            bool retn = check_XML_valid(xmlStateData_final_);
            if (!retn) {
                throw Electrode_Error("Electrode_Integrate() 2", "xmlStateData_final_ corrupted");
            }
        }
        if (xmlStateData_init_) {
            bool retn = check_XML_valid(xmlStateData_init_);
            if (!retn) {
                throw Electrode_Error("Electrode_Integrate() 2", "xmlStateData_init_ corrupted");
            }
        }
#endif
    }
    
    // *****************************************************************************************************************
    //  Loop over any possible subcycling of the time step
    do {
        /*
         *  Increment the local subcycle counter
         */
        iterSubCycle++;
        /*
         *  Increment the grand counter, which can be used to isolate one bad integration amongst the thousands
         */
        counterNumberSubIntegrations_++;

        /*
         *  Set the initial time interval
         */
        tinit_ = tfinal_;

        /*
         *   Update initial values of the state vector with the final values from the last subcycle
         *   iteration, if we are beyond the first subcycle
         */
        if (iterSubCycle > 1) {
            setInitStateFromFinal(false);
        }

	/*
	 * Set up the yvalNLS_init 
	 */
	check_yvalNLS_init(false);

        /*
         *  Calculate the time interval from the previous step or from the initial value coming into the routine
         */
        deltaTsubcycle_ = deltaTsubcycleNext_;
#ifdef DEBUG_MODE
	if (s_printLvl_DEBUG_SPECIAL && deltaTsubcycle_ < 0.04) {
	    // printf("WARNING deltaTubcycle_ = %g\n", deltaTsubcycle_ );
	}
#endif
	/*
	 *   Or if we have specified to follow a time step history wo get the next step size from the history counter
	 */
        //	if (subIntegrationType == FVDELTA_TIMEINTEGRATION_SIR || subIntegrationType == FIXEDSUBCYCLENUMBER_TIMEINTEGRATION_SIR) {
	//   tfinal_ = timeHistory_base_.getNextRegularTime();
	// deltaTsubcycle_ = tfinal_ - tinit_;
	//}

#ifdef DEBUG_MODE_NOT
        if (counterNumberIntegrations_ >= 27502) {
            if (electrodeCellNumber_ == -1 && electrodeDomainNumber_ == 2) {
                //enableExtraPrinting_ =15;
                //detailedResidPrintFlag_ =15;
                //printLvl_ =10;
                printf("we are here c6 cNI = %d cNSI = %d\n", counterNumberIntegrations_, counterNumberSubIntegrations_);
            }
        }

        if (counterNumberSubIntegrations_ >= 134632) {
            if (electrodeCellNumber_ == -1) {
                //enableExtraPrinting_ =15;
                //detailedResidPrintFlag_ =15;
                //printLvl_ =10;
                printf("we are here 3\n");
            }
        }
#endif

        /*
         * Set the final time of the interval, but bound the time interval by the global time step
         */
        tfinal_ = tinit_ + deltaTsubcycle_;
        if (tfinal_ > t_init_init_ + deltaT) {
            tfinal_ = t_init_init_ + deltaT;
            deltaTsubcycle_ = tfinal_ - tinit_;
            nsteps_est =  iterSubCycle;
        }

        /*
         * Set the calculated deltaTsubcycle to the current deltaTsubcycle
         */
        deltaTsubcycleCalc_ = deltaTsubcycle_;
        /*
         *  Check internal consistency of the init solution
         *    -> failures here produce error exits
         *    -> check
         */
        check_init_consistency();

        /*
         *   Prepare the problem for predicting the solution
         *      -> determine the largest mole fraction in a phase
         *      -> determine if a surface phase can have reactions turned on
         *      -> do stuff that affects the equation system that will be solved at the
         *         current step.
         */
        prepareProblemStatement();

        /*
         *  Set a flag indicating whether the nonlinear problem has converged or not
         */
        int nonlinConverged = 0;

        /*
         *  Set a counter to determine how many times the subproblem has been attempted
         */
        int attemptSubCycle = 0;

        /*
         *  Loop until we have convergence of an acceptable step
         */
        while (nonlinConverged <= 0) {
topConvergence:
            /*
             *    Increment the subcycle attempts
             */
            attemptSubCycle++;
            if (attemptSubCycle > 10) {
                throw Electrode_Error("Electrode_Integrator::integrate()", "FAILURE too many nonlinear convergence failures");
            }
	    tfinal_ = tinit_ + deltaTsubcycle_;
	    tfinal_start = tfinal_;
            /*
             *   Zero needed counters
             */
            std::fill(justDiedPhase_.begin(), justDiedPhase_.end(), 0);
            std::fill(justBornPhase_.begin(), justBornPhase_.end(), 0);
            /*
             *  Ok at this point we have a time step deltaTsubcycle_ and initial conditions consisting of 
             *  phaseMoles_init_ and spMF_init_. We now calculate predicted solution components from these conditions.
             *  Additionally, we evaluate whether any multispecies phases are going to pop into existence,
             *  adding in potential seed values and we set the phaseExistence flags in the kinetics solver
             */
            int info = predictSoln();
#ifdef DEBUG_MODE
	 //   if (s_printLvl_DEBUG_SPECIAL && deltaTsubcycle_ < 0.04) {
	 //       printf("WARNING deltaTubcycle_ = %g\n", deltaTsubcycle_ );
	 //   }
#endif
            if (info != 1) {
                conseqFailures++;
                nonlinConverged = 0;
                if (printCSVLvl_) {
                    writeCSVData(-2);
                }
                tfinal_ = tinit_;
                deltaTsubcycle_ = deltaTsubcycle_ * 0.25;
                deltaTsubcycleNext_ = deltaTsubcycle_;
                deltaTsubcycleCalc_ = deltaTsubcycle_;
                tfinal_ += deltaTsubcycle_;
                /*
                 * Revert to the old solution -> copy _init_ to _final_
                 */
                setFinalStateFromInit();
                goto topConvergence;
            }

	    /*
             *  We've done an explicit predictor. However, for some systems this isn't a good predictor, or it's too expenseive
	     *  and time consuming to do a good job. Therefore, we will back it up with a traditional predictor -corrector method.
	     *  For some other systems, this step can be bypassed. Whatever is more accurate is used as the error indicator
             *
             *  If the predictDotBetter_ variable is true, then we actually use the predicted solution from predictSolnDot as 
             *  the initial guess to the nonlinear equation solver. If that prediction produces an out of bounds result, we 
             *  fail the solution step. and try again with a smaller time step after setting predictDotBetter_ to false;
             */
	    info = predictSolnDot();
	    if (info != 1) {
                conseqFailures++ ;
                nonlinConverged = 0;
                if (printCSVLvl_) {
                    writeCSVData(-2);
                }
                tfinal_ = tinit_;
                deltaTsubcycle_ = deltaTsubcycle_ * 0.25;
                deltaTsubcycleNext_ = deltaTsubcycle_;
                deltaTsubcycleCalc_ = deltaTsubcycle_;
                tfinal_ += deltaTsubcycle_;
                /*
                 * Revert to the old solution -> copy _init_ to _final_
                 */
                setFinalStateFromInit();
                goto topConvergence;
            }


            /*
             *  Gather the predicted integrated src prediction and the solution prediction
             */
            gatherIntegratedSrcPrediction();

            nsteps_est =  iterSubCycle + (t_init_init_ + deltaT - tinit_) / deltaTsubcycle_;
            if (nsteps_est > 100) {
                nsteps_est = 100;
            }
            if (nsteps_est < 1) {
                nsteps_est = 1;
            }

            //	rtolNLS_ = GlobalRtolSrcTerm/ nsteps_est;
            rtolNLS_ = GlobalRtolSrcTerm;

            /*
             *  We have a successful prediction: go and prep the nonlinear problem
             */
            num_newt_its = 0;
            int num_linear_solves = 0;
            int numBacktracks = 0;
            int solverPrintLvl = 0;
            if (s_printLvl_SOLVER) {
                solverPrintLvl = MAX(0, printLvl_ - 2);
                solverPrintLvl = MAX(solverPrintLvl, s_printLvl_SOLVER - 1);
            }

            /*
             *  formulate the nonlinear solver problem to be solved.
             *     Fields to be filled in
             *             yvalNLS_
             *             ylowNLS_
             *             yhighNLS_
             *             atolNLS_
             *             deltaBoundsMagnitudesNLS_
             */
            initialPackSolver_nonlinFunction();
            /*
             *  Calculate consistent ydotNLS_ going into the nonlinear solver
             */
            calc_ydotNLS_final();

            jacPtr_->m_useReturnErrorCode = 1;

            /*
             *  Check the nonlinear residual equations for completeness and the ability to be solved
             */
            int retn = check_nonlinResidConditions();
#ifdef DEBUG_MODE
	    if (s_printLvl_DEBUG_SPECIAL && deltaTsubcycle_ < 0.04) {
		//	printf("WARNING deltaTubcycle_ = %g\n", deltaTsubcycle_ );
	    }
#endif
            if (retn < 0) {
                conseqFailures++ ;
                nonlinConverged = 0;
                writeCSVData(-2);
                tfinal_ = tinit_;
                deltaTsubcycle_ = deltaTsubcycle_ * 0.25;
                deltaTsubcycleNext_ = deltaTsubcycle_;
                deltaTsubcycleCalc_ = deltaTsubcycle_;
                tfinal_ += deltaTsubcycle_;
                setFinalStateFromInit();
                goto topConvergence;
            }


            int solnType = NSOLN_TYPE_STEADY_STATE;
            //int solnType = NSOLN_TYPE_TIME_DEPENDENT;

            /*
             *  Set the tolerances on the solution variables
             */
           pSolveJAC_->setAtol(DATA_PTR(atolNLS_));
           pSolveJAC_->setRtol(0.1 *  extraRtolNonlinearSolver_ * rtolNLS_);
           // } else {
           //     pSolve_->setAtol(DATA_PTR(atolNLS_));
           //     pSolve_->setRtol(0.1 * extraRtolNonlinearSolver_ * rtolNLS_);
           // }
            /*
             *  Set the tolerances on the residual evaluation of the nonlinear solver
             *    We set the relative tolerance of the residual to GlobalRtolElectronSrcTerm/ nsteps_est. This is the
             *    same value as the relative tolerance on the electrode state variables.
             *
             *    Set the absolute tolerances for the residual calculation to a vector residAtolNLS_[] which is calculated
             *    for the problem and is associated with the degrees of freedom in the problem. residAtolNLS_[] has
             *    units of the residual.
             *
             *    The residual tolerance is given by the minimum of the row weight sums and the user tolerances given below.
             */
            double extraR = 1.0E-2;
           pSolveJAC_->setResidualTols(extraR * extraRtolNonlinearSolver_ * rtolResidNLS_,  &atolResidNLS_[0]);
           pSolveJAC_->setBoundsConstraints(&ylowNLS_[0], &yhighNLS_[0]);
           pSolveJAC_->setDeltaBoundsMagnitudes(DATA_PTR(deltaBoundsMagnitudesNLS_));
           // } else {
           //     pSolve_->setResidualTols(extraR * extraRtolNonlinearSolver_ * rtolResidNLS_,  &atolResidNLS_[0]);
           //     pSolve_->setBoundsConstraints(&ylowNLS_[0], &yhighNLS_[0]);
           //     pSolve_->setDeltaBoundsMagnitudes(DATA_PTR(deltaBoundsMagnitudesNLS_));
           // }
	 
            num_newt_its = 0;

#ifdef DEBUG_HKM_1DELECTRODE
            if (electrodeCellNumber_ == 0) {
                if (counterNumberIntegrations_ == 40) {
                    solverPrintLvl = 10;
                    detailedResidPrintFlag_ = 10;
                    enableExtraPrinting_ = 1;
                } else {
                    solverPrintLvl = 0;
                    detailedResidPrintFlag_ = 0;
                    enableExtraPrinting_ = 0;
                }
            }
#endif
            //   How to turn on the jacobian printing if debugging
            //pSolve_->s_print_NumJac = 1;
            //solverPrintLvl = 10;
	    //pSolve_->m_min_newt_its = 4;
#ifdef DEBUG_MODE
	    if (s_printLvl_DEBUG_SPECIAL && deltaTsubcycle_ < 0.04) {
		//	printf("WARNING deltaTubcycle_ = %g  iterSubCycle = %d countSub = %d\n", 
		//     deltaTsubcycle_ , iterSubCycle, counterNumberIntegrations_ );
	    }
#endif
            int nonlinearFlag;
            nonlinearFlag = pSolveJAC_->solve_nonlinear_problem(solnType, &yvalNLS_[0], &ydotNLS_[0], 0.0,
			                                        tfinal_, *jacMng_,  num_newt_its, num_linear_solves,
						                numBacktracks, solverPrintLvl);
            if (nonlinearFlag < 0) {
                if (printLvl_ > 2) {
                    printf("Electrode_Integrate::integrate(): Unsuccessful Nonlinear Solve flag = %d\n", nonlinearFlag);
                }
            }
	    //bool specialSolve = false;
otherFailureType:
            if (nonlinearFlag < 0) {
                conseqFailures++ ;
                nonlinConverged = 0;
                if (printCSVLvl_) {
                    writeCSVData(-1);
                }
                tfinal_ = tinit_;
                deltaTsubcycle_ = deltaTsubcycle_ * 0.25;
                deltaTsubcycleNext_ = deltaTsubcycle_;
                tfinal_ += deltaTsubcycle_;
                
                // Revert to the old solution -> copy _init_ to _final_
                setFinalStateFromInit();
            } else {
                //  HKM -> Don't think this can happen, but checking anyway
                info = unpackNonlinSolnVector(DATA_PTR(yvalNLS_));
                if (info != 1) {
                    if (printLvl_ > 1) {
                        printf("Electrode_Integrate::integrate(): WARNING "
                               "Nonlinear problem produced a solution that had a bounds problem!\n");
                    }
                   nonlinearFlag = -1;
                   goto otherFailureType;
                }

                //  Correct the final time
                if (deltaTsubcycleCalc_ < deltaTsubcycle_ * (1.0 - 1.0E-10)) {
                    if (printLvl_ > 1) {
                        printf("deltaT reduced to %g from %g due to phase death capture\n",
                               deltaTsubcycleCalc_,  deltaTsubcycle_);
                    }
                    // Pick the next delta T to be equal to the current delta T
                    deltaTsubcycleNext_ = deltaTsubcycle_;
                    tfinal_ = tinit_ + deltaTsubcycleCalc_;
		    //specialSolve = true;
                } else  if (deltaTsubcycleCalc_ > deltaTsubcycle_ * (1.0 + 1.0E-10)) {
                    if (printLvl_ > 1) {
                        printf("deltaT increased to %g from %g due to phase death capture\n",
                               deltaTsubcycleCalc_,  deltaTsubcycle_);
                    }
                    // Pick the next delta T to be equal to the current delta T
                    deltaTsubcycleNext_ = deltaTsubcycle_;
                    tfinal_ = tinit_ + deltaTsubcycleCalc_;
		    //specialSolve = true;
                }

                conseqFailures--;
                conseqFailures = MAX(0, conseqFailures);
                nonlinConverged = 1;
            }
        }  // End of convergence of the nonlinear iteration for the current subcycle -> we pass here if nonlinConverged=1

        //
        //  ACTIVITIES THAT ARE NEEDED TO DETERMINE IF THE SUBCYCLE IS SUCCESSFUL EVEN IF NONLINEAR SOLVER CONVERGES
        //
        bool stepAcceptable = true;
        stepAcceptable = checkSubIntegrationStepAcceptable();
        /*
         *  Massage the results due to birth and death, before gauging the accuracy of the time stepping.
         */
        bool sa = changeSolnForBirthDeaths();
        if (!sa) {
            stepAcceptable = false;
        }

        /*
         * Print out intermediate results if the print level is high enough
         */
       if (printLvl_ > 3) {
            calcSrcTermsOnCompletedStep();
            accumulateSrcTermsOnCompletedStep();
            printElectrode(true, true);
            accumulateSrcTermsOnCompletedStep(true);
        }
#ifdef DEBUG_MODE
       if (s_printLvl_DEBUG_SPECIAL && iterSubCycle > 10) {
	   if (iterSubCycle >= 21) {
	       printf("iterSubCycle = %d\n", iterSubCycle);
	   } else {
	       printf("iterSubCycle = %d\n", iterSubCycle);
	   }
       }
#endif
        //
        // ACTIVITIES THAT ARE NEEDED TO DETERMINE IF THE SUBCYCLE IS ACCURATE ENOUGH
        // 
	double pnormSrc = 0.0;
	double pnormSoln = 0.0;
        if (stepAcceptable) {
	    /*
	     * Calculate solnDot_final_ so that we can use it to evaluate the next step as a predictor.
	     */
	    calc_solnDot_final();
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
            pnormSrc = predictorCorrectorGlobalSrcTermErrorNorm();

            pnormSoln = predictorCorrectorWeightedSolnNorm(yvalNLS_);

            pnorm = fmax(pnormSrc, pnormSoln);

            /*
             * Adjust the next time step according to a predictor-corrector critera
             */

            /*
             * If we are way over the allowed local subcycle tolerance requirement, reject the step
             */
            if (pnorm > 2.0) {
                stepAcceptable = decide_normHighLogic(pnorm, num_newt_its, iterSubCycle);
            }

            /*
             * Print out a table of predictions versus corrections. This is a good way to debug the
             * predictor -corrector process
             */
            if (s_printLvl_PREDICTOR_CORRECTOR) {
#ifdef DEBUG_MODE_PREDICTION
                predictorCorrectorPrint(yvalNLS_, pnormSrc, pnormSoln);
#else
                if (printLvl_ > 4) {
                    predictorCorrectorPrint(yvalNLS_, pnormSrc, pnormSoln);
                }
#endif
            }
        }
        /*
         *  HANDLE THE CASE WHERE WE HAVE NOT FOUND THE CURRENT STEP TO BE ACCEPTABLE
         *   - we decrease the time step and then redo the calculation
         */
        if (! stepAcceptable) {
#ifdef DEBUG_MODE_NOT
            printLvl_ = 9;
            enableExtraPrinting_ = 9;
            detailedResidPrintFlag_ = 9;
#endif
            conseqFailures++ ;
            nonlinConverged = 0;
            if (printCSVLvl_) {
                writeCSVData(-1);
            }
            tfinal_ = tinit_;
            deltaTsubcycle_ = deltaTsubcycle_ * 0.25;
	    if (deltaTsubcycle_ < relativeLocalToGlobalTimeStepMinimum_ * deltaT) {
		deltaTsubcycle_ = relativeLocalToGlobalTimeStepMinimum_ * deltaT;
	    }
            deltaTsubcycleNext_ = deltaTsubcycle_;
            tfinal_ += deltaTsubcycle_;
            /*
             * Revert to the old solution -> copy _init_ to _final_
             */
            setFinalStateFromInit();

            goto topConvergence;
        }

        //
        //  ACTIVITIES THAT ARE NEEDED AT THE END OF A SUCCESSFUL SUBCYCLE
        //

        /*
         *  Massage the results now that we have a successful step. This will mean zeroing out small phases
         */
        manageBirthDeathSuccessfulStep();

        /*
         * Accumulate local results into global vectors
         */
        accumulateSrcTermsOnCompletedStep();

        /*
         * Accumulate local error values into global error vectors on a completed Step
         */
        accumulateLocalErrorToGlobalErrorOnCompletedStep();

	int timeTypeSoln = 0;
	if (fabs( tfinal_start - tfinal_) > 1.0E-13) {
	    timeTypeSoln = 2;
	}
	 
	/*
	 *   Find the local current for this step using a virtual function so that it is extensible. Then find the
	 *   src term using this value of the current
	 *       units = current = coulombs/sec = amps
	 *               srcTerm = kmol / sec
	 */
	double currentStep = integratedLocalCurrent();
	double srcElectronsStep = currentStep * (tfinal_ - tinit_) / Faraday;
        timeHistory_current_.addTimeStep(tinit_, tfinal_start, tfinal_, timeTypeSoln,
					 num_newt_its, pnormSoln, pnormSrc, deltaVoltage_, srcElectronsStep, currentStep);

        /*
         * Adjust the next time step according to a predictor-corrector critera
         */
        if (choiceDeltaTsubcycle_init_ != 2) {
            if (pnorm < 0.05) {
                if (num_newt_its < 14 &&  conseqFailures == 0) {
                    if (deltaTsubcycleNext_ <= deltaTsubcycle_) {
                        // unusual case here
                        if ((1.000000001 * tfinal_) < t_init_init_ + deltaT) {
                            deltaTsubcycleNext_ = 1.5 * deltaTsubcycle_;
                        }
                    }
                }
            } else {
                if (pnorm < 0.10 && num_newt_its < 8 && conseqFailures == 0) {
                    if ((1.000000001 * tfinal_) < t_init_init_ + deltaT) {
                        deltaTsubcycleNext_ = 1.5 * deltaTsubcycle_;
                    }
                }
                if (deltaTsubcycle_ > relativeLocalToGlobalTimeStepMinimum_ * deltaT) {
                    if (pnorm > 0.10) {
                        if (deltaTsubcycleNext_ > deltaTsubcycle_) {
                            deltaTsubcycleNext_ = deltaTsubcycle_;
                        }
                    }
                    if (pnorm > 0.30) {
                        if (deltaTsubcycleNext_ > 0.5 * deltaTsubcycle_) {
                            deltaTsubcycleNext_ = 0.5 * deltaTsubcycle_;
                            if (deltaTsubcycleNext_ < relativeLocalToGlobalTimeStepMinimum_ * deltaT) {
                                deltaTsubcycleNext_ = 1.2 * deltaTsubcycle_;
                                deltaTsubcycleNext_ = relativeLocalToGlobalTimeStepMinimum_ * deltaT;
                            }
                        }
                    }
                }
            }
        }

        /*
         *  Enforce the minimum time step requirement
         */
        if (deltaTsubcycleNext_ < relativeLocalToGlobalTimeStepMinimum_ * deltaT) {
            deltaTsubcycleNext_ = 1.2 * deltaTsubcycle_;
            if (pnorm < 0.2) {
                deltaTsubcycleNext_ = 1.5 * deltaTsubcycle_;
            }
            if (pnorm < 0.05) {
                deltaTsubcycleNext_ = 2.0 * deltaTsubcycle_;
            }
            if (deltaTsubcycleNext_ > relativeLocalToGlobalTimeStepMinimum_ * deltaT) {
                deltaTsubcycleNext_ = relativeLocalToGlobalTimeStepMinimum_ * deltaT;
            }
        }
        /*
         *  Find the deltaTsubcycleNext_ for delta jacobians calculations 
         */
       	if (subIntegrationType == FVDELTA_TIMEINTEGRATION_SIR || subIntegrationType == FIXEDSUBCYCLENUMBER_TIMEINTEGRATION_SIR) {
            // Check this out
	    double tbase = timeHistory_base_.getNextRegularTime(tfinal_);
	    if (tfinal_ > tbase * 1.00000000001) {
		timeHistory_base_.advanceTimeStepCounter();
	    }
	    double tfinalnew = timeHistory_base_.getNextRegularTime(tfinal_);
	    double dtNew = tfinalnew - tfinal_;
	    
	    // we may compare the calculated deltaTnew with the base deltaTnew here if we want.
	    if (dtNew >  deltaTsubcycleNext_ * 1.001) {
		//printf("Warning\n");
	    }
	    // Right now, let's just assign it.
	    deltaTsubcycleNext_ = dtNew;
	}

        /*
         *  If this is the first subcycle time step, calculate the first time step for the next global iteration
         */
        if ((choiceDeltaTsubcycle_init_ == 1) && (iterSubCycle == 1)) {
            deltaTsubcycle_init_next_ = MIN(deltaTsubcycle_init_next_, deltaTsubcycleNext_);
        }
        /*
         * Signal that we have a good solnDot_init_ to work with
         */
        haveGood_solnDot_init_ = true;
        /*
         * Determine if we are at the end of the global time step
         */
        if ((1.0000000000001) * tfinal_ >= t_init_init_ + deltaT) {
            notDone = false;
            tfinal_ = t_init_init_ + deltaT;
        }

	/*
	 *  Determine if we have reached the end of number of allowed subcycles
         *  If we have, we adjust the value of t_final_final_ to the current time so that we return a reduced time step length
	 */
	if (iterSubCycle >= maxNumberSubCycles_) {
	    notDone = false;
	    t_final_final_ = tfinal_;
	}
	
        /*
         * Potentially write out intermediate data to a CSV file
         */
        if (printCSVLvl_ > 0) {
            if (notDone) {
                writeCSVData(1);
            } else {
                writeCSVData(2);
            }
        }
        /*
         *  Check the final state for any error conditions and possibly fix the results
         *  - one possible use of this is to enforce mass balances up to round-off error.
         */
        check_final_state();
        /*
         *  Calculate the first delta time step value for the next subgrid iteration of the next global iteration
         */
        if (choiceDeltaTsubcycle_init_ == 0) {
            deltaTsubcycle_init_next_ = deltaTsubcycleNext_;
        }
        /*
         *  Save the state of the Electrode object into the XML final spot
         */
        if (eState_save_) {
#ifdef DEBUG_CHECK_XML
            if (xmlStateData_final_) {
                bool retn = check_XML_valid(xmlStateData_final_);
                if (!retn) {
                    throw Electrode_Error("Electrode_Integrate() 3", "xmlStateData_final_ corrupted");
                }
            }
            if (xmlStateData_init_) {
                bool retn = check_XML_valid(xmlStateData_init_);
                if (!retn) {
                    throw Electrode_Error("Electrode_Integrate() 3", "xmlStateData_init_ corrupted");
                }
            }
#endif
            eState_save_->copyElectrode_intoState(this, true);
            SAFE_DELETE(xmlStateData_final_);
            xmlStateData_final_ = eState_save_->write_electrodeState_ToXML();
            //makeXML_TI_intermediate();
            addtoXML_TI_final(notDone);
#ifdef DEBUG_CHECK_XML
            if (xmlStateData_final_) {
                bool retn = check_XML_valid(xmlStateData_final_);
                if (!retn) {
                    throw Electrode_Error("Electrode_Integrate() 4", "xmlStateData_final_ corrupted");
                }
            }
            if (xmlStateData_init_) {
                bool retn = check_XML_valid(xmlStateData_init_);
                if (!retn) {
                    throw Electrode_Error("Electrode_Integrate() 4", "xmlStateData_init_ corrupted");
                }
            }
#endif
        }

        //-------------------------------- End of Subcycle ---------------------------------------------------------------

    } while (notDone);

    /*
     *  Calcualte the first time step for the next global iteration
     */
    if (choiceDeltaTsubcycle_init_ == 0) {
        deltaTsubcycle_init_next_ = MIN(deltaTsubcycle_init_next_, 2.0 * deltaT);
    }

    /*
     *  Store the base time integration strategy
     */
    if (subIntegrationType == BASE_TIMEINTEGRATION_SIR || subIntegrationType == BASE_REEVALUATION_TIMEINTEGRATION_SIR) {
	setTimeHistoryBaseFromCurrent();
    }
    /*
     *  Copy the results into the final_final holding pens.
     */
    setFinalFinalStateFromFinal();
    numIntegrationSubCycles_final_final_ = iterSubCycle;

    return iterSubCycle;
}
//==================================================================================================================
//   Calculate the largest mole fraction in each of the phases
/*
 *  We use this to determine the equation system
 */
void Electrode_Integrator::determineBigMoleFractions()
{
}
//==================================================================================================================
// Check the consistency of the initial state
/*
 *  (virtual from Electrode_Integrator)
 *
 *   This is a checker routine. Therefore, it can be taken out in nondebug mode
 */
void  Electrode_Integrator::check_init_consistency() const
{

}
//==================================================================================================================
//!   Prepare the problem for predicting the solution
/*!
 *      -> determine the largest mole fraction in a phase
 *      -> determine if a surface phase can have reactions turned on
 *      -> do stuff that affects the equation system that will be solved at the
 *         current step.
 */
void  Electrode_Integrator::prepareProblemStatement()
{

}
//==================================================================================================================
// Set the base tolerances for the nonlinear solver within the integrator
/*
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
 *   @param atolResid  Absolute tolerance cutoff for the electron source term over the interval.
 *                     Below this value we do not care about the results.
 *                     atol has units of kmol.
 */
void  Electrode_Integrator::setNLSGlobalSrcTermTolerances(double rtolResid)
{
    throw Electrode_Error("Electrode_Integrator::setNLSGlobalSrcTermTolerances","unimplemented");
}
//======================================================================================================================
void  Electrode_Integrator::zeroGlobalStepAccumulationTerms()
{
    /*
     *  Zero the global error vectors
     */
    std::fill(errorGlobalNLS_.begin(), errorGlobalNLS_.end(), 0.0);
    /*
     *   Zero the integrated source term - This will be accumulated over the subcycling
     */
    std::fill(spMoleIntegratedSourceTerm_.begin(), spMoleIntegratedSourceTerm_.end(), 0.0);
    std::fill(spMoleIntegratedSourceTermLast_.begin(), spMoleIntegratedSourceTermLast_.end(), 0.0);
    /*
     *  Zero the ThermalEnergySourceTerm
     */
    integratedThermalEnergySourceTerm_ = 0.0;
    integratedThermalEnergySourceTermLast_ = 0.0;
    integratedThermalEnergySourceTerm_overpotential_ = 0.0;
    integratedThermalEnergySourceTerm_overpotential_Last_ = 0.0;
    integratedThermalEnergySourceTerm_reversibleEntropy_ = 0.0;
    integratedThermalEnergySourceTerm_reversibleEntropy_Last_ = 0.0;
}
//==================================================================================================================
size_t Electrode_Integrator::numIntegratedSrcTerms() const
{
    return numIntegratedSrc_;
}
//==================================================================================================================
//   Set the Residual absolute error tolerances
/*
 *  (virtual from Electrode_Integrator)
 *
 *   Set the absolute error tolerances for the residuals for the nonlinear solvers. This is called at the top
 *   of the integrator() routine.
 *
 *   Calculates residAtolNLS_[]
 *   Calculates atolNLS_[]
 */
void  Electrode_Integrator::setResidAtolNLS()
{
    throw Electrode_Error("Electrode_Integrator::setResidAtolNLS()", "unimplemented");
}
//==================================================================================================================================
int Electrode_Integrator::predictSoln()
{
    throw Electrode_Error(" Electrode_Integrator::predictSoln()","unimplemented");
}
//==================================================================================================================================
int Electrode_Integrator::predictSolnDot()
{
    int info = 1;
    soln_predict_fromDot_[0] = deltaTsubcycleCalc_;
    for (size_t i = 1; i < yvalNLS_.size(); i++) {
	soln_predict_fromDot_[i] = yvalNLS_init_[i] + deltaTsubcycleCalc_ * solnDot_init_[i];
    }
    // We only use the prediction if the predictSolnDot method is deemed the best predictor!
    if (predictDotBetter_) {
        // HKM -> This can happen if the prediction is bad. It must be checked.
	info = unpackNonlinSolnVector(&soln_predict_fromDot_[0]);
        if (info != 1) {
            predictDotBetter_ = false;
            // Go back to the previous soln_predict_[] created from the previous estimation method to get an initial solution.
            info = unpackNonlinSolnVector(&soln_predict_[0]);
            if (info != 1) {
                return -1;
            }
        }
	updateState();
	extractInfo();
	updateSpeciesMoleChangeFinal();
    }
    return info;
}
//==================================================================================================================================
void  Electrode_Integrator::extractInfo()
{
    throw Electrode_Error("Electrode_Integrator::extractInfo()",  "unimplemented");
}
//==================================================================================================================================
void Electrode_Integrator::updateSpeciesMoleChangeFinal()
{
    throw Electrode_Error("Electrode_Integrator::updateSpeciesMoleChangeFinal()",  "unimplemented");
}
//==================================================================================================================================
int Electrode_Integrator::calcResid(double* const resid, const ResidEval_Type evalType)
{
    throw Electrode_Error("Electrode_Integrator::calcResid()",  "unimplemented");
    return 0;
}
//==================================================================================================================================
int Electrode_Integrator::GFCEO_evalResidNJ(const double t, const double delta_t, const double* const y, const double* const ydot,
                                            double* const resid, const ResidEval_Type evalType,
                                            const int id_x, const double delta_x)
{
    throw Electrode_Error("Electrode_Integrator::GFCE_evalResidNJ()",  "unimplemented");
    return 0;
}
//==================================================================================================================================
int Electrode_Integrator::GFCEO_calcResid(double* const resid, const ResidEval_Type evalType)
{
    throw Electrode_Error("Electrode_Integrator::GFCE_calcResid()",  "unimplemented");
    return 0;
}
//==================================================================================================================================
void Electrode_Integrator::gatherIntegratedSrcPrediction()
{
    throw Electrode_Error("Electrode_Integrator::gatherIntegratedSrcPrediction()",  "unimplemented");
}
//==================================================================================================================================
void Electrode_Integrator::calcSrcTermsOnCompletedStep()
{
    /*
     *  Calculate the integrated source term
     *       An alternative would be to redo the residual calculation. However, here we assume that
     *       the residual calculation has been done and the results are in _final_
     */
    for (size_t i = 0; i < m_NumTotSpecies; i++) {
        spMoleIntegratedSourceTermLast_[i] = spMoles_final_[i] - spMoles_init_[i];
    }
    if (doThermalPropertyCalculations_) {
        integratedThermalEnergySourceTermLast_ = thermalEnergySourceTerm_EnthalpyFormulation_SingleStep();
        /*
         * these last two are needed for informative output only
         */
        integratedThermalEnergySourceTerm_overpotential_Last_ = thermalEnergySourceTerm_Overpotential_SingleStep();
        integratedThermalEnergySourceTerm_reversibleEntropy_Last_ = thermalEnergySourceTerm_ReversibleEntropy_SingleStep();
    }
    if (doPolarizationAnalysis_) {
        polarSrc_list_Last_.clear();
        double icurrAccount = polarizationAnalysisSurf(polarSrc_list_Last_);
        double icurrLast = - spMoleIntegratedSourceTermLast_[kElectron_];
        if (fabs (icurrAccount - icurrLast) < 1.0E-10) {
            throw Electrode_Error("Electrode_Integrator::calcSrcTermsOnCompletedStep()",
                                  " error in accounting for electrons");
        }
    }
}
//==================================================================================================================================
void Electrode_Integrator::accumulateSrcTermsOnCompletedStep(bool remove)
{
    if (remove) {
        for (size_t i = 0; i < m_NumTotSpecies; i++) {
            spMoleIntegratedSourceTerm_[i] -= spMoleIntegratedSourceTermLast_[i];
        }
	if (doThermalPropertyCalculations_) {
	    integratedThermalEnergySourceTerm_ -= integratedThermalEnergySourceTermLast_;
	}
    } else {
        for (size_t i = 0; i < m_NumTotSpecies; i++) {
            spMoleIntegratedSourceTerm_[i] += spMoleIntegratedSourceTermLast_[i];
        }
	if (doThermalPropertyCalculations_) {
	    integratedThermalEnergySourceTerm_ += integratedThermalEnergySourceTermLast_;
            integratedThermalEnergySourceTerm_overpotential_ += integratedThermalEnergySourceTerm_overpotential_Last_;
            integratedThermalEnergySourceTerm_reversibleEntropy_ += integratedThermalEnergySourceTerm_reversibleEntropy_Last_;
	}
        if (doPolarizationAnalysis_) {
              integratedPolarizationCalc();
        }
    }
}
//==================================================================================================================================
void  Electrode_Integrator::check_yvalNLS_init(bool doOthers)
{
}
//==================================================================================================================================
void Electrode_Integrator::setInitStateFromFinal(bool setInitInit)
{
    Electrode::setInitStateFromFinal(setInitInit);
    if (solnDot_init_.empty()) {
	create_solvers();
    }
    // We carry over the time derivative of the solution from the old step to the new step here.
    for (size_t i = 0; i < neq_; i++) {
	solnDot_init_[i] = solnDot_final_[i];
	yvalNLS_init_[i] = yvalNLS_[i];
    }
    if (setInitInit) {
	for (size_t i = 0; i < neq_; i++) {
	    solnDot_init_init_[i] =  solnDot_final_[i];
	    yvalNLS_init_init_[i] =  yvalNLS_[i];
	}
    }
}
//==================================================================================================================================
void Electrode_Integrator::setInitInitStateFromFinalFinal()
{
    Electrode::setInitInitStateFromFinalFinal();

    if( solnDot_init_.empty() ) create_solvers();
    for (size_t i = 0; i < neq_; i++) {
	solnDot_init_init_[i] = solnDot_final_final_[i];
	solnDot_init_[i]      = solnDot_final_final_[i];
	solnDot_final_[i]     = solnDot_final_final_[i];
	yvalNLS_init_init_[i] = yvalNLS_final_final_[i];
	yvalNLS_init_[i]      = yvalNLS_final_final_[i];
	yvalNLS_[i]           = yvalNLS_final_final_[i];
    }
}
//==================================================================================================================
// Set the internal final intermediate and from the internal init state
/*
 *  (virtual function from Electrode)
 *
 *  Set the final state from the init state. This is commonly called during a failed time step
 *
 */
void Electrode_Integrator::setFinalStateFromInit()
{
    Electrode::setFinalStateFromInit();
    int neqNLS = nEquations();
    if (solnDot_init_.empty()) {
	create_solvers();
    }
    for (int i = 0; i < neqNLS; i++) {
	solnDot_final_[i] = solnDot_init_[i];
	yvalNLS_[i] = yvalNLS_init_[i];
    }
}
//==================================================================================================================================
void Electrode_Integrator::setInitStateFromInitInit(bool setFinal)
{
    size_t i, neqNLS;
    Electrode::setInitStateFromInitInit(setFinal);
    neqNLS = nEquations();
    if (solnDot_init_.empty()) {
	create_solvers();
    }
    for (i = 0; i < neqNLS; i++) {
	solnDot_init_[i] = solnDot_init_init_[i];
	yvalNLS_init_[i] = yvalNLS_init_init_[i];
    }
    if (setFinal) {
        for (i = 0; i < neqNLS; i++) {
            solnDot_final_[i] = solnDot_init_init_[i];
	    yvalNLS_[i] = yvalNLS_init_init_[i];
        }
    }
}
//==================================================================================================================================
void Electrode_Integrator::setFinalFinalStateFromFinal()
{
    size_t i, neqNLS;
    Electrode::setFinalFinalStateFromFinal();
    neqNLS = nEquations();
    if (solnDot_init_.empty() ) {
         create_solvers();
    }
    for (i = 0; i < neqNLS; i++) {
	 solnDot_final_final_[i] = solnDot_final_[i];
	 yvalNLS_final_final_[i] = yvalNLS_[i];
     }
 }
//==================================================================================================================================
int Electrode_Integrator::calcDeltaSolnVariables(const double t, const double* const ySoln, const double* const ySolnDot, 
                                                 double* const deltaYSoln, const double* const solnWeights)
{
    int retn = ResidJacEval::calcDeltaSolnVariables(t, ySoln, ySolnDot, deltaYSoln, solnWeights);
    return retn;
}
//==================================================================================================================================
int Electrode_Integrator::unpackNonlinSolnVector(const double* const y)
{
    throw Electrode_Error(" Electrode_Integrator::unpackNonlinSolnVector()", "unimplemented");
    return 1;
}
//==================================================================================================================================
bool Electrode_Integrator::checkSubIntegrationStepAcceptable() const
{
    return true;
}
//==================================================================================================================================
void Electrode_Integrator::calc_solnDot_final()
{
    const int m_order = 1;
    doublevalue c1;
    switch (m_order) {
    case 0:
        for (size_t i = 0; i < neq_; i++) {
	    solnDot_final_[i] = c1 * (yvalNLS_[i] - yvalNLS_init_[i]); 
        }
    case 1:             /* First order forward Euler/backward Euler */
        c1 = 1.0 / deltaTsubcycleCalc_;
        for (size_t i = 0; i < neq_; i++) {
	    solnDot_final_[i] = c1 * (yvalNLS_[i] - yvalNLS_init_[i]); 
        }
        return;
    case 2:             /* Second order Adams-Bashforth / Trapezoidal Rule */
        c1 = 2.0 / deltaTsubcycleCalc_;
        for (size_t i = 0; i < neq_; i++) {
            solnDot_final_[i] = c1 * (yvalNLS_[i] - yvalNLS_init_[i]) - solnDot_init_[i];
        }
        return;
    default:
        throw Electrode_Error("calc_solnDot_final()", "Case not covered");
    }

    solnDot_final_[0] = 0.0;
    if (!haveGood_solnDot_init_) {
        for (size_t i = 0; i < neq_; i++) {
            solnDot_init_[i] = solnDot_final_[i]; 
        }
    }
}
//==================================================================================================================================
// This routine must agree with calc_ydot() within the Zuzax' nonlinear solver, NonlinearSolver
// ydotNLS_[] is updated within the nonlinear solver as yvalNLS_ is updated. This routine serves to initialize the values on
// going into the nonlinear solver.
double Electrode_Integrator::calc_ydotNLS_final()
{
    const int m_order = 1;
    doublevalue c1;
    switch (m_order) {
    case 0:             
        c1 = 0.0;
        for (size_t i = 0; i < neq_; i++) {
            ydotNLS_[i] = solnDot_init_[i];
        }
        break;
    case 1:             /* First order forward Euler/backward Euler */
        c1 = 1.0 / deltaTsubcycleCalc_;
        for (size_t i = 0; i < neq_; i++) {
            ydotNLS_[i] = c1 * (yvalNLS_[i] - yvalNLS_init_[i]);
        }
        break;
    case 2:             /* Second order Adams-Bashforth / Trapezoidal Rule */
        c1 = 2.0 / deltaTsubcycleCalc_;
        for (size_t i = 0; i < neq_; i++) {
            ydotNLS_[i] = c1 * (yvalNLS_[i] - yvalNLS_init_[i]) - solnDot_init_[i];
        }
        break;
    default:
        throw Electrode_Error("calc_ydotNLS_final()", "Case not covered");
    }
    return c1;
}
//==================================================================================================================================
double Electrode_Integrator::predictorCorrectorWeightedSolnNorm(const std::vector<double>& yvalNLS)
{
    double pnorm     = l0norm_PC_NLS(soln_predict_,         yvalNLS, neq_, atolNLS_, rtolNLS_);
    double pnorm_dot = l0norm_PC_NLS(soln_predict_fromDot_, yvalNLS, neq_, atolNLS_, rtolNLS_);
    if (pnorm_dot < pnorm) {
       pnorm = pnorm_dot;
    } else {
       // Have to redo this again, in order to assign the errorLocalNLS_[] back to the first way to solution prediction method
       pnorm         = l0norm_PC_NLS(soln_predict_,         yvalNLS, neq_, atolNLS_, rtolNLS_);
    }
    return pnorm;
}
//==================================================================================================================================
// On a completed local step, accumulate local errors into the global error vectors for the global time step
/*
 */
void Electrode_Integrator::accumulateLocalErrorToGlobalErrorOnCompletedStep()
{
    for (size_t k = 0; k < neq_; k++) {
        errorGlobalNLS_[k] += errorLocalNLS_[k];
    }
}
//==================================================================================================================================
// When the norm is high this routine decides what to do with the step
/*
 *  @param  pnorm  Norm of the step
 *
 *  @return   Returns whether the step is accepted or not
 */
bool Electrode_Integrator::decide_normHighLogic(double pnorm, int num_newt_its, int iterSubCycle)
{
    double deltaT = t_final_final_ - t_init_init_;
    /*
     *  The fundamental decision is to reject the step out of time step truncation error if pnorm is greater than 2
     */
    if (pnorm > 2.0) {
	/*
	 *  Calculate the deltaT below which we will not check for errors
	 */
	double deltaT_min = deltaT /  maxNumberSubGlobalTimeSteps_;
	double timeLeft =  t_final_final_ - tinit_ + 1.0E-20;
	double itsLeft = maxNumberSubGlobalTimeSteps_ - iterSubCycle + 1;
	double deltaT_min_2 = deltaT;
	if (itsLeft > 0) {
	    deltaT_min_2 = timeLeft / itsLeft;
	}
	double deltaT_min_3 = relativeLocalToGlobalTimeStepMinimum_ * deltaT * 2.0;
	deltaT_min = std::max( deltaT_min, deltaT_min_2);
	deltaT_min = std::max( deltaT_min , deltaT_min_3);
	/*
	 *  If we are below the deltaT subcycle number limit we will accept the step anyway.
	 *  (note if there are convergence issues we should modify this)
	 */
        if ((deltaTsubcycle_ > deltaT_min) || num_newt_its >= 10) {
            if (choiceDeltaTsubcycle_init_ != 2) {
                if (printLvl_ > 2) {
                    printf("\t\tdecide_normHighLogic(): WARNING Step rejected because time step"
			   " truncation error norm, %g, is greater than 2\n", pnorm);
                }
                return false;
            }
        } else {
	    if (printLvl_ > 2) {
		printf("\t\tdecide_normHighLogic(): WARNING Step accepted, time step"
		       " truncation error norm, %g, but deltaTsubcycle, %g is below deltat limit, %g\n",
		       pnorm, deltaTsubcycle_, deltaT_min);
	    }
	}
    } else if (pnorm > 1.0) {
	if (num_newt_its >= 10) {
	    if (choiceDeltaTsubcycle_init_ != 2) {
                if (printLvl_ > 2) {
		    printf("\t\tdecide_normHighLogic(): WARNING Step rejected because time step"
			   " truncation error norm, %g, is greater than 1 and newt its are high, %d\n", pnorm, num_newt_its);
                  
                }
                return false;
            }
	}
    }
    return true;
}
//====================================================================================================================
// Calculate the vector of predicted errors in the source terms that this integrator is responsible for
/*
 *  (virtual from Electrode_Integrator)
 *
 *    In the base implementation we assume that the there are just one source term, the electron
 *    source term.
 *    However, this will be wrong in almost all cases.
 *    The number of source terms is unrelated to the number of unknowns in the nonlinear problem.
 *    Source terms will have units associated with them.
 *    For example the integrated source term for electrons will have units of kmol
 */
void Electrode_Integrator::predictorCorrectorGlobalSrcTermErrorVector()
{
}
//====================================================================================================================
//  Calculate the normalized norm of the errors in the global source terms
/*
 *  (virtual from Electrode_Integrator)
 *
 *   This routine make use of the source term error vector along with rtols and atols for the
 *   individual source terms to calculated a normalized error measure. This is the single number
 *   that the integration routine will try to control as it calculates a time stepping strategy.
 *
 *   @return  Returns a single nondimensional number representing the normalized error
 *            for the calculation of the source term
 */
double Electrode_Integrator::predictorCorrectorGlobalSrcTermErrorNorm()
{
    return 0.0;
}
//==================================================================================================================================
/*
 *  @param yval           Vector of corrector values
 *  @param pnormSrc       Norm of the predictor-corrector comparison for the source vector.
 *  @param pnormSoln      Norm of the predictor-corrector comparison for the solution vector.
 */
void Electrode_Integrator::predictorCorrectorPrint(const std::vector<double>& yval, double pnormSrc, double pnormSoln) const
{
    double denom, tmp, diff1, diff2, sdiff;
    printf(" -------------------------------------------------------------------------------------------------------------------\n");
    printf(" PREDICTOR_CORRECTOR  SubIntegrationCounter = %7d       t_init = %12.5E,       t_final = %12.5E       deltaT = %12.5E\n",
           counterNumberSubIntegrations_, tinit_, tfinal_, tfinal_ - tinit_);
    printf("                         IntegrationCounter = %7d  t_init_init = %12.5E, t_final_final = %12.5E   deltaTGlob = %12.5E\n",
           counterNumberIntegrations_, t_init_init_, t_final_final_, t_final_final_ - t_init_init_);
    printf(" -------------------------------------------------------------------------------------------------------------------\n");
    printf("             Initial |     Prediction   Prediction_Dot    Actual   |   SDifference   |    ATol   |"
           "   Contrib  | RTOL = %g\n", rtolNLS_);

    denom =  rtolNLS_ * MAX(fabs(yval[0]), fabs(soln_predict_[0]));
    denom = MAX(denom, atolNLS_[0]);
    tmp = fabs((yval[0] - soln_predict_[0])/ denom);
    diff1 = yval[0] - soln_predict_[0];
    diff2 = yval[0] - soln_predict_fromDot_[0];
    sdiff = (fabs(diff2) < fabs(diff1)) ? diff2 : diff1;
    printf(" DeltaT   |% 14.7E | % 14.7E % 14.7E % 14.7E | % 14.7E | % 10.3E | % 10.3E |\n",
           deltaTsubcycleCalc_, soln_predict_[0], soln_predict_fromDot_[0], yval[0], sdiff, atolNLS_[0], tmp);
    for (int i = 1; i < (int) yval.size(); i++) {
        denom = rtolNLS_ * MAX(fabs(yval[i]), fabs(soln_predict_[i]));
        denom = MAX(denom, atolNLS_[i]);
        diff1 = yval[i] - soln_predict_[i];
        diff2 = yval[i] - soln_predict_fromDot_[i];
        sdiff = (fabs(diff2) < fabs(diff1)) ? diff2 : diff1;
        tmp = fabs((yval[i] - soln_predict_[i])/ denom);
        printf(" soln %3d | % 14.7E | % 14.7E % 14.7E % 14.7E | % 14.7E | % 10.3E | % 10.3E | \n",
               i, yvalNLS_init_[i], soln_predict_[i], soln_predict_fromDot_[i], yval[i], sdiff, atolNLS_[i], tmp);
    }
    printf(" -----------------------------------------------------------------------------------------------------"
           "--------------\n");
    printf("                                                                                        %10.3E\n",
           pnormSoln);
}
//==================================================================================================================================
bool Electrode_Integrator::changeSolnForBirthDeaths()
{
    return true;
}
//==================================================================================================================================
void Electrode_Integrator::manageBirthDeathSuccessfulStep()
{
}
//==================================================================================================================================
void Electrode_Integrator::check_final_state()
{
#ifdef DEBUG_NEWMODELS
    // While in debug mode, check the inventory of capacity to account for all electrons
    checkCapacityBalances_final();
#endif
}
//==================================================================================================================================
void Electrode_Integrator::printResid_TimeHeader()
{
}
//==================================================================================================================================
int  Electrode_Integrator::residEval_BaseChecks()
{
    return 1;
}
//==================================================================================================================================
void  Electrode_Integrator::printResid_ResidSatisfaction()
{
}
//==================================================================================================================================
void Electrode_Integrator::initialPackSolver_nonlinFunction()
{
    throw Electrode_Error(" Electrode_Integrator::initialPackSolver_nonlinFunction()",  "unimplemented");
}
//==================================================================================================================================
int Electrode_Integrator::check_nonlinResidConditions()
{
    return 0;
}
//==================================================================================================================================
double Electrode_Integrator::reportStateVariableIntegrationError(int& maxSV, double* const errorVector) const
{
    maxSV = 0;
    double vmax = errorGlobalNLS_[0];
    if (errorVector) {
        errorVector[0] = errorGlobalNLS_[0];
        for (size_t k = 1; k < neq_; k++) {
            errorVector[k] = errorGlobalNLS_[k];
            if (errorVector[k] > vmax) {
                maxSV = k;
                vmax = errorVector[k];
            }
        }
    } else {
        for (size_t k = 1; k < neq_; ++k) {
            if (errorGlobalNLS_[k] > vmax) {
                maxSV = k;
                vmax = errorGlobalNLS_[k];
            }
        }
    }
    return vmax;
}
//==================================================================================================================================
int Electrode_Integrator::evalResidNJ(const double t, const double delta_t,
                                      const double* const y, const double* const ySolnDot,
                                      double* const resid,
                                      const ResidEval_Type evalType, const int id_x,
                                      const double delta_x)
{
    int retn = 1;
    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
        printf("\t\t=========================================================================================================="
               "=====================\n");
        printf("\t\t  EXTRA PRINTING FROM NONLINEAR RESIDUAL: ");
        if (evalType ==  ResidEval_Type::Base_ResidEval) {
            printf(" BASE RESIDUAL");
        } else if (evalType == ResidEval_Type::JacBase_ResidEval) {
            printf(" BASE JAC RESIDUAL");
        } else  if (evalType == ResidEval_Type::JacDelta_ResidEval) {
            printf(" DELTA JAC RESIDUAL");
            printf(" var = %d delta_x = %12.4e Y_del = %12.4e Y_base = %12.4e", id_x, delta_x, y[id_x], y[id_x] - delta_x);
        } else  if (evalType == ResidEval_Type::Base_ShowSolution) {
            printf(" BASE RESIDUAL - SHOW SOLUTION");
        }
        printf(" DomainNumber = %2d , CellNumber = %2d , IntegrationCounter = %d\n",
               electrodeDomainNumber_, electrodeCellNumber_, counterNumberIntegrations_);
    }
    /*
     *  UNPACK THE SOLUTION VECTOR
     */
    retn = unpackNonlinSolnVector(y);
    if (retn != 1) {
        if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
          printf(" Bounds problem unpacking solution vector. Bailing on residual calculation\n");
          return -1;
        }
    }

    if (evalType != ResidEval_Type::JacDelta_ResidEval && (evalType != ResidEval_Type::Base_LaggedSolutionComponents)) {
        //mdp::mdp_copy_dbl_1(DATA_PTR(phaseMoles_final_lagged_),(const double *)DATA_PTR(phaseMoles_final_), m_NumTotPhases);
    }
    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
        printResid_TimeHeader();

    }

    /*
     *  Take the unpacked solution and calculate the consistent and full final state
     */
    updateState();

    /*
     *  Query Zuzax for all of the rate information at the final state (and the initial state if we are doing higher order)
     */
    extractInfo();

    /*
     *   We take the ROP_inner_[] and ROP_outer_[] rates of progress, combine them with the surface area calculation,
     *   and the stoichiometric coefficients to calculate the DspMoles_final_[], which is production rate for
     *   all species in the electrode.
     */
    updateSpeciesMoleChangeFinal();

    /*
     * Calculate the residual
     */
    calcResid(resid, evalType);

    /*
     * Change the problem specification for the nonlinear solve in certain circumstances when the solver
     *  is calculating the base residual of a jacobian
     */
    if (evalType == ResidEval_Type::JacBase_ResidEval) {
        /*
         *  Perform some checks that would lead to the return flag indicating an error condition
         */
        retn = residEval_BaseChecks();
    }


    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
        printResid_ResidSatisfaction();
    }
    int rr = 1;
    if (retn < 0) {
        rr = retn;
    }
    return rr;
}
//==================================================================================================================================
void Electrode_Integrator::printElectrode(int pSrc, bool subTimeStep)
{
    int iph;
    double egv = TotalVol();
    double tsm = SolidTotalMoles();
    printf("   ==============================================================================================\n");
    if (subTimeStep) {
        printf("      Electrode at intermediate-step time final = %12.5E\n", tfinal_);
        printf("                   intermediate-step time init  = %12.5E     deltaT = %12.5E\n",
               tinit_, deltaTsubcycleCalc_);
        printf("                ChemModel Type = %3d , DomainNumber = %2d , CellNumber = %2d , SubIntegrationCounter = %d\n",
               electrodeChemistryModelType_, electrodeDomainNumber_, electrodeCellNumber_, counterNumberSubIntegrations_);
    } else {
        printf("      Electrode at time final_final = %12.5E\n", t_final_final_);
        printf("                   time init_init   = %12.5E        deltaTglobal = %12.5E\n", t_init_init_,
               t_final_final_ - t_init_init_);
        printf("                ChemModel Type = %3d , DomainNumber = %2d , CellNumber = %2d , IntegrationCounter = %d\n",
               electrodeChemistryModelType_, electrodeDomainNumber_, electrodeCellNumber_, counterNumberIntegrations_);
	printf("                numIntegrationSubCycles = %d, SubIntegrationCounter = %d\n",
	       static_cast<int>(numIntegrationSubCycles_final_final_), counterNumberSubIntegrations_);
    }
    printf("   ==============================================================================================\n");
    printf("\n");
    printf("          Number of external surfaces = %d\n", (int) numExternalInterfacialSurfaces_);
    printf("          Solid Volume = %11.4E m**3\n", ElectrodeSolidVolume_);
    printf("          Total Volume = %11.4E m**3\n", egv);
    if (egv > 0.0) {
        printf("          Porosity     = %10.3E\n", (egv - ElectrodeSolidVolume_) / egv);
    } else {
        printf("          Porosity     =       NA\n");
    }
    printf("          Temperature = %g K\n", temperature_);
    printf("          Pressure = %g Pa\n", pressure_);
    printf("          Total Solid Moles = %11.4E kmol\n", tsm);
    if (tsm > 0.0) {
        printf("          Molar Volume of Solid = %11.4E cm3 gmol-1\n",  ElectrodeSolidVolume_ / tsm * 1.0E3);
    } else {
    }
    printf("          Particle Number to Follow = %11.4E\n", particleNumberToFollow_);
    if (subTimeStep) {
        double curr = integratedLocalCurrent();
        printf("          Current = %12.5E amps\n", curr);
    } else {
        double curr = integratedCurrent();
        printf("          Current = %12.5E amps\n", curr);
    }
    printf("          Voltage (phiMetal - phiElectrolyte) = %12.5E volts\n", deltaVoltage_);

    printf("          followElectrolyteMoles = %d\n", followElectrolyteMoles_);
    printf("          ElectrolytePseudoMoles = %g\n", electrolytePseudoMoles_);
    printf("\n");

    printElectrodeCapacityInfo(pSrc, subTimeStep);
    printf("\n");
    printElectrodePhaseList(pSrc, subTimeStep);

    int m = m_NumTotPhases;
    if ((size_t) numSurfaces_ > m_NumSurPhases) {
        m = numSurfaces_ + m_NumVolPhases;
    }
    for (iph = 0; iph < m; iph++) {
        printElectrodePhase(iph, pSrc, subTimeStep);
    }

}
//==================================================================================================================================
void Electrode_Integrator::printElectrodePhase(size_t iph, int pSrc, bool subTimeStep)
{
    size_t isph = npos;
    double* netROP = new double[m_NumTotSpecies];
    ThermoPhase& tp = thermo(iph);
    std::string pname = tp.id();
    size_t istart = m_PhaseSpeciesStartIndex[iph];
    size_t nsp = tp.nSpecies();
    printf("     ===============================================================\n");
    printf("          Phase %d %s \n", static_cast<int>(iph), pname.c_str());
    printf("                Total moles = %g\n", phaseMoles_final_[iph]);
    if (iph == metalPhase_) {
        double deltaT = t_final_final_ - t_init_init_;
        if (subTimeStep) {
            deltaT = tfinal_ - tinit_;
        }
        if (deltaT > 1.0E-200) {
            double amps = spMoleIntegratedSourceTerm_[istart] / deltaT * Faraday;
            if (subTimeStep) {
                amps = spMoleIntegratedSourceTermLast_[istart] / deltaT * Faraday;
            }
            printf("                Current = %g amps \n", amps);
        } else {
            printf("                Current = NA amps \n");
        }
    }
    if (iph == metalPhase_ || iph == solnPhase_) {
        printf("                  Voltage = %g\n", tp.electricPotential());
    }
    if (iph >= m_NumVolPhases) {
        isph = iph - m_NumVolPhases;
        printf("                surface area (final) = %11.5E m2\n",  surfaceAreaRS_final_[isph]);
        printf("                surface area (init)  = %11.5E m2\n",  surfaceAreaRS_init_[isph]);
        int ddd =  isExternalSurface_[isph];
        printf("                IsExternalSurface = %d\n", ddd);
        double oc = openCircuitVoltage(isph);
        if (oc != 0.0) {
            printf("                 Open Circuit Voltage = %g\n", oc);
        }
    }
    printf("\n");
    printf("                Name               MoleFrac_final  kMoles_final kMoles_init SrcTermLastStep(kMoles)\n");
    for (size_t k = 0; k < nsp; k++) {
        std::string sname = tp.speciesName(k);
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
    if (iph >= m_NumVolPhases) {
        const std::vector<double>& rsSpeciesProductionRates = RSD_List_[isph]->calcNetSurfaceProductionRateDensities();
        RSD_List_[isph]->getNetRatesOfProgress(netROP);

        double* spNetProdPerArea = (double*) spNetProdPerArea_List_.ptrColumn(isph);
        std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
        size_t nphRS = RSD_List_[isph]->nPhases();
        size_t kIndexKin = 0;
        for (size_t kph = 0; kph < nphRS; kph++) {
            size_t jph = RSD_List_[isph]->KinToPL_PhaseIndex_[kph];
            size_t istart = m_PhaseSpeciesStartIndex[jph];
            size_t nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
            for (size_t k = 0; k < nsp; ++k) {
                spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
                kIndexKin++;
            }
        }
        printf("\n");
        printf("                           spName                  Source (kmol/m2/s) \n");
        for (size_t k = 0; k <  m_NumTotSpecies; k++) {
            std::string ss = speciesName(k);
            printf("                           %-22s %10.3E\n", ss.c_str(), spNetProdPerArea[k]);
        }
    }
    delete [] netROP;
}
//==================================================================================================================================
double Electrode_Integrator::l0norm_PC_NLS(const std::vector<double>& v1, const std::vector<double>& v2, size_t num,
                                           const std::vector<double>& atolVec, const double rtol)
{
    double max0 = 0.0;
    double denom, diff, ee;
    for (size_t k = 0; k < num; k++) {
        diff = fabs(v1[k] - v2[k]);
        denom = rtol * MAX(fabs(v1[k]), fabs(v2[k]));
        denom = MAX(denom, atolVec[k]);
        ee = diff / denom;
        if (ee > max0) {
            max0 = ee;
        }
        // Unnormalize the error so that the real relative level of error is reflected in the storred quantity
        errorLocalNLS_[k] = ee * rtol;
    }
    return max0;
}
//==================================================================================================================================
void Electrode_Integrator::setState_EState(const EState& es)
{
    Electrode::setState_EState(es);
    if (es.solnDot_.size() >= solnDot_final_.size()) {
        solnDot_final_ = es.solnDot_;
        solnDot_init_ = es.solnDot_;
        solnDot_init_init_ = es.solnDot_;
        solnDot_final_final_ = es.solnDot_;
    }
}
//==================================================================================================================================
int Electrode_Integrator::setTimeHistoryBaseFromCurrent()
{
    timeHistory_base_ = timeHistory_current_;
    return timeHistory_base_.nTimeStepsRegular_;
}
//==================================================================================================================================
int Electrode_Integrator::setTimeHistoryBase(const SubIntegrationHistory &timeHistory)
{
    timeHistory_base_ = timeHistory;
    return timeHistory.nTimeStepsRegular_;
}
//==================================================================================================================================
void Electrode_Integrator::setMaxNumberSubCycles(int maxN)
{
    maxNumberSubCycles_ = maxN;
}
//==================================================================================================================================
void Electrode_Integrator::setMaxNumberSubGlobalTimeSteps(int maxN)
{   
    // default number is 1000
    maxNumberSubGlobalTimeSteps_ = maxN;
}
//==================================================================================================================================
SubIntegrationHistory& Electrode_Integrator::timeHistory(bool returnBase)
{
    if (returnBase) {
	return timeHistory_base_;
    }
    return timeHistory_current_;
}
//==================================================================================================================================
} // End of namespace
//----------------------------------------------------------------------------------------------------------------------------------
