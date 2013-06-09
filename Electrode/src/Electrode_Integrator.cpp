/*
 * $Id: Electrode_Integrator.cpp 591 2013-05-09 22:06:20Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"



#include "cantera/base/mdp_allo.h"

#include "Electrode_Integrator.h"
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

namespace Cantera
{

//======================================================================================================================
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
//======================================================================================================================
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
//======================================================================================================================
SubIntegrationHistory::~SubIntegrationHistory()
{

}
//======================================================================================================================
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
//======================================================================================================================
void  
SubIntegrationHistory::clear()
{
    nTimeStepsRegular_ = 0;
    nTimeSteps_ = 0;
    TimeStepList_.clear();
    GsolnErrorNorm_ = 0.0;
    GwdotErrorNorm_ = 0.0;
    iCounter = 0;
    time_step_next = 0.0;
}
//======================================================================================================================

//! Add a step to the time history
void  
SubIntegrationHistory::addTimeStep(double t_init, double t_final, double t_final_calc, int timeTypeSoln,int numNonLinSolves,
				   double solnErrorNorm, double wdotErrorNorm, double volts,
				   double srcElectronsStep, double currentStep)
{
    /*
     * check for consistency
     */
    if (iCounter > 0) {
	TimeStepHistory& tshPrevious =  TimeStepList_[iCounter-1];
	double t_final_calc_Previous =  tshPrevious.t_final_calc_;
	if (fabs(t_init - t_final_calc_Previous) > 1.0E-13) {
	    throw Cantera::CanteraError("Shouldn't be here","");
	}
    }
    if (timeTypeSoln != 2) {
	if (fabs(t_final - t_final_calc) > 1.0E-13) {
	    throw CanteraError("Shouldn't be here", "");
	}
    } else {
	if (fabs(t_final - t_final_calc) < 1.0E-16) {
	    throw CanteraError("warning t_final == t_final_calc", "");
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
//======================================================================================================================
void  
SubIntegrationHistory::zeroTimeStepCounter()
{
    iCounter = 0;
}
//======================================================================================================================
void  
SubIntegrationHistory::advanceTimeStepCounter()
{
    if (iCounter + 1 < nTimeSteps_) {
	iCounter++;
    }
}
//======================================================================================================================
double
SubIntegrationHistory::getNextRegularTime(double currentTime) const
{
    // put in extra assignments for debugging
    //int iC = iCounter;

    if (iCounter >=  nTimeSteps_) {
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
    if (t_final_calc <= currentTime) {
	t_final_calc = currentTime + time_step_next;
    }
 
    return (t_final_calc);
}
//======================================================================================================================
int
SubIntegrationHistory::assureTimeInterval(double gtinit, double gtfinal) 
{
    int iC = 0;
    TimeStepHistory* tshCurrent_ptr = & TimeStepList_[iC];
    double tinit =  tshCurrent_ptr->t_init_;
    if (fabs(gtinit - tinit) > 1.0E-13) {
	throw CanteraError("shouldn't be here", "");
    }
    tshCurrent_ptr =  & TimeStepList_[nTimeSteps_-1];
    double tfinal =  tshCurrent_ptr->t_final_;
    if (tshCurrent_ptr->timeTypeSoln_ == 2) {
	printf("weird condition\n");
	tshCurrent_ptr->timeTypeSoln_ = 0;
    }
    if (fabs(gtfinal - tfinal) > 1.0E-13) {
	throw CanteraError("shouldn't be here", "");
    }
    int iCsave = iCounter;
    for (int i = 0; i < nTimeStepsRegular_; i++) {
	tfinal = getNextRegularTime(tinit);
	if (tfinal <= tinit) {
	   throw CanteraError("shouldn't be here", "");
	}
	advanceTimeStepCounter();
	tinit = tfinal;
    }
    if (fabs(tfinal - gtfinal) > 1.0E-13) {
	throw CanteraError("shouldn't be here", "");
    }
    iCounter = iCsave;
    return iCounter;
}
//======================================================================================================================
double
SubIntegrationHistory::globalStartTime() const
{
    const TimeStepHistory& tshCurrent =  TimeStepList_[0];
    return tshCurrent.t_init_;
}
//======================================================================================================================
double
SubIntegrationHistory::globalEndTime() const
{
    const TimeStepHistory& tshCurrent =  TimeStepList_[nTimeSteps_-1];
    return tshCurrent.t_final_;
}
//======================================================================================================================
void
SubIntegrationHistory::print(int lvl) const
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
//======================================================================================================================

void SubIntegrationHistory::
setConstantStepSizeHistory(double gtinit, double gtfinal, int nsteps)
{
    zeroTimeStepCounter();
    double delta = (gtfinal - gtinit) / nsteps;
    for (int i = 0; i < nsteps; i++) {
	double tinit = gtinit + i * delta;
	double tfinal = gtinit + (i+1) * delta;
	addTimeStep(tinit, tfinal, tfinal, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
}
//======================================================================================================================
bool  SubIntegrationHistory::operator==(const SubIntegrationHistory& other) const
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
//======================================================================================================================
bool  SubIntegrationHistory::operator!=(const SubIntegrationHistory& other) const
{
    return ! (*this == other);
}
//======================================================================================================================
/*
 *  ELECTRODE_INPUT: constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.
 */
Electrode_Integrator::Electrode_Integrator() :
    Electrode(),
    neq_(0),
    deltaTsubcycleCalc_(0.0),
    atolResidNLS_(0),
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
    relativeLocalToGlobalTimeStepMinimum_(1.0E-3)
{

}
//======================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
Electrode_Integrator::Electrode_Integrator(const Electrode_Integrator& right) :
    Electrode(),
    neq_(0),
    deltaTsubcycleCalc_(0.0),
    atolResidNLS_(0),
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
    relativeLocalToGlobalTimeStepMinimum_(1.0E-3)
{
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
Electrode_Integrator&
Electrode_Integrator::operator=(const Electrode_Integrator& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    Electrode::operator=(right);
    neq_                                = right.neq_;
    deltaTsubcycleCalc_                 = right.deltaTsubcycleCalc_;

    rtolResidNLS_                       = right.rtolResidNLS_;
    atolResidNLS_                       = right.atolResidNLS_;
    atolNLS_                            = right.atolNLS_;
    rtolNLS_                            = right.rtolNLS_;
    ylowNLS_                            = right.ylowNLS_;
    yhighNLS_                           = right.yhighNLS_;
    yvalNLS_                            = right.yvalNLS_;
    ydotNLS_                            = right.ydotNLS_;
    errorLocalNLS_                      = right.errorLocalNLS_;
    errorGlobalNLS_                     = right.errorGlobalNLS_;

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
Electrode_Integrator::~Electrode_Integrator()
{
    SAFE_DELETE(jacPtr_);
    SAFE_DELETE(pSolve_);
}
//======================================================================================================================
int
Electrode_Integrator::electrode_model_create(ELECTRODE_KEY_INPUT* ei)
{

    Electrode::electrode_model_create(ei);


    setupIntegratedSourceTermErrorControl();

    return 0;
}
//======================================================================================================================
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

    create_solvers();

    return 0;
}
//===================================================================================================================
int
Electrode_Integrator::create_solvers()
{

    int neqNLS = nEquations();

    if (pSolve_) {
        int m = nEquations();
        if (m != (int) yvalNLS_.size()) {
            printf("shouldn't be here\n");
            exit(-1);
        }
        return nEquations();
    }



    yvalNLS_.resize(neqNLS, 0.0);
    ydotNLS_.resize(neqNLS, 0.0);
    ylowNLS_.resize(neqNLS, 0.0);
    yhighNLS_.resize(neqNLS, 0.0);
    errorLocalNLS_.resize(neqNLS, 0.0);
    errorGlobalNLS_.resize(neqNLS, 0.0);

    atolNLS_.resize(neqNLS, 1.0E-12);
    atolResidNLS_.resize(neqNLS, 1.0E-12);
    deltaBoundsMagnitudesNLS_.resize(neqNLS, 1.0E300);

    phaseJustDied_.resize(m_NumTotPhases, 0);
    phaseJustBorn_.resize(m_NumTotPhases, 0);

    // Add a couple of extra doubles, to the predictor, because some objects store extra info in those slots
    soln_predict_.resize(neqNLS + 2, 0.0);

    IntegratedSrc_Predicted.resize(numIntegratedSrc_, 0.0);
    IntegratedSrc_final_.resize(numIntegratedSrc_, 0.0);
    IntegratedSrc_Errors_local_.resize(numIntegratedSrc_, 0.0);
    IntegratedSrc_Errors_globalStep_.resize(numIntegratedSrc_, 0.0);
    atol_IntegratedSrc_global_.resize(numIntegratedSrc_, 0.0);


    jacPtr_ = new SquareMatrix(neqNLS, 0.0);

    pSolve_ = new NonlinearSolver(this);
    return neqNLS;
}
//===================================================================================================================
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
    for (int i = 0; i <  numIntegratedSrc_; i++) {
        atol_IntegratedSrc_global_[i] = 1.0E-14 * sm;
    }

    return numDofs;
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
void  Electrode_Integrator::resetStartingCondition(doublereal Tinitial, bool doResetAlways)
{
    /*
     * If the initial time is input, then the code doesn't advance
     */
    double tbase = MAX(t_init_init_, 1.0E-50);
    if (fabs(Tinitial - t_init_init_) < (1.0E-9 * tbase) && !doResetAlways) {
        return;
    }
    /*
     *  Zero the global error vectors
     */
    mdp::mdp_zero_dbl_1(&errorGlobalNLS_[0], neq_);

    Electrode::resetStartingCondition(Tinitial, doResetAlways);

}
//====================================================================================================================


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

    timeHistory_current_.clear();
   

    double pnorm = 0.0;
    int num_newt_its = 0;


    // Create the solvers and other vectors if we haven't already
    //  -> Note we haven't malloced the memory initially because we need to know the
    //     number of equations.
    if (pSolve_ == 0) {
        int neqNLS = create_solvers();
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
     *   Set up the global tolerances for the source term vector. Here we firm up what the source term vector
     *   is
     */
    setNLSGlobalSrcTermTolerances(GlobalRtolSrcTerm);

    /*
     *   Set the relative tolerance on the nonlinear solver
     */
    rtolNLS_      = GlobalRtolSrcTerm;
    rtolResidNLS_ = GlobalRtolSrcTerm * 1.0E-2;

    /*
     *   Set the absolute tolerance vector for the nonlinear solver, both the residual and the solution tolerances
     *   We control both because we have had trouble making sure that the equations are solved
     *   to a sufficient accuracy.
     */
    setResidAtolNLS();

    /*
     *   Zero vectors that are accumulated over local time step to represent global time step quantities
     *                   the integrated source term - This will be accumulated over the subcycling
     *                   Estimated global time step errors
     */
    zeroGlobalStepAccumulationTerms();

    /*
     *   Update the state
     */
    updateState();

    /*
     *  Determine the big mole fractions of multispecies phases
     */
    determineBigMoleFractions();

    /*
     *  Save the state into an XML state object
     */
    if (eState_final_) {
        eState_final_->copyElectrode_intoState(this);
        SAFE_DELETE(xmlStateData_init_);
        xmlStateData_init_ = eState_final_->writeStateToXML();
        startXML_TI_final();
    }
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
         *  Calculate the time interval from the previous step or from the initial value coming into the routine
         */
        deltaTsubcycle_ = deltaTsubcycleNext_;
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
         * Set the final time of the interval
         */
        tfinal_ = tinit_ + deltaTsubcycle_;
	tfinal_start = tinit_ + deltaTsubcycle_;

        /*
         *  Bound the time interval by the global time step
         */
        if (tfinal_ > t_init_init_ + deltaT) {
            tfinal_ = t_init_init_ + deltaT;
            deltaTsubcycle_ = tfinal_ - tinit_;
            nsteps_est =  iterSubCycle;
        }


        // HKM debugging point
#ifdef DEBUG_MODE_NOT
        if (fabs(tfinal_ - 0.02) < 0.001) {
            printf("we are here\n");
        }
#endif

        /*
         * Set the calced deltaT to the current deltaT.
         */
        deltaTsubcycleCalc_ = deltaTsubcycle_;

        /*
         *  Check internal consistency of the init solution
         *    -> failures here produce error exits
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
                throw CanteraError("Electrode_Integrator::integrate()",
                                   "FAILURE too many nonlinear convergence failures");
            }
	    tfinal_ = tinit_ + deltaTsubcycle_;
            /*
             *   Zero needed counters
             */
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
            int loglevelInput = MAX(0, printLvl_ - 2);

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
            jacPtr_->m_useReturnErrorCode = 1;

            /*
             *  Check the nonlinear residual equations for completeness and the ability to be solved
             */
            int retn = check_nonlinResidConditions();
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

            /*
             *  Set the tolerances on the solution variables
             */
            pSolve_->setAtol(DATA_PTR(atolNLS_));
            pSolve_->setRtol(0.1 * rtolNLS_);
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
            pSolve_->setResidualTols(rtolResidNLS_,  &atolResidNLS_[0]);

            pSolve_->setBoundsConstraints(&ylowNLS_[0], &yhighNLS_[0]);
            pSolve_->setDeltaBoundsMagnitudes(DATA_PTR(deltaBoundsMagnitudesNLS_));
            num_newt_its = 0;

#ifdef DEBUG_HKM_1DELECTRODE
            if (electrodeCellNumber_ == 0) {
                if (counterNumberIntegrations_ == 40) {
                    loglevelInput = 10;
                    detailedResidPrintFlag_ = 10;
                    enableExtraPrinting_ = 1;
                } else {
                    loglevelInput = 0;
                    detailedResidPrintFlag_ = 0;
                    enableExtraPrinting_ = 0;
                }
            }
#endif
            //   How to turn on the jacobian printing if debugging
            // pSolve_->s_print_NumJac = 1;
	    tfinal_start = tfinal_;

            int nonlinearFlag = pSolve_->solve_nonlinear_problem(solnType, &yvalNLS_[0], &ydotNLS_[0], 0.0,
								 tfinal_, *jacPtr_,  num_newt_its, num_linear_solves,
								 numBacktracks, loglevelInput);
            if (nonlinearFlag < 0) {
                if (printLvl_ > 2) {
                    printf("Electrode_Integrate::integrate(): Unsuccessful Nonlinear Solve flag = %d\n", nonlinearFlag);
                }
            }
	    bool specialSolve = false;
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
                        printf("deltaT reduced to %g from %g due to phase death capture\n",
                               deltaTsubcycleCalc_,  deltaTsubcycle_);
                    }
                    // Pick the next delta T to be equal to the current delta T
                    deltaTsubcycleNext_ = deltaTsubcycle_;
                    tfinal_ = tinit_ + deltaTsubcycleCalc_;
		    specialSolve = true;
                } else  if (deltaTsubcycleCalc_ > deltaTsubcycle_ * (1.0 + 1.0E-10)) {
                    if (printLvl_ > 1) {
                        printf("deltaT increased to %g from %g due to phase death capture\n",
                               deltaTsubcycleCalc_,  deltaTsubcycle_);
                    }
                    // Pick the next delta T to be equal to the current delta T
                    deltaTsubcycleNext_ = deltaTsubcycle_;
                    tfinal_ = tinit_ + deltaTsubcycleCalc_;
		    specialSolve = true;
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

        //
        // ACTIVITIES THAT ARE NEEDED TO DETERMINE IF THE SUBCYCLE IS ACCURATE ENOUGH
        // 
	double pnormSrc = 0.0;
	double pnormSoln = 0.0;
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
                stepAcceptable = decide_normHighLogic(pnorm);
            }

#ifdef DEBUG_MODE_PREDICTION
            predictorCorrectorPrint(yvalNLS_, pnormSrc, pnormSoln);
#else
            if (printLvl_ > 4) {
                predictorCorrectorPrint(yvalNLS_, pnormSrc, pnormSoln);
            }
#endif
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

       	if (subIntegrationType == FVDELTA_TIMEINTEGRATION_SIR || subIntegrationType == FIXEDSUBCYCLENUMBER_TIMEINTEGRATION_SIR) {
	    double tbase = timeHistory_base_.getNextRegularTime(tfinal_);
	    if (tfinal_ > tbase * 1.00000000001) {
		timeHistory_base_.advanceTimeStepCounter();
	    }
	    double tfinalnew = timeHistory_base_.getNextRegularTime(tfinal_);
	    double dtNew = tfinalnew - tfinal_;
	    
	    // we may compare the calculated deltaTnew with the base deltaTnew here if we want.
	    if (dtNew >  deltaTsubcycleNext_ * 1.001) {
		printf("Warning\n");
	    }
	    // Right now, let's just assign it.
	    deltaTsubcycleNext_ = dtNew;
	}


        /*
         *  If this is the first subcycle time step, calcualte the first time step for the next global iteration
         */
        if ((choiceDeltaTsubcycle_init_ == 1) && (iterSubCycle == 1)) {
            deltaTsubcycle_init_next_ = MIN(deltaTsubcycle_init_next_, deltaTsubcycleNext_);
        }

        /*
         * Determine if we are at the end of the global time step
         */
        if ((1.0000000000001) * tfinal_ >= t_init_init_ + deltaT) {
            notDone = false;
            tfinal_ = t_init_init_ + deltaT;
        }

        // HKM debugging point
#ifdef DEBUG_MODE_NOT
        if (fabs(tfinal_ - 0.02) < 0.001) {
            printf("we are here\n");
        }
#endif

        /*
         * Potential write out intermediate data to a CSV file
         */
        if (printCSVLvl_ > 0) {
            if (notDone) {
                writeCSVData(1);
            } else {
                writeCSVData(2);
            }
        }

        /*
         *  Check the final state for any error conditions
         */
        check_final_state();

        /*
         *  Calcualte the first time step for the next global iteration
         */
        if (choiceDeltaTsubcycle_init_ == 0) {
            deltaTsubcycle_init_next_ = deltaTsubcycleNext_;
        }

        /*
         *  Save the state of the object
         */
        if (eState_final_) {
            eState_final_->copyElectrode_intoState(this);
            SAFE_DELETE(xmlStateData_final_);
            xmlStateData_final_ = eState_final_->writeStateToXML();
            makeXML_TI_intermediate();
            addtoXML_TI_final(notDone);
        }

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
    throw CanteraError("Electrode_Integrator::setNLSGlobalSrcTermTolerances","unimplemented");
}

//======================================================================================================================
void  Electrode_Integrator::zeroGlobalStepAccumulationTerms()
{
    /*
     *  Zero the global error vectors
     */
    mdp::mdp_zero_dbl_1(&errorGlobalNLS_[0], neq_);
    /*
     *   Zero the integrated source term - This will be accumulated over the subcycling
     */
    mdp::mdp_zero_dbl_1(&spMoleIntegratedSourceTerm_[0], m_NumTotSpecies);
}
//==================================================================================================================
// Number of degrees of freedom in the integrated source terms that constitute the output
int Electrode_Integrator::numIntegratedSrcTerms() const
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
    throw CanteraError("Electrode_Integrator::setResidAtolNLS()", "unimplemented");
}
//==================================================================================================================
//  Predict the solution
/*
 *  (virtual from Electrode_Integrator)
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
int  Electrode_Integrator::predictSoln()
{
    throw CanteraError(" Electrode_Integrator::predictSoln()","unimplemented");
}
//==================================================================================================================
// Extract information from cantera
/*
  *  (virtual fucntion from Electrode_Integrator)
 */
void  Electrode_Integrator::extractInfo()
{
    throw CanteraError("Electrode_Integrator::extractInfo()",  "unimplemented");
}
//==================================================================================================================
// Collect mole change information
/*
 *  (virtual fucntion from Electrode_Integrator)
 */
void Electrode_Integrator::updateSpeciesMoleChangeFinal()
{
    throw CanteraError(" Electrode_Integrator::updateSpeciesMoleChangeFinal()",  "unimplemented");

}
//==================================================================================================================
// calculate the residual
/*
 *
 *  (virtual fucntion from Electrode_Integrator)
 *
 */
int Electrode_Integrator::calcResid(doublereal* const resid, const ResidEval_Type_Enum evalType)
{
    throw CanteraError(" Electrode_Integrator::calcResid()",  "unimplemented");
    return 0;
}
//==================================================================================================================
//  Gather the predicted solution values and the predicted integrated source terms
/*
 *  (virtual from Electrode_Integrator)
 *
 *  Both the predicted solution values and the predicted integrated source terms are used
 *  in the time step control
 */
void Electrode_Integrator::gatherIntegratedSrcPrediction()
{
    throw CanteraError(" Electrode_Integrator::gatherIntegratedSrcPrediction()",  "unimplemented");
}

//==================================================================================================================
//   Calculate the integrated source terms and do other items now that we have a completed time step
/*
 *  Calculate source terms on completion of a step. At this point we have solved the nonlinear problem
 *  for the current step, and we are calculating post-processed quantities like source terms.
 */
void Electrode_Integrator::calcSrcTermsOnCompletedStep()
{
    /*
     *  Calculate the integrated source term
     *       An alternative would be to redo the residual calculation. However, here we assume that
     *       the residual calculation has been done and the results are in _final_
     */
    for (int i = 0; i < m_NumTotSpecies; i++) {
        spMoleIntegratedSourceTermLast_[i] = spMoles_final_[i] - spMoles_init_[i];

    }
}
//==================================================================================================================
//   Accumulate src terms and other results from the local step into the global holding bins.
/*
 *  Accumulate source terms on completion of a step. At this point we have solved the nonlinear problem
 *  for the current step and we have satisfied all accuracy requirements.
 *  The step is good. We now accumulate the results before going on to a new local step.
 */
void Electrode_Integrator::accumulateSrcTermsOnCompletedStep(bool remove)
{
    if (remove) {
        for (int i = 0; i < m_NumTotSpecies; i++) {
            spMoleIntegratedSourceTerm_[i] -= spMoleIntegratedSourceTermLast_[i];
        }
    } else {
        for (int i = 0; i < m_NumTotSpecies; i++) {
            spMoleIntegratedSourceTerm_[i] += spMoleIntegratedSourceTermLast_[i];
        }
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
void Electrode_Integrator::setFinalStateFromInit_Oin()
{
    Electrode::setFinalStateFromInit_Oin();
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
    Electrode_Integrator::setFinalStateFromInit_Oin();
}
//==================================================================================================================
// Return the number of equations in the equation system that is used to solve the ODE integration
int Electrode_Integrator::nEquations() const
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
int Electrode_Integrator::calcDeltaSolnVariables(const doublereal t, const doublereal* const ySoln,
        const doublereal* const ySolnDot, doublereal* const deltaYSoln,
        const doublereal* const solnWeights)
{
    int retn = ResidJacEval::calcDeltaSolnVariables(t, ySoln, ySolnDot,deltaYSoln,solnWeights);
    return retn;
}
//==================================================================================================================
//   Unpack the soln vector
/*
 *  (virtual from Electrode_Integrator)
 *
 *  This function unpacks the solution vector into  phaseMoles_final_,  spMoles_final_, and spMf_final_[]
 */
void Electrode_Integrator::unpackNonlinSolnVector(const double* const y)
{
    throw CanteraError(" Electrode_Integrator::unpackNonlinSolnVector()", "unimplemented");
}
//====================================================================================================================
// Check to see that the preceding step is a successful one
/*
 *   We check to see if the preceding step is a successful one.
 *
 *  @return Returns a bool true if the step is acceptable, and false if it is unacceptable.
 */
bool Electrode_Integrator::checkSubIntegrationStepAcceptable() const
{
    return true;
}
//====================================================================================================================
//  Calculate the norm of the difference between the predicted answer and the final converged answer
//  for the current time step
/*
 *   The norm calculated by this routine is used to determine whether the time step is accurate enough.
 *
 *  @return    Returns the norm of the difference. Normally this is the L2 norm of the difference
 */
double Electrode_Integrator::predictorCorrectorWeightedSolnNorm(const std::vector<double>& yvalNLS)
{
    double pnorm = l0normM(soln_predict_, yvalNLS, neq_, atolNLS_, rtolNLS_);
    return pnorm;
}
//====================================================================================================================
// On a completed local step, accumulate local errors into the global error vectors for the global time step
/*
 */
void Electrode_Integrator::accumulateLocalErrorToGlobalErrorOnCompletedStep()
{
    for (int k = 0; k < neq_; k++) {
        errorGlobalNLS_[k] += errorLocalNLS_[k];
    }
}
//====================================================================================================================
// When the norm is high this routine decides what to do with the step
/*
 *  @param  pnorm  Norm of the step
 *
 *  @return   Returns whether the step is accepted or not
 */
bool Electrode_Integrator::decide_normHighLogic(double pnorm)
{
    double deltaT = t_final_final_ - t_init_init_;
    if (pnorm > 2.0) {
        if (deltaTsubcycle_ > relativeLocalToGlobalTimeStepMinimum_ * deltaT) {
            if (choiceDeltaTsubcycle_init_ != 2) {
                if (printLvl_ > 4) {
                    printf("\t\tdecide_normHighLogic(): WARNING Step rejected because norm is greater than 2\n");
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
//====================================================================================================================
// Print table representing prediction vs. corrector information
/*
 *  @param yval           Vector of corrector values
 *  @param pnormSrc       Norm of the predictor-corrector comparison for the source vector.
 *  @param pnormSoln      Norm of the predictor-corrector comparison for the solution vector.
 */
void  Electrode_Integrator::predictorCorrectorPrint(const std::vector<double>& yval,
        double pnormSrc, double pnormSoln)
{
    double atolVal = 1.0E-8;
    double denom;
    double tmp;
    printf(" -------------------------------------------------------------------------------------------------------------------\n");
    printf(" PREDICTOR_CORRECTOR  SubIntegrationCounter = %7d       t_init = %12.5E,       t_final = %12.5E       deltaT = %12.5E\n",
           counterNumberSubIntegrations_, tinit_, tfinal_, tfinal_ - tinit_);
    printf("                         IntegrationCounter = %7d  t_init_init = %12.5E, t_final_final = %12.5E   deltaTGlob = %12.5E\n",
           counterNumberIntegrations_, t_init_init_, t_final_final_, t_final_final_ - t_init_init_);
    printf(" -------------------------------------------------------------------------------------------------------------------\n");
    printf("                          Initial      Prediction      Actual          Difference         Tol   Contrib      |\n");

    denom = MAX(fabs(yval[0]), fabs(soln_predict_[0]));
    denom = MAX(denom, atolVal);
    tmp = fabs((yval[0] - soln_predict_[0])/ denom);
    printf("DeltaT                   | %14.7E %14.7E %14.7E | %14.7E | %10.3E | %10.3E |\n",
           deltaTsubcycleCalc_, soln_predict_[0],  yval[0], yval[0] - soln_predict_[0], atolVal, tmp);
    for (int i = 1; i < (int) yval.size(); i++) {
        denom = MAX(fabs(yval[i]), fabs(soln_predict_[i]));
        denom = MAX(denom, atolVal);
        tmp = fabs((yval[i] - soln_predict_[i])/ denom);
        printf(" soln %3d |              %14.7E %14.7E | %14.7E | %10.3E | %10.3E | \n",
               i, soln_predict_[i],  yval[i], yval[i] - soln_predict_[i], atolVal, tmp);
    }
    printf(" -------------------------------------------------------------------------------------------------------------------\n");
    printf("                                                                                                        %10.3E\n",
           pnormSoln);
}
//====================================================================================================================
// Possibly change the solution due to phase births and deaths.
/*
 *   (virtual from Electrode_Integrator)
 *
 */
bool Electrode_Integrator::changeSolnForBirthDeaths()
{
    return true;
}
//====================================================================================================================
// Possibly change the solution due to phase births and deaths.
/*
 *   (virtual from Electrode_Integrator)
 *
 *  This routine is carried out after the step is deemed a success. Massaging of the solution
 *  must be carried out within strict tolerances.
 */
void Electrode_Integrator::manageBirthDeathSuccessfulStep()
{
}
//====================================================================================================================
//   Error check on the routine step
/*
 *    (virtual from Electrode_Integrator)
 *
 *   Error checks go here. All errors are fatal exits.
 */
void Electrode_Integrator::check_final_state()
{
}
//====================================================================================================================
// Print a header for the residual printout
/*
 *  (virtual from Eelctrode_Integrator)
 */
void Electrode_Integrator::printResid_TimeHeader()
{
}
//====================================================================================================================
//   Check for problems with the residual
/*
 *  (virtual from Electrode_Integrator)
 *
 *  Checks here will cause the current nonlinear solve to fail
 *
 *    @return Return a negative value if there is a fatal problem.
 */
int  Electrode_Integrator::residEval_BaseChecks()
{
    return 1;
}
//====================================================================================================================
//     Print details about the satisfaction of the residual
/*
 *  (virtual from Electrode_Integrator)
 */
void  Electrode_Integrator::printResid_ResidSatisfaction()
{
}
//====================================================================================================================
//  Pack the nonlinear solver proplem
/*
 *  formulate the nonlinear solver problem to be solved.
 *     Fields to be filled in
 *             yvalNLS_
 *             ylowNLS_
 *             yhighNLS_
 *             atolNLS_
 *             deltaBoundsMagnitudesNLS_
 *             ysType
 */
void Electrode_Integrator::initialPackSolver_nonlinFunction()
{
    throw CanteraError(" Electrode_Integrator::initialPackSolver_nonlinFunction()",  "unimplemented");
}
//====================================================================================================================
// Check the nonlinear residual equations for completeness and the ability to be solved
/*
 *
 */
int Electrode_Integrator::check_nonlinResidConditions()
{
    return 0;
}
//====================================================================================================================
// Report the number of state variables and their relative integration errors during the
// current global integration step
/*
 *   Note rtol doesn't factor into this immediately. Therefore, a value or 1E-3
 *                                  would mean the error in the value is 1 part in 1000.
 *
 *  @param[out] numSV               Returns the number of state variables
 *  @param[out] errorVector         Returns a vector of errors in the state variables for the global step
 *                                  Note rtol doesn't factor into this immediately. Therefore, a value or 1E-3
 *                                  would mean the error in the value is 1 part in 1000.
 *  @return     Returns the large value of the errors in the errorVector.
 */
double Electrode_Integrator::reportStateVariableIntegrationError(int& numSV, double* const errorVector) const
{
    numSV = neq_;
    double vmax = 0.0;
    if (errorVector) {
        for (size_t k = 0; k < (size_t) neq_; k++) {
            errorVector[k] = errorGlobalNLS_[k];
            vmax = MAX(vmax, errorVector[k]);
        }
    } else {
        for (size_t k = 0; k < (size_t) neq_; k++) {
            vmax = MAX(vmax, errorGlobalNLS_[k]);
        }
    }
    return vmax;
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
int Electrode_Integrator::evalResidNJ(const doublereal t, const doublereal delta_t,
                                      const doublereal* const y, const doublereal* const ySolnDot,
                                      doublereal* const resid,
                                      const ResidEval_Type_Enum evalType, const int id_x,
                                      const doublereal delta_x)
{
    int retn = 1;
    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
        printf("\t\t===============================================================================================================================\n");
        printf("\t\t  EXTRA PRINTING FROM NONLINEAR RESIDUAL: ");
        if (evalType ==  Base_ResidEval) {
            printf(" BASE RESIDUAL");
        } else if (evalType == JacBase_ResidEval) {
            printf(" BASE JAC RESIDUAL");
        } else  if (evalType == JacDelta_ResidEval) {
            printf(" DELTA JAC RESIDUAL");
            printf(" var = %d delta_x = %12.4e Y_del = %12.4e Y_base = %12.4e", id_x, delta_x, y[id_x], y[id_x] - delta_x);
        } else  if (evalType == Base_ShowSolution) {
            printf(" BASE RESIDUAL - SHOW SOLUTION");
        }
        printf(" DomainNumber = %2d , CellNumber = %2d , IntegrationCounter = %d\n",
               electrodeDomainNumber_, electrodeCellNumber_, counterNumberIntegrations_);
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


    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
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
int  Electrode_Integrator::getInitialConditions(const doublereal t0, doublereal* const y, doublereal* const ydot)
{
    return 0;
}
//====================================================================================================================
// Print conditions of the electrode for the current integration step to stdout
/*
 *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
 *                       to the final_final time.
 *                       The default is to print out the source terms
 *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
 *                       time step. The default is to print out the global values
 */
void Electrode_Integrator::printElectrode(int pSrc, bool subTimeStep)
{
    int iph;
    double egv = TotalVol();
    double tsm = SolidTotalMoles();
    printf("   ==============================================================================================\n");
    if (subTimeStep) {
        printf("      Electrode at intermediate-step time final = %g\n", tfinal_);
        printf("                   intermediate-step time init  = %g                deltaT = %g\n",
               tinit_, deltaTsubcycleCalc_);
        printf("                    Model Type = %3d , DomainNumber = %2d , CellNumber = %2d , SubIntegrationCounter = %d\n",
               electrodeModelType_, electrodeDomainNumber_, electrodeCellNumber_, counterNumberSubIntegrations_);
    } else {
        printf("      Electrode at time final = %g\n", t_final_final_);
        printf("                   time init  = %g                         deltaTglobal = %g\n", t_init_init_,
               t_final_final_ - t_init_init_);
        printf("                    Model Type = %3d , DomainNumber = %2d , CellNumber = %2d , IntegrationCounter = %d\n",
               electrodeModelType_, electrodeDomainNumber_, electrodeCellNumber_, counterNumberIntegrations_);
    }
    printf("   ==============================================================================================\n");
    printf("\n");
    printf("          Number of external surfaces = %d\n", numExternalInterfacialSurfaces_);
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
    if (numSurfaces_ > m_NumSurPhases) {
        m = numSurfaces_ + NumVolPhases_;
    }
    for (iph = 0; iph < m; iph++) {
        printElectrodePhase(iph, pSrc, subTimeStep);
    }

}
//===================================================================================================================

void Electrode_Integrator::printElectrodePhase(int iph, int pSrc, bool subTimeStep)
{
    int isph;
    double* netROP = new double[m_NumTotSpecies];
    ThermoPhase& tp = thermo(iph);
    string pname = tp.id();
    int istart = m_PhaseSpeciesStartIndex[iph];
    int nsp = tp.nSpecies();
    printf("     ===============================================================\n");
    printf("          Phase %d %s \n", iph,pname.c_str());
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
    if (iph >= NumVolPhases_) {
        isph = iph - NumVolPhases_;
        printf("                surface area (final) = %g\n",  surfaceAreaRS_final_[isph]);
        printf("                surface area (init)  = %g\n",  surfaceAreaRS_init_[isph]);
        int ddd =  isExternalSurface_[isph];
        printf("                IsExternalSurface = %d\n", ddd);
        double oc = openCircuitVoltage(isph);
        if (oc != 0.0) {
            printf("                 Open Circuit Voltage = %g\n", oc);
        }
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
        const vector<double>& rsSpeciesProductionRates = RSD_List_[isph]->calcNetProductionRates();
        RSD_List_[isph]->getNetRatesOfProgress(netROP);

        doublereal* spNetProdPerArea = (doublereal*) spNetProdPerArea_List_.ptrColumn(isph);
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

//====================================================================================================================
double Electrode_Integrator::l0normM(const std::vector<double>& v1, const std::vector<double>& v2, int num,
                                     const std::vector<double>& atolVec, const double rtol)
{
    double max0 = 0.0;
    double denom, diff, ee;

    for (int k = 0; k < num; k++) {

        diff = fabs(v1[k] - v2[k]);
        denom = rtol * MAX(fabs(v1[k]), fabs(v2[k]));
        denom = MAX(denom, atolVec[k]);
        ee = diff / denom;
        if (ee > max0) {
            max0 = ee;
        }
        errorLocalNLS_[k] = ee * rtol;
    }
    return max0;
}

//====================================================================================================================
int Electrode_Integrator::setTimeHistoryBaseFromCurrent()
{
    timeHistory_base_ = timeHistory_current_;
    return timeHistory_base_.nTimeStepsRegular_;
}
//====================================================================================================================
int Electrode_Integrator::setTimeHistoryBase(const SubIntegrationHistory &timeHistory)
{
    timeHistory_base_ = timeHistory;
    return timeHistory.nTimeStepsRegular_;
}
//====================================================================================================================
SubIntegrationHistory& Electrode_Integrator::timeHistory(bool returnBase)
{
    if (returnBase) {
	return timeHistory_base_;
    }
    return timeHistory_current_;
}

} // End of namespace Cantera
//======================================================================================================================
