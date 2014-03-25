/**
 *  @file BEulerInt.cpp
 *
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "BEulerInt.h"

#include "cantera/numerics/SquareMatrix.h"

#include "mdp_allo.h"

#include "m1d_SolNonlinear.h"
#include "m1d_ProblemResidEval.h"
#include "m1d_EpetraJac.h"
#include "m1d_Comm.h"

#include "Epetra_Vector.h"
#include "m1d_ProblemStatement.h"

#include <iostream>

#ifdef DEBUG_HKM
#include <cstdio>
extern int iDebug_HKM;
extern int iTimeStep_HKM;
#endif

#include <fstream>

using namespace std;
using namespace Cantera;
using namespace m1d;

namespace beuler {

static void print_line(const char *str, int n)
{
    for (int i = 0; i < n; i++) {
        printf("%s", str);
    }
    printf("\n");
}
//=====================================================================================================================
/*
 * Exception thrown when a BEuler error is encountered. We just call the
 * Cantera Error handler in the initialization list
 */
BEulerErr::BEulerErr(std::string msg) :
        CanteraError("BEulerInt", msg)
{
}
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
/*
 *  Constructor. Default settings: dense jacobian, no user-supplied
 *  Jacobian function, Newton iteration.
 */
BEulerInt::BEulerInt() :
        m_nonlin(0),
        Comm_ptr_(0),
        mypid_(0),
        m_method(BEulerVarStep),
        m_doSpecialStartCalc(0),
        m_jacFormMethod(BEULER_JAC_NUM),
        m_rowScaling(true),
        m_colScaling(true),
        m_matrixConditioning(false),
        m_reltol(1.e-4),
        m_abstol(0),
        m_ewt(0),
        m_order(1),
        m_maxord(1),
        delta_t_max(1.0E300),
        delta_t_min(1.0E-300),
        m_min_newt_its(0),
        m_time_step_num(0),
        m_time_step_attempts(0),
        m_max_time_step_attempts(11000000),
        m_dumpJacobians(false),
        doResidSolnDamping_(true),
        doDeltaDamping_(true),
        m_numInitialConstantDeltaTSteps(0),
        m_failure_counter(0),
        m_print_flag(3),
        m_printSolnStepInterval(1),
        m_printSolnNumberToTout(1),
        m_printSolnFirstSteps(0),
        m_y_n(0),
        m_y_nm1(0),
        m_y_pred_n(0),
        m_ydot_n(0),
        m_ydot_nm1(0),
        m_y_n_owned(0),
        m_ydot_n_owned(0),
        m_resid(0),
        m_isAlgebraic(0),
        m_t0(0.0),
        m_time_final(0.0),
        time_n(0.0),
        time_nm1(0.0),
        time_nm2(0.0),
        delta_t_n(0.0),
        delta_t_nm1(0.0),
        delta_t_nm2(0.0),
        delta_t_np1(1.0E-8),
        m_residWts(0),
        m_wksp(0),
        m_func(0),
        m_rowScales(0),
        m_colScales(0),
        tdjac_ptr(0),
        m_nfe(0),
        m_nJacEval(0),
        m_numTotalNewtIts(0),
        m_numTotalLinearSolves(0),
        m_numTotalConvFails(0),
        m_numTotalTruncFails(0),
        num_failures(0),
        m_DAEDOT_Predicted_Bad(1),
        m_DOT_Predicted_Bad(0),
        m_DAE_Predicted_Bad(0),
        m_isArithmeticScaled(0),
        m_numTimeRegions(1),
        m_currentTimeRegion(0),
        m_timeRegionBoundaries(0)
{
}
//=====================================================================================================================
// Copy constructor
/*
 *
 * @param r  Object to be copied
 */
BEulerInt::BEulerInt(const BEulerInt &r)
{
    *this = r;
}
//=====================================================================================================================
BEulerInt &
BEulerInt::operator=(const BEulerInt &r)
{
    if (this == &r) {
        return *this;
    }
    Cantera::Integrator::operator=(r);

    return *this;
}
//=====================================================================================================================
// Destructor
BEulerInt::~BEulerInt()
{
    safeDelete(m_nonlin);
    safeDelete(tdjac_ptr);

    safeDelete(m_y_n_owned);
    safeDelete(m_ydot_n_owned);
    safeDelete(m_y_n);
    safeDelete(m_y_nm1);
    safeDelete(m_y_pred_n);
    safeDelete(m_ydot_n);
    safeDelete(m_ydot_nm1);
    safeDelete(m_resid);
    safeDelete(m_isAlgebraic);
    safeDelete(m_residWts);
    safeDelete(m_wksp);
    safeDelete(m_ewt);
    safeDelete(m_abstol);
    safeDelete(m_rowScales);
    safeDelete(m_colScales);
    safeDelete(m_isArithmeticScaled);
}
//=====================================================================================================================
//  Function called by BEuler to evaluate the Jacobian matrix and the
//  current residual at the current time step.
/*
 *  @param J = Jacobian matrix to be filled in
 *  @param res This routine returns the current value of the residiaul
 *  @param time_curr current time
 *  @param CJ   coefficient for the time derivative term
 *  @param y_curr    current solution
 *  @param y_curr_dot  Current solution derivative
 *  @param num_newt_its  number of newton iterations
 *  @param jacType        0  Normal jacobian
 *                        1  DAE initialization problem
 */
void BEulerInt::beuler_jac(m1d::EpetraJac &tjac, m1d::Epetra_Vector_Owned * const res, double time_curr, double CJ,
                           m1d::Epetra_Vector_Ghosted * const y_curr, m1d::Epetra_Vector_Ghosted * const y_curr_dot,
                           int num_newt_its, int jacType)
{
    m_nfe++;
    m_nJacEval++;
    Solve_Type_Enum solnType = TimeDependentAccurate_Solve;
    if (jacType == 1) {
        solnType = DAESystemInitial_Solve;
    }

    tjac.matrixResEval(true, y_curr, y_curr_dot, res, time_curr, CJ, solnType);
}

//=====================================================================================================================
void BEulerInt::setTolerancesEpetra(double reltol, const m1d::Epetra_Vector_Owned &abstol)
{
    for (int i = 0; i < m_NumLcOwnedEqns; i++) {
        (*m_abstol)[i] = abstol[i];
    }
    m_reltol = reltol;
}
//=====================================================================================================================
void BEulerInt::setAlgebraicConstraints(const Epetra_IntVector &algebraicConstraint)
{
    for (int i = 0; i < m_NumLcOwnedEqns; i++) {
        (*m_isAlgebraic)[i] = algebraicConstraint[i];
    }
}
//=====================================================================================================================
void BEulerInt::setProblemType(int jacFormMethod)
{
    m_jacFormMethod = jacFormMethod;
}
//=====================================================================================================================
void BEulerInt::setTimeRegionBoundaries(const std::vector<double>& timeRegionBoundaries)
{
    m_timeRegionBoundaries = timeRegionBoundaries;
    m_numTimeRegions = m_timeRegionBoundaries.size() - 1;
    m_currentTimeRegion = 0;
    m_func->m_numTimeRegions = m_numTimeRegions;
    m_func->m_currentTimeRegion = m_currentTimeRegion;

}
//=====================================================================================================================
void BEulerInt::setTimeRegion(double time, int iregion)
{
    if (m_timeRegionBoundaries.size() > 0) {
        if (m_timeRegionBoundaries[iregion] > time || m_timeRegionBoundaries[iregion + 1] < time) {
            throw m1d_Error("BEulerInt::setTimeRegion()",
                    "time, region combination outside of previously specified time boundaries");
        }
    }
    m_currentTimeRegion = iregion;
    m_func->m_currentTimeRegion = m_currentTimeRegion;
}
//=====================================================================================================================
//  Specify the nonlinear solver to be used in the problem
/*
 *   Note, a default nonlinear solver is used, so this step is not necessary necessarily
 *   Note, the BEulerInt object owns the nonlinear solver, so nothing needs to be freed.
 */
void BEulerInt::specifyNonLinearSolver(m1d::SolGlobalNonlinear *nonlin)
{
    if (m_nonlin) {
        safeDelete(m_nonlin);
    }
    if (!nonlin) {
        m_nonlin = new SolNonlinear();
    }

    m_nonlin = nonlin;

    // Setup the problem for solution.
    if (tdjac_ptr) {
        Solve_Type_Enum stype = TimeDependentAccurate_Solve;
        m_nonlin->setup_problem(stype, m_y_n, m_ydot_n, 0.0, *m_func, *tdjac_ptr);
    }

}
//=====================================================================================================================
void BEulerInt::setMethodBEMT(BEulerMethodType t)
{
    m_method = t;
}
//=====================================================================================================================
void BEulerInt::setMaxStep(doublereal hmax)
{
    delta_t_max = hmax;
}
//=====================================================================================================================
void BEulerInt::setMinStep(doublereal hmin)
{
    delta_t_min = hmin;
}
//=====================================================================================================================
void BEulerInt::setMaxNumTimeSteps(int maxNumTimeSteps)
{
    m_max_time_step_attempts = maxNumTimeSteps;
}
//=====================================================================================================================
void BEulerInt::setNumInitialConstantDeltaTSteps(int num)
{
    m_numInitialConstantDeltaTSteps = num;
}
//=====================================================================================================================
/*
 * setPrintSolnOptins():
 *
 * This routine controls when the solution is printed
 *
 * @param printStepInterval If greater than 0, then the
 *                     soln is printed every printStepInterval
 *                     steps.
 *
 * @param printNumberToTout The solution is printed at
 *                  regular invervals a total of
 *                  "printNumberToTout" times.
 *
 * @param printSolnFirstSteps The solution is printed out
 *                   the first "printSolnFirstSteps"
 *                   steps. After these steps the other
 *                   parameters determine the printing.
 *                   default = 0
 *
 * @param dumpJacobians Dump jacobians to disk.
 *
 *                   default = false
 *
 */
void BEulerInt::setPrintSolnOptions(int printSolnStepInterval, int printSolnNumberToTout, int printSolnFirstSteps,
                                    bool dumpJacobians)
{
    m_printSolnStepInterval = printSolnStepInterval;
    m_printSolnNumberToTout = printSolnNumberToTout;
    m_printSolnFirstSteps = printSolnFirstSteps;
    m_dumpJacobians = dumpJacobians;
}
//=====================================================================================================================
/*
 *
 * setNonLinOptions()
 *
 *  Set the options for the nonlinear method
 *
 *  Defaults are set in the .h file. These are the defaults:
 *    min_newt_its = 0
 *    matrixConditioning = false
 *    colScaling = false
 *    rowScaling = true
 */
//=====================================================================================================================
void BEulerInt::setNonLinOptions(int min_newt_its, bool matrixConditioning, bool colScaling, bool rowScaling)
{
    m_min_newt_its = min_newt_its;
    m_matrixConditioning = matrixConditioning;
    m_colScaling = colScaling;
    m_rowScaling = rowScaling;
    if (m_colScaling) {
        if (!m_colScales) {
            m_colScales = new Epetra_Vector(m_resid->Map());
            m_colScales->PutScalar(1.0);
        }
    }
    if (m_rowScaling) {
        if (!m_rowScales) {
            m_rowScales = new Epetra_Vector(m_resid->Map());
            m_rowScales->PutScalar(1.0);
        }
    }
    /*
     *  Pass parameters down to the nonlinear solver
     */
    m_nonlin->setNonLinOptions(min_newt_its, matrixConditioning, colScaling, rowScaling);
}
//=====================================================================================================================
/*
 * setInitialTimeStep():
 *
 * Set the initial time step. Right now, we set the
 * time step by setting delta_t_np1.
 */
void BEulerInt::setInitialTimeStep(double deltaT)
{
    delta_t_n = deltaT;
    delta_t_np1 = deltaT;
    if (delta_t_np1 < 1.0E-300) {
        throw BEulerErr("BEulerInt::setInitialTimeStep(double deltaT) ERROR: input below zero");
    }
}
//=====================================================================================================================
// Set the level of printing that occurs for each time step
/*
 *   0 -> absolutely nothing is printed for a single time step.
 *   1 -> One line summary per time step
 *   2 -> short description, points of interest
 *   3 -> More printed per time step -> major algorithm issues are displayed
 *   4 -> Additional time step error control information is printed out
 *   5 -> Summaries of the nonlinear solve iterates are printed out
 *   6 -> Algorithm information on the nonlinear iterates are printed out
 *   7 -> Additional info on the nonlinear iterates are printed out
 *   8 -> Additional info on the linear solve is printed out.
 *   9 -> Info on a per iterate of the linear solve is printed out.
 */
void BEulerInt::setPrintFlag(int print_flag, int nonlinear_print_flag)
{
    m_print_flag = print_flag;
    if (nonlinear_print_flag < 0) {
        if (m_print_flag <= 3) {
            m_nonlin->setPrintFlag(0);
        } else if (m_print_flag == 4) {
            m_nonlin->setPrintFlag(1);
        } else if (m_print_flag == 5) {
            m_nonlin->setPrintFlag(3);
        } else if (m_print_flag == 6) {
            m_nonlin->setPrintFlag(5);
        } else if (m_print_flag == 7) {
            m_nonlin->setPrintFlag(6);
        } else if (m_print_flag == 8) {
            m_nonlin->setPrintFlag(7);
        } else if (m_print_flag == 9) {
            m_nonlin->setPrintFlag(8);
        } else if (m_print_flag > 9) {
            m_nonlin->setPrintFlag(m_print_flag);
        }
    } else {
        m_nonlin->setPrintFlag(nonlinear_print_flag);
    }
}
//=====================================================================================================================
//  Setup the problem by specifying the object that generates the problem to be solved,
//  mallocing internal objects and memory and by calling the problem residual to establish
//  the initial conditions.
/*
 * Find the initial conditions for y and ydot.
 *
 * @param t0     Initial time for the simulation
 * @param func   func is the object that supplies the problem.
 */
void BEulerInt::initializePRE(m1d::ProblemResidEval &func )
{
    m_NumGlEqns = func.m_neq;
    m_NumLcEqns = func.m_NumLcEqns;
    m_NumLcOwnedEqns = func.m_NumLcOwnedEqns;
    //m_t0 = t0;
    m_ghostedMap = func.ghostedMap();
    m_ownedMap = func.ownedMap();

    m_func = &func;

    internalMalloc();

    Comm_ptr_ = & (m_ghostedMap->Comm());
    mypid_ = Comm_ptr_->MyPID();

    /*
     * Get the initial conditions.
     */
    //  m_func->initialConditions(true, m_y_n, m_ydot_n, m_t0, delta_t_n);
    ProblemStatement *ps = m_func->psInput_ptr_;
    m_func->createMatrix(ps->I_LinearSolverBlock);

    /*
     *  Get a pointer to the time dependent jacobian for the problem
     */
    tdjac_ptr = &func.jacobian();
    /*
     * Initialize the nonlinear solver, if not done already
     */
    if (!m_nonlin) {
        m_nonlin = new SolNonlinear();
    }
    // Setup the problem for solution.
    Solve_Type_Enum stype = TimeDependentAccurate_Solve;
    m_nonlin->setup_problem(stype, m_y_n, m_ydot_n, 0.0, func, *tdjac_ptr);

    /*
     * Initialize the various time counters in the object
     */
    time_n = m_t0;
    time_nm1 = time_n;
    time_nm2 = time_nm1;
    delta_t_nm1 = 0.0;

    /*
     *  zero out the solution file
     */
    string sname = m_func->getBaseFileName() + ".xml";
    fstream ff(sname.c_str(), fstream::in | fstream::out | fstream::trunc);
    ff.close();

    /*
     *
     */
    m_func->fillIsAlgebraic(*m_isAlgebraic);
}
//=====================================================================================================================
// We input t0 here. However, the initialConditions function may override.
void BEulerInt::determineInitialConditions(double t0, double delta_t)
{
    m_t0 = t0;
    /*
     *  Get the initial conditions. The initial conditions will depend on delta_t_n, as source term may be integrals
     *  over time. This means that their values depend on the time increment.
     *  We may get a delta_t_np1 from this procedure.
     */
    delta_t_n = delta_t;
    m_func->initialConditions(true, m_y_n, m_ydot_n, m_t0, delta_t_n, delta_t_np1);
    //
    //   Set up the initial time step
    //
    ProblemStatement *ps = m_func->psInput_ptr_;
    if (ps->initialTimeStep_ > 0.0) {
       setInitialTimeStep(ps->initialTimeStep_);
    } else {
       setInitialTimeStep(delta_t_np1);
    }
   
    //
    // Setup the inital delta T constant field
    //
    m_numInitialConstantDeltaTSteps = ps->m_numInitialConstantDeltaTSteps;

    /*
     *  This is necessary to calculate the "old" variables
     */
    m_func->advanceTimeBaseline(true, m_y_n, m_ydot_n, m_y_n, m_t0, m_t0);

    /*
     * Initialize the various time counters in the object
     */
    time_n = m_t0;
    time_nm1 = time_n;
    time_nm2 = time_nm1;
    delta_t_nm1 = delta_t_n;

    if (m_doSpecialStartCalc) {
        calcConsistentInitialDerivs();
    }
}
//=====================================================================================================================
/*
 * reinitialize():
 *
 */
void BEulerInt::reinitializePRE(m1d::ProblemResidEval& func)
{
    m_NumGlEqns = func.m_neq;
    m_NumLcEqns = func.m_NumLcEqns;
    m_NumLcOwnedEqns = func.m_NumLcOwnedEqns;
    m_ghostedMap = func.ghostedMap();
    m_ownedMap = func.ownedMap();

    m_func = &func;

    internalMalloc();

    /*
     * Initialize the various time counters in the object
     */
    time_n = m_t0;
    time_nm1 = time_n;
    time_nm2 = time_nm1;
    delta_t_n = 0.0;
    delta_t_nm1 = 0.0;

    m_ghostedMap = func.ghostedMap();
    m_ownedMap = func.ownedMap();

    /**
     * Set up the internal weights that are used for testing convergence
     */
    setSolnWeights();

    // Store a pointer to the function
    m_func = &func;
}
//=====================================================================================================================
// Returns the next time of a printout
/*
 * The variable m_printSolnNumberToTout causes the program to print out the solution at even intervals
 * from the beginning of the simulation to the end. This program will return the time for the next printout.
 *   @param   time_current  Current time
 *   @return  returns the next print out time
 */
double BEulerInt::getPrintTime(const double time_current) const
{
    double tnext;
    if (m_printSolnNumberToTout > 0) {
        double dt = (m_time_final - m_t0) / m_printSolnNumberToTout;
        for (int i = 0; i <= m_printSolnNumberToTout; i++) {
            tnext = m_t0 + dt * i;
            if (tnext >= time_current)
                return tnext;
        }
    }
    return 1.0E300;
}
//=====================================================================================================================
/*
 * nEvals():
 * Return the total number of function evaluations
 */
int BEulerInt::nEvals() const
{
    return m_nfe;
}
//=====================================================================================================================
/*
 *
 * internalMalloc():
 *
 *  Internal routine that sets up the fixed length storage based on
 *  the size of the problem to solve.
 */
void BEulerInt::internalMalloc()
{
    Epetra_DataAccess eee = View;
    double *V;

    m_ewt = new Epetra_Vector(*m_ownedMap);
    m_ewt->PutScalar(0.0);

    m_y_n = new Epetra_Vector(*m_ghostedMap);
    m_y_n->PutScalar(0.0);

    m_y_n->ExtractView(&V);
    m_y_n_owned = new Epetra_Vector(eee, *m_ownedMap, V);

    m_y_nm1 = new Epetra_Vector(*m_ghostedMap);
    m_y_nm1->PutScalar(0.0);

    m_y_pred_n = new Epetra_Vector(*m_ghostedMap);
    m_y_pred_n->PutScalar(0.0);

    m_ydot_n = new Epetra_Vector(*m_ghostedMap);
    m_ydot_n->PutScalar(0.0);

    m_ydot_n->ExtractView(&V);
    m_ydot_n_owned = new Epetra_Vector(eee, *m_ownedMap, V);

    m_ydot_nm1 = new Epetra_Vector(*m_ghostedMap);
    m_ydot_nm1->PutScalar(0.0);

    m_resid = new Epetra_Vector(*m_ownedMap);
    m_resid->PutScalar(0.0);

    m_isAlgebraic = new Epetra_IntVector(*m_ghostedMap);
    m_isAlgebraic->PutValue(0);

    m_residWts = new Epetra_Vector(*m_ownedMap);
    m_residWts->PutScalar(0.0);

    m_wksp = new Epetra_Vector(*m_ownedMap);
    m_wksp->PutScalar(0.0);

    m_rowScales = new Epetra_Vector(*m_ownedMap);
    m_rowScales->PutScalar(1.0);

    m_colScales = new Epetra_Vector(*m_ownedMap);
    m_colScales->PutScalar(1.0);

    m_abstol = new Epetra_Vector(*m_ownedMap);
    m_abstol->PutScalar(1.0E-9);

    m_isArithmeticScaled = new Epetra_IntVector(*m_ownedMap);
    for (int i = 0; i < m_NumLcOwnedEqns; i++) {
        (*m_isArithmeticScaled)[i] = 0;
    }
    m_func->fillIsArithmeticScaled(*m_isArithmeticScaled);
}
//=====================================================================================================================
/*
 * setSolnWeights():
 *
 * Set the solution weights
 *  This is a very important routine as it affects quite a few
 *  operations involving convergence.
 */
void BEulerInt::setSolnWeights()
{
    Epetra_IntVector &isS = *m_isArithmeticScaled;
    for (int i = 0; i < m_NumLcOwnedEqns; i++) {
        if (isS[i] == 1) {
            (*m_ewt)[i] = m_reltol * (*m_abstol)[i];
        } else {
            (*m_ewt)[i] = (*m_abstol)[i] + m_reltol * 0.5 * (fabs( (*m_y_nm1)[i]) + fabs( (*m_y_pred_n)[i]));
        }
    }
}
//=====================================================================================================================
/*
 * setColumnScales():
 *
 * Set the column scaling vector at the current time
 */
void BEulerInt::setColumnScales()
{
    m_nonlin->getColumnScaleVector(*m_colScales);
}
//=====================================================================================================================
void BEulerInt::boundStep_y_pred()
{

}
//=====================================================================================================================
//    L2 Weighted Norm of a delta in the solution
/*
 *   The vector m_ewt[i]'s are always used to weight the solution errors in
 *   the calculation.
 *
 *   The second argument has a default of false. However,
 *   if true, then a table of the largest values is printed
 *   out to standard output.
 *
 *   @param delta_y  Norm of a delta of the solution vector
 *   @param printLargest if True a table is printed of the largest contributors.
 *   @param typeYsoln  Parameter indicating whether m_y_curr[] is currently
 *                     evaluated before the delta or after the delta has been implemented.
 *   @param dampFactor Damping factor that will be applied to delta_y before creating a new ysoln
 */
double BEulerInt::soln_error_norm(const Epetra_Vector &delta_y, const int printLargest, const int typeYsoln,
                                  const double dampFactor) const
{
    int i, idLocalEqnMax=-1;
    double sum_norm = 0.0, error;
    stream0 ss;
    double gbSum = 0.0, gmax1;
    for (i = 0; i < m_NumLcOwnedEqns; i++) {
        error = delta_y[i] / (*m_ewt)[i];
        sum_norm += (error * error);
    }
    Comm_ptr_->SumAll(&sum_norm, &gbSum, 1);
    sum_norm = sqrt(gbSum / m_NumGlEqns);
    if (printLargest) {
        if (typeYsoln != 1) {
            printf("not implemented\n");
            exit(-1);
        }
        const int num_entries = 8;
        double dmax1, normContrib;
        int j;
        int *imax = mdpUtil::mdp_alloc_int_1(num_entries, -1);
        print0_sync_start(false, ss, *Comm_ptr_);
        if (!mypid_) {
            printf("\t  ");
            print_line("-", 90);
            printf("\t  soln_error_norm(): delta soln L2 norm = %12.4g\n", sum_norm);
            printf("\t\t      Printout of Largest Contributors:\n");
            printf("\t\t                                                      (damp = %g)\n", dampFactor);
            printf("\t\t      I VarName LcNode weightdeltaY/sqtN|     deltaY      ysolnOld     ysolnNew   Soln_Weights\n");
            printf("\t\t   ");
            print_line("-", 80);
        }
        print0_sync_end(false, ss, *Comm_ptr_);
        for (int jnum = 0; jnum < num_entries; jnum++) {
            dmax1 = -1.0;
            // pick out the contribution for this processor
            for (i = 0; i < m_NumLcOwnedEqns; i++) {
                bool used = false;
                for (j = 0; j < jnum; j++) {
                    if (imax[j] == i)
                        used = true;
                }
                if (!used) {
                    error = delta_y[i] / (*m_ewt)[i];
                    normContrib = sqrt(error * error);
                    if (normContrib > dmax1) {
                        idLocalEqnMax = i;
                        dmax1 = normContrib;
                    }
                }
            }
            int procWinner = procChoice_Max(dmax1, *Comm_ptr_, gmax1);
            print0_sync_start(false, ss, *Comm_ptr_);
            if (procWinner == mypid_) {
                imax[jnum] = idLocalEqnMax;
                i = idLocalEqnMax;
                int idGlobalEqnMax = m_func->LI_ptr_->IndexGbEqns_LcEqns[idLocalEqnMax];
                if (i >= 0) {
                    error = delta_y[i] / (*m_ewt)[i];
                    int iLcNode;
                    int iGbNode;
                    int iNodeEqnNum;
                    VarType var;
                    VAR_TYPE_SUBNUM vtsub;
                    std::string vstring = m_func->variableID(i, iLcNode, iGbNode, iNodeEqnNum, var, vtsub);
                    string v16 = var.VariableName(16);
                    ss.print0("\t\t   %4d %24s-%-4d    %12.4e  | %12.4e %12.4e %12.4e %12.4e\n", idGlobalEqnMax, v16.c_str(),
                            iLcNode, error / sqrt(m_NumGlEqns), delta_y[idLocalEqnMax], (*m_y_n)[idLocalEqnMax],
                            (*m_y_n)[idLocalEqnMax] + delta_y[idLocalEqnMax] * dampFactor, (*m_ewt)[idLocalEqnMax]);
                }
            }
            print0_sync_end(false, ss, *Comm_ptr_);
        }
        print0_sync_start(false, ss, *Comm_ptr_);
        if (!mypid_) {
            printf("\t\t   ");
            print_line("-", 80);
            printf("\t  ");
            print_line("-", 90);
        }
        print0_sync_end(false, ss, *Comm_ptr_);
        mdpUtil::mdp_safe_free((void **) &imax);
    }
    return sum_norm;
}
//=====================================================================================================================
//    L2 Weighted Norm of the residual
/*
 *   The vector m_residWts[i]'s are always used to weight the residual errors in
 *   the calculation.
 *
 *   The second argument has a default of false. However,
 *   if true, then a table of the largest values is printed
 *   out to standard output.
 *
 *   @param delta_y  Norm of a delta of the solution vector
 *   @param printLargest if True a table is printed of the largest contributors.
 */
double BEulerInt::res_error_norm(const Epetra_Vector_Owned &resid, const char *title, const int printLargest) const
{
    return m_nonlin->res_error_norm(resid, title, printLargest);
}
//====================================================================================================================
// computeResidWts():
/*
 *  The residual weights are defined here to be equal to the inverse of the row scaling factors used to
 *  row scale the matrix, after column scaling is used. They are multiplied by 10-3 because the column
 *  weights are also multiplied by that same quantity.
 *
 *  The basic idea is that a change in the solution vector on the order of the convergence tolerance
 *  multiplied by  [RJC] which is of order one after row scaling should give you the relative weight
 *  of the row. Values of the residual for that row can then be normalized by the value of this weight.
 *  When the tolerance in delta x is achieved, the tolerance in the residual is also achieved.
 *
 * Here a small weighting indicates that the change in solution is
 * very sensitive to that equation.
 *
 * There are several prerequisites.
 *    Row and column scaling must be used on the matrix.
 *    The Jacobian must have already been formed.
 */
void BEulerInt::computeResidWts()
{
    m_nonlin->getResidWts(*m_residWts);
}
//====================================================================================================================
/*
 * filterNewStep():
 *
 * void BEulerInt::
 *
 */
double BEulerInt::filterNewStep(double timeCurrent, Epetra_Vector_Ghosted &y_current, Epetra_Vector_Ghosted &ydot_current)
{
    /*
     *  Call the residual function's filter program to calculate a
     *  delta_y. It's stored in m_wksp.
     */
    m_func->applyFilter(timeCurrent, delta_t_n, y_current, ydot_current, *m_wksp);
    double s1 = soln_error_norm(*m_wksp, false);
    if (s1 > 0.0) {
        if (m_print_flag > 2) {
            //m_func->writeSolutionFilter(1, timeCurrent, y_current, &ydot_current, *m_wksp, "SkewNegSectFilter");
        } else {
            if (m_print_flag > 1 && s1 >= 1.0) {
                //  m_func->writeSolutionFilter(1, timeCurrent, y_current, &ydot_current, *m_wksp, "SkewNegSectFilter");
            }
        }
        for (int j = 0; j < m_NumLcEqns; j++) {
            y_current[j] += (*m_wksp)[j];
        }
        calc_ydot(m_order, y_current, ydot_current);
    }
    return s1;
}
//====================================================================================================================
void BEulerInt::print_time_step1(int order, int n_time_step, double time, double delta_t_n, double delta_t_nm1, bool step_failed,
                                 int num_failures) const
                                 /*
                                  * Print out for relevant time step information
                                  */
{
    if (!mypid_) {
        const char *string = 0;
        if (order == 0)
            string = "Backward Euler";
        else if (order == 1)
            string = "Forward/Backward Euler";
        else if (order == 2)
            string = "Adams-Bashforth/TR";
        printf("\n");
        print_line("=", 80);
        printf("\nStart of Time Step: %5d       Time_n = %9.5g Time_nm1 = %9.5g\n", n_time_step, time, time - delta_t_n);
        printf("\tIntegration method = %s\n", string);
        if (step_failed)
            printf("\tPreviously attempted step was a failure\n");
        if (delta_t_n > delta_t_nm1)
            string = "(Increased from previous iteration)";
        else if (delta_t_n < delta_t_nm1)
            string = "(Decreased from previous iteration)";
        else {
            string = "(same as previous iteration)";
        }
        printf("\tdelta_t_n        = %8.5e %s", delta_t_n, string);
        if (num_failures > 0)
            printf("\t(Bad_History Failure Counter = %d)", num_failures);
        printf("\n\tdelta_t_nm1      = %8.5e\n", delta_t_nm1);
    }
}
//====================================================================================================================
/*
 * Print out for relevant time step information
 */
void BEulerInt::print_time_step2(int time_step_num, int order, double time, double time_error_factor, double delta_t_n,
                                 double delta_t_np1) const
{
    if (!mypid_) {
        printf("\tTime Step Number %5d was a success: time = %10g\n", time_step_num, time);
        printf("\t\tEstimated Error\n");
        printf("\t\t--------------------   =   %8.5e\n", time_error_factor);
        printf("\t\tTolerated Error\n\n");
        printf("\t- Recommended next delta_t (not counting history) = %g\n", delta_t_np1);
        printf("\n");
        print_line("=", 80);
        printf("\n");
    }
}
//====================================================================================================================
/*
 * Print Out descriptive information on why the current step failed
 */
void BEulerInt::print_time_fail(bool convFailure, int time_step_num, double time, double delta_t_n, double delta_t_np1,
                                double time_error_factor) const
{
    if (!mypid_) {
        printf("\n");
        print_line("=", 80);
        if (convFailure) {
            printf("\tTime Step Number %5d experienced a convergence "
                    "failure\n", time_step_num);
            printf("\tin the non-linear or linear solver\n");
            printf("\t\tValue of time at failed step           = %g\n", time);
            printf("\t\tdelta_t of the   failed step           = %g\n", delta_t_n);
            printf("\t\tSuggested value of delta_t to try next = %g\n", delta_t_np1);
        } else {
            printf("\tTime Step Number %5d experienced a truncation error "
                    "failure!\n", time_step_num);
            printf("\t\tValue of time at failed step           = %g\n", time);
            printf("\t\tdelta_t of the   failed step           = %g\n", delta_t_n);
            printf("\t\tSuggested value of delta_t to try next = %g\n", delta_t_np1);
            printf("\t\tCalculated truncation error factor  = %g\n", time_error_factor);
        }
        printf("\n");
        print_line("=", 80);
    }
}
//====================================================================================================================
/*
 * Print out the final results and counters
 */
void BEulerInt::print_final(double time, int step_failed, int time_step_num, int num_newt_its, int total_linear_solves,
                            int numConvFails, int numTruncFails, int nfe, int nJacEval) const
{
    if (!mypid_) {
        printf("\n");
        print_line("=", 80);
        printf("TIME INTEGRATION ROUTINE HAS FINISHED: ");
        if (step_failed)
            printf(" IT WAS A FAILURE\n");
        else
            printf(" IT WAS A SUCCESS\n");
        printf("\tEnding time                   = %g\n", time);
        printf("\tNumber of time steps          = %d\n", time_step_num);
        printf("\tNumber of newt its            = %d\n", num_newt_its);
        printf("\tNumber of linear solves       = %d\n", total_linear_solves);
        printf("\tNumber of convergence failures= %d\n", numConvFails);
        printf("\tNumber of TimeTruncErr fails  = %d\n", numTruncFails);
        printf("\tNumber of Function evals      = %d\n", nfe);
        printf("\tNumber of Jacobian evals/solvs= %d\n", nJacEval);
        printf("\n");
        print_line("=", 80);
    }
}
//====================================================================================================================
/*
 * Header info for one line comment about a time step
 */
void BEulerInt::print_lvl1_Header(int nTimes) const
{
    if (!mypid_) {
        printf("\n");
        if (nTimes) {
            print_line("-", 80);
        }
        printf("time       Time              Time                     Time  ");
        if (nTimes == 0) {
            printf("     START");
        } else {
            printf("    (continued)");
        }
        printf("\n");

        printf("step      (sec)              step  Newt   Aztc bktr  trunc  ");
        printf("\n");

        printf(" No.               Rslt      size    Its  Its  stps  error     |");
        printf("  comment");
        printf("\n");
        print_line("-", 80);
    }
}
//====================================================================================================================
void BEulerInt::print_lvl3_Header() const
{
    print_SVector("initial Atol vector", *m_abstol);
}
//====================================================================================================================
/*
 * One line entry about time step
 *   rslt -> 4 letter code
 */
void BEulerInt::print_lvl1_summary(int time_step_num, double time, const char *rslt, double delta_t_n, int newt_its, int aztec_its,
                                   int bktr_stps, double time_error_factor, const char *comment) const
{
    if (!mypid_) {
        printf("%6d %11.6g %4s %10.4g %4d %4d %4d %11.4g", time_step_num, time, rslt, delta_t_n, newt_its, aztec_its, bktr_stps,
                time_error_factor);
        if (comment)
            printf(" | %s", comment);
        printf("\n");
    }
}
//====================================================================================================================
void BEulerInt::print_SVector(std::string header, const Epetra_MultiVector &v) const
{
    if (!mypid_) {
        printf("%s\n", header.c_str());
    }

    Epetra_VbrMatrix *A = (*tdjac_ptr).A_;
    const Epetra_BlockMap &domainMap = A->DomainMap();
    Epetra_Vector_Owned *vv = new Epetra_Vector(domainMap);

    for (int i = 0; i < m_NumLcOwnedEqns; i++) {
        double *a = v[0];
        (*vv)[i] = a[i];
    }

    m_func->showSolutionVector(header, time_n, delta_t_n, *vv);

    delete vv;
    // stream0 ss;
    //Print0_epMultiVector(ss, v);
}
//====================================================================================================================
void BEulerInt::print_IntSVector(std::string header, const Epetra_IntVector &v) const
{
    if (!mypid_) {
        printf("%s\n", header.c_str());
    }
    m_func->showSolutionIntVector(header, time_n, delta_t_n, v);

    // stream0 ss;
    //Print0_epMultiVector(ss, v);
}
//====================================================================================================================
/*
 * Function to calculate the predicted solution vector, m_y_pred_n for the
 * (n+1)th time step.
 *
 * This routine can be used by a first order - forward
 * Euler / backward Euler predictor / corrector method or for a second order
 * Adams-Bashforth / Trapezoidal Rule predictor / corrector method.  See Nachos
 * documentation Sand86-1816 and Gresho, Lee, Sani LLNL report UCRL - 83282 for
 * more information.
 *
 * variables:
 *
 * on input:
 *
 *     N          - number of unknowns
 *     order      - indicates order of method
 *                  = 1 -> first order forward Euler/backward Euler
 *                         predictor/corrector
 *                  = 2 -> second order Adams-Bashforth/Trapezoidal Rule
 *                         predictor/corrector
 *
 *    delta_t_n   - magnitude of time step at time n     (i.e., = t_n+1 - t_n)
 *    delta_t_nm1 - magnitude of time step at time n - 1 (i.e., = t_n - t_n-1)
 *    y_nm1[]       - solution vector at time nm1
 *    y_dot_nm1[]   - acceleration vector from the predictor at time nm1
 *    y_dot_nm2[] - acceleration vector from the predictor at time nm1 - 1
 *
 * on output:
 *
 *    m_y_pred_n[]    - predicted solution vector at time n
 */
void BEulerInt::calc_y_pred(const int order)
{
    int i;
    double c1;
    if (!m_DOT_Predicted_Bad) {
        switch (order) {
        case 0:
        case 1:
            c1 = delta_t_n;
            for (i = 0; i < m_NumLcEqns; i++) {
                if ( (*m_isAlgebraic)[i] == 1) {
                    if (!m_DAEDOT_Predicted_Bad) {
                        (*m_y_pred_n)[i] = (*m_y_nm1)[i] + c1 * (*m_ydot_nm1)[i];
                    } else {
                        (*m_y_pred_n)[i] = (*m_y_nm1)[i];
                    }
                } else {
                    (*m_y_pred_n)[i] = (*m_y_nm1)[i] + c1 * (*m_ydot_nm1)[i];
                }
            }
            break;
        case 2:
            exit(-1);
        }
    } else {
        for (i = 0; i < m_NumLcEqns; i++) {
            (*m_y_pred_n)[i] = (*m_y_nm1)[i];
        }
    }
    /*
     * Filter the predictions.
     */
    m_func->filterSolnPrediction(time_n, *m_y_pred_n);
}
//====================================================================================================================
/* Function to calculate the acceleration vector ydot for the first or
 * second order predictor/corrector time integrator.  This routine can be
 * called by a first order - forward Euler / backward Euler predictor /
 * corrector or for a second order Adams - Bashforth / Trapezoidal Rule
 * predictor / corrector.  See Nachos documentation Sand86-1816 and Gresho,
 * Lee, Sani LLNL report UCRL - 83282 for more information.
 *
 *  variables:
 *
 *    on input:
 *
 *       N          - number of local unknowns on the processor
 *                    This is equal to internal plus border unknowns.
 *       order      - indicates order of method
 *                    = 1 -> first order forward Euler/backward Euler
 *                           predictor/corrector
 *                    = 2 -> second order Adams-Bashforth/Trapezoidal Rule
 *                           predictor/corrector
 *
 *      delta_t_n   - Magnitude of the current time step at time n
 *                    (i.e., = t_n - t_n-1)
 *      y_curr[]    - Current Solution vector at time n
 *      y_nm1[]     - Solution vector at time n-1
 *      ydot_nm1[] - Acceleration vector at time n-1
 *
 *   on output:
 *
 *      ydot_curr[]   - Current acceleration vector at time n
 *
 * Note we use the current attribute to denote the possibility that
 * y_curr[] may not be equal to m_y_n[] during the nonlinear solve
 * because we may be using a look-ahead scheme.
 */
void BEulerInt::calc_ydot(int order, m1d::Epetra_Vector_Ghosted & y_curr, m1d::Epetra_Vector_Ghosted & ydot_curr)
{
    int i;
    double c1;
    switch (order) {
    case 0:
    case 1: /* First order forward Euler/backward Euler */
        c1 = 1.0 / delta_t_n;
        for (i = 0; i < m_NumLcEqns; i++) {
            ydot_curr[i] = c1 * (y_curr[i] - (*m_y_nm1)[i]);
        }
        return;
    case 2: /* Second order Adams-Bashforth / Trapezoidal Rule */
        c1 = 2.0 / delta_t_n;
        for (i = 0; i < m_NumLcEqns; i++) {
            ydot_curr[i] = c1 * (y_curr[i] - (*m_y_nm1)[i]) - (*m_ydot_nm1)[i];
        }
        return;
    }
}
//=====================================================================================================================
// This function calculates the time step error norm using an L2 formula
/*
 *    on input:
 *
 *      m_y_n        -  Actual value of the solution
 *      m_y_pred_n   -  Predicted value of the solution
 *      m_ewt        - Weight vector
 *
 *   Returns the l2 norm of (m_y_n[i] - m_y_pred_n[i]) / m_ewt[i]
 */

double BEulerInt::time_error_norm() const
{
    stream0 ss;
    double gbSum = 0.0, gmax1;
    int i;
    double sum_norm, error, tfac;
    double dd = MAX(delta_t_nm1, delta_t_nm2);
    tfac = delta_t_n / dd;
    if (tfac > 1.0)
        tfac = 1.0;
    double afac = tfac * tfac;
    sum_norm = 0.0;
    for (i = 0; i < m_NumLcOwnedEqns; i++) {
        error = ( (*m_y_n)[i] - (*m_y_pred_n)[i]) / (*m_ewt)[i];
        if ( (*m_isAlgebraic)[i]) {
            error *= afac;
        }
        sum_norm += (error * error);
    }
    Comm_ptr_->SumAll(&sum_norm, &gbSum, 1);
    sum_norm = sqrt(gbSum / m_NumGlEqns);

    if (m_print_flag > 2) {
        if (!mypid_) {
            printf("\tTime step truncation error = %g\n", sum_norm);
        }
        if (m_print_flag > 4) {

            const int num_entries = 8;
            double dmax1, normContrib;
            int j;
            int idLocalEqnMax=-1;
            int *imax = mdpUtil::mdp_alloc_int_1(num_entries, -1);
            print0_sync_start(false, ss, *Comm_ptr_);
            if (!mypid_) {
                printf("\t\tTime step truncation error contributors\n");
                printf(
                        "\t\t    I               VarName  LcNode       entry    |   y_actual     y_predicted      weight       ydot\n");
                printf("\t\t");
                print_line("-", 80);
            }
            print0_sync_end(false, ss, *Comm_ptr_);
            for (int jnum = 0; jnum < num_entries; jnum++) {
                dmax1 = -1.0;
                // pick out the contribution for this processor
                for (i = 0; i < m_NumLcOwnedEqns; i++) {
                    bool used = false;
                    for (j = 0; j < jnum; j++) {
                        if (imax[j] == i)
                            used = true;
                    }
                    if (!used) {
                        error = ( (*m_y_n)[i] - (*m_y_pred_n)[i]) / (*m_ewt)[i];
                        if ( (*m_isAlgebraic)[i]) {
                            error *= afac;
                        }
                        normContrib = sqrt(error * error);
                        if (normContrib > dmax1) {
                            idLocalEqnMax = i;
                            dmax1 = normContrib;
                        }
                    }
                }
                int procWinner = procChoice_Max(dmax1, *Comm_ptr_, gmax1);
                print0_sync_start(false, ss, *Comm_ptr_);
                if (procWinner == mypid_) {
                    imax[jnum] = idLocalEqnMax;
                    i = idLocalEqnMax;
                    //int i_gb = m_func->LI_ptr_->IndexGbNode_LcNode[lcmax];
                    int idGlobalEqnMax = m_func->LI_ptr_->IndexGbEqns_LcEqns[idLocalEqnMax];
                    if (i >= 0) {
                        error = ( (*m_y_n)[idLocalEqnMax] - (*m_y_pred_n)[idLocalEqnMax]) / (*m_ewt)[idLocalEqnMax];
                        string saa = "   ";
                        double aafac = 1.0;
                        if ( (*m_isAlgebraic)[i]) {
                            saa = "alg";
                            aafac = afac;
                        }
                        int iLcNode;
                        int iGbNode;
                        int iNodeEqnNum;
                        VarType var;
                        VAR_TYPE_SUBNUM vtsub;
                        std::string vstring = m_func->variableID(i, iLcNode, iGbNode, iNodeEqnNum, var, vtsub);
                        string v16 = var.VariableName(16);
                        ss.print0("\t\t  %4d %24s-%-4d  %12.4e | %12.4e %12.4e %12.4e %3s %12.4e\n", idGlobalEqnMax, v16.c_str(),
                                iLcNode, error / sqrt(m_NumGlEqns) * aafac, (*m_y_n)[idLocalEqnMax], (*m_y_pred_n)[idLocalEqnMax],
                                (*m_ewt)[idLocalEqnMax], saa.c_str(), (*m_ydot_n)[idLocalEqnMax]);
                    }
                }
                print0_sync_end(false, ss, *Comm_ptr_);
            }
            print0_sync_start(false, ss, *Comm_ptr_);
            if (!mypid_) {
                printf("\t\t   ");
                print_line("-", 80);
                printf("\t  ");
                print_line("-", 90);
            }
            print0_sync_end(false, ss, *Comm_ptr_);
            mdpUtil::mdp_safe_free((void **) &imax);
        }
    }
    return sum_norm;
}
//=====================================================================================================================

// Time step control function for the selection of the time step size
/*
 * The algorithm is based on a desired accuracy of time integration and on an estimate of the relative
 * error of the time integration process. This routine can be called for a
 * first order - forward Euler/backward Euler predictor/ corrector and for a
 * second order Adams- Bashforth/Trapezoidal Rule predictor/corrector. See
 * Nachos documentation Sand86-1816 and Gresho, Lee, Sani LLNL report UCRL -
 * 83282 for more information.
 *
 *  variables:
 *
 *    on input:
 *
 *       order      - indicates order of method
 *                    = 1 -> first order forward Euler/backward Euler
 *                           predictor/corrector
 *                    = 2 -> second order forward Adams-Bashforth/Trapezoidal
 *                          rule predictor/corrector
 *
 *      delta_t_n   - Magnitude of time step at time t_n
 *      delta_t_nm1 - Magnitude of time step at time t_n-1
 *      rel_error   - Generic relative error tolerance
 *      time_error_factor   - Estimated value of the time step truncation error
 *                           factor. This value is a ratio of the computed
 *                           error norms. The premultiplying constants
 *                           and the power are not yet applied to normalize the
 *                           predictor/corrector ratio. (see output value)
 *
 *   on output:
 *
 *      return - delta_t for the next time step
 *               If delta_t is negative, then the current time step is
 *               rejected because the time-step truncation error is
 *               too large.  The return value will contain the negative
 *               of the recommended next time step.
 *
 *      time_error_factor  - This output value is normalized so that
 *                           values greater than one indicate the current time
 *                           integration error is greater than the user
 *                           specified magnitude.
 */
double BEulerInt::time_step_control(int order, double time_error_factor) const
{
    double factor = 0.0, power = 0.0, delta_t;
    const char *yo = "time_step_control";

    /*
     * Special case time_error_factor so that zeroes don't cause a problem.
     */
    time_error_factor = MAX(1.0E-50, time_error_factor);

    /*
     * Calculate the factor for the change in magnitude of time step.
     */
    switch (order) {
    case 1:
        factor = 1.0 / (2.0 * (time_error_factor));
        power = 0.5;
        break;
    case 2:
        factor = 1.0 / (3.0 * (1.0 + delta_t_nm1 / delta_t_n) * (time_error_factor));
        power = 0.3333333333333333;
        break;
    }
    factor = pow(factor, power);
    if (factor < 0.5) {
        if (m_print_flag > 1) {
            printf("\t%s: WARNING - Current time step will be chucked\n", yo);
            printf("\t\tdue to a time step truncation error failure.\n");
        }
        delta_t = -0.5 * delta_t_n;
    } else {
        factor = MIN(factor, 1.5);
        delta_t = factor * delta_t_n;
    }
    return delta_t;
}
//=====================================================================================================================

/*
 * @param   val
 *  @param start Is the input a starting point. If so, we return the greater region
 *  @return returns the region
 */
int BEulerInt::findTimeRegion(double val, bool start)
{
    int i;
    if (m_timeRegionBoundaries.size() < 2) {
        return m_currentTimeRegion;
    }
    if (val < m_timeRegionBoundaries[0]) {
        throw m1d_Error("BEulerInt::findTimeRegion()", "out of bounds");
    }
    for (i = 1; i < (int) m_timeRegionBoundaries.size(); i++) {
        if (start) {
            if (val < m_timeRegionBoundaries[i]) {
                return i - 1;
            }
        } else {
            if (val <= m_timeRegionBoundaries[i]) {
                return i - 1;
            }
        }
    }
    throw m1d_Error("BEulerInt::findTimeRegion()", "out of bounds");
    return i;
}

//=====================================================================================================================
/*
 *
 * integrate():
 *
 *  defaults are located in the .h file. They are as follows:
 *     time_init = 0.0
 */
double BEulerInt::integratePRE(double tout)
{
    double time_current;
    bool weAreNotFinished = true;
    m_time_final = tout;
    int flag = BE_SUCCESS;

    /**
     * Initialize the time step number to zero. step will increment so that
     * the first time step is number 1
     */
    m_time_step_num = 0;
    m_time_step_attempts = 0;

    /*
     * Do the integration a step at a time
     */
    int istep = 0;
    int printStep = 0;
    double print_time = m_t0;
    bool doPrintSoln = false;
    time_current = m_t0;
    time_n = m_t0;
    time_nm1 = m_t0;
    time_nm2 = m_t0;

    /*
     *   Check for output time crossing a time region
     */
    if (m_timeRegionBoundaries.size() > 0) {
        int startRegion = findTimeRegion(m_t0, true);
        int endRegion = findTimeRegion(tout, false);
        m_currentTimeRegion = startRegion;
        if (endRegion != startRegion) {
            tout = m_timeRegionBoundaries[m_currentTimeRegion + 1];
        }
    }
    m_func->m_currentTimeRegion = m_currentTimeRegion;
    m_func->m_numTimeRegions = m_numTimeRegions;


    /*
     * Install the atol vector from the ProblemResidFunc.
     */
    const Epetra_Vector_Owned & aabstol = m_func->atolVector();
    setTolerancesEpetra(m_reltol, aabstol);

    /*
     *  HKM -> Not sure that this call is necessary. Consider taking it out
     */
    m_func->advanceTimeBaseline(true, m_y_n, m_ydot_n, m_y_nm1, time_current, time_nm1);
    m_func->evalTimeTrackingEqns(0, time_current, 0.0, *m_y_n, m_ydot_n);

    /*
     *   Always call writeSolution to write out the initial conditions
     *    -> we use a special flag for the solution type, because we may not have solved a problem 
     *       at t_nm1
     */
    if (m_print_flag > 0) {
        m_func->writeSolution(0, true, time_current, delta_t_n, istep, *m_y_n, m_ydot_n, 
			      TimeDependentInitial, delta_t_np1);
    }

    /*
     * We print out column headers here for the case of
     */
    if (m_print_flag == 1) {
        print_lvl1_Header(0);
    }
    /*
     * Call a different user routine at the end of each step,
     * that will probably print to a file.
     */
    m_func->user_out(0, time_current, 0.0, istep, *m_y_n, m_ydot_n);

#ifdef DEBUG_MODE
    print_lvl3_Header();
#else
    if (m_print_flag > 2) {
        print_lvl3_Header();
    }
#endif
    //
    //   MAIN LOOP OVER THE TIME STEPS
    //
    do {

        print_time = getPrintTime(time_current);
        if (print_time >= tout) {
            print_time = tout;
        }

        // Call the time baseline advancement routine
        /*
         *  This function provides a hook for a residual that gets called whenever a
         *  time step has been accepted and we are about to move on to the next time step.
         *  The call is made with the current time as the time that is accepted. The old time may be
         *  obtained from t and rdelta_t_accepted.
         *  After this call interrogation of the previous time step's results will not be valid.
         */
        m_func->advanceTimeBaseline(true, m_y_n, m_ydot_n, m_y_nm1, time_current, time_nm1);

        m_func->m_currentTimeRegion = m_currentTimeRegion;

        /******************************************************************************************************
         *                                  Step the solution
         ******************************************************************************************************/
        istep++;
        printStep++;
        // Increment the step number in the function object
        m_func->m_StepNumber++;

        time_current = step(tout);

        /*******************************************************************************************************/

        if (time_current < 0.0) {
            if (time_current == -1234.) {
                time_current = 0.0;
            } else {
                time_current = -time_current;
            }
            flag = BE_FAILURE;
        }

        if (flag != BE_FAILURE) {
            bool retn = m_func->evalStoppingCritera(time_current, delta_t_n, *m_y_n, *m_ydot_n);
            if (retn) {
                weAreNotFinished = false;
                doPrintSoln = true;
            }
        }

        /*
         * Determine conditional printing of soln
         */
        if (time_current >= print_time) {
            doPrintSoln = true;
        }
        if (m_printSolnStepInterval == printStep) {
            doPrintSoln = true;
        }
        if (m_printSolnFirstSteps > istep) {
            doPrintSoln = true;
        }

        /*
         * Evaluate time integrated quantities that are calculated at the
         * end of every successful time step.
         */
        if (flag != BE_FAILURE) {
            m_func->evalTimeTrackingEqns(1, time_current, delta_t_n, *m_y_n, m_ydot_n);
        }

        /*
         * Call the printout routine.
	 *    -> We included the expected delta t for the next step, delta_t_np1, into the call so that
	 *       we can have good restart capability and reproducibility.
         */
        if (doPrintSoln && flag != BE_FAILURE) {
            m_func->writeSolution(1, true, time_current, delta_t_n, istep, *m_y_n, m_ydot_n,
				  TimeDependentAccurate_Solve, delta_t_np1);
            printStep = 0;
            doPrintSoln = false;
            if (m_print_flag == 1) {
                print_lvl1_Header(1);
            }
        }
        /*
         * Call a different user routine at the end of each step,
         * that will probably print to a file.
         */
        if (flag == BE_FAILURE) {
            m_func->user_out(-1, time_current, delta_t_n, istep, *m_y_n, m_ydot_n);
        } else {
            m_func->user_out(1, time_current, delta_t_n, istep, *m_y_n, m_ydot_n);
        }

    } while (time_current < tout && m_time_step_attempts < m_max_time_step_attempts && flag == BE_SUCCESS && weAreNotFinished);

    /*
     * Check current time against the max solution time.
     */
    if (!mypid_) {
        if (time_current >= tout) {
            printf("Simulation completed time integration in %d time steps\n", m_time_step_num);
            printf("Final Time: %e\n\n", time_current);
        } else if (m_time_step_attempts >= m_max_time_step_attempts) {
            printf("Simulation ran into time step attempt limit in"
                    "%d time steps\n", m_time_step_num);
            printf("Final Time: %e\n\n", time_current);
        } else if (flag == BE_FAILURE) {
            printf("ERROR: time stepper failed at time = %g\n", time_current);
        }
    }

    /*
     * Print out the final results and counters.
     */
    print_final(time_n, flag, m_time_step_num, m_numTotalNewtIts, m_numTotalLinearSolves, m_numTotalConvFails, m_numTotalTruncFails,
            m_nfe, m_nJacEval);

    /*
     * Call a different user routine at the end of each step,
     * that will probably print to a file.
     */
    if (flag == BE_SUCCESS) {
        m_func->writeSolution(2, false, time_current, delta_t_n, istep, *m_y_n, m_ydot_n, 
			      TimeDependentAccurate_Solve, delta_t_np1);
        m_func->user_out(2, time_current, delta_t_n, istep, *m_y_n, m_ydot_n);
    }
    if (flag != BE_SUCCESS) {
        throw BEulerErr(" BEuler error encountered.");
    }
    return time_current;
}
//=====================================================================================================================
/*
 * step():
 *
 * This routine advances the calculations one step using a predictor
 * corrector approach. We use an implicit algorithm here.
 *
 */
doublereal BEulerInt::step(double t_max)
{
    double CJ = 0.0;
    bool step_failed = false;
    bool giveUp = false;
    bool convFailure = false;
    const char *rslt = "";
    double time_error_factor = 0.0;
    double normFilter = 0.0;
    int numTSFailures = 0;
    int bktr_stps = 0;
    int num_newt_its = 0;
    int aztec_its = 0;
    string comment;
    /*
     * Increment the time counter - May have to be taken back,
     * if time step is found to be faulty.
     */
    m_time_step_num++;

    if (1.E-15 > t_max - time_n) {
        return t_max;
    }

    /*
     * Save the old solution, before overwriting with the new solution
     * - use
     */
    mdpUtil::mdp_copy_dbl_1(& (*m_y_nm1)[0], & ( (*m_y_n)[0]), m_NumLcEqns);

    /*
     * Save the old time derivative, if necessary, before it is
     * overwritten.
     * This overwrites ydot_nm1, losing information from the previous time
     * step.
     */
    mdpUtil::mdp_copy_dbl_1(& (*m_ydot_nm1)[0], & ( (*m_ydot_n)[0]), m_NumLcEqns);

    /*
     * Loop here until we achieve a successful step or we set the giveUp
     * flag indicating that repeated errors have occurred.
     */

    delta_t_nm2 = delta_t_nm1;
    delta_t_nm1 = delta_t_n;
    delta_t_n = delta_t_np1;
    time_nm2 = time_nm1;
    time_nm1 = time_n;
    do {

        m_time_step_attempts++;
        comment.clear();

        if (delta_t_n <= delta_t_min) {
            delta_t_n = delta_t_min;
        }

        /*
         * Possibly adjust the delta_t_n value for this time step from the
         * recommended delta_t_np1 value determined in the previous step
         * due to maximum time step constraints or other occurences,
         * known to happen at a given time.
         */
        if ( (time_n + delta_t_n) >= t_max) {
            delta_t_n = t_max - time_n;
        }

        if (delta_t_n >= delta_t_max) {
            delta_t_n = delta_t_max;
        }
        //delta_t_n_old = delta_t_n;

        /*
         * Increment the delta_t counters and the time for the current
         * time step.
         */
        time_n = time_nm1 + delta_t_n;

        m_DAEDOT_Predicted_Bad = 0;
        m_DOT_Predicted_Bad = 0;
        m_DAE_Predicted_Bad = 0;
        if (m_time_step_num == 1) {
            if (m_doSpecialStartCalc) {
                m_DAEDOT_Predicted_Bad = 1;
            } else {
                m_DAEDOT_Predicted_Bad = 1;
                m_DOT_Predicted_Bad = 1;
                m_DAE_Predicted_Bad = 1;
            }
        } else if (m_time_step_num == 2) {
            if (!m_doSpecialStartCalc) {
                m_DAEDOT_Predicted_Bad = 1;
                // -> not true but first want to check for algorithm validity
                m_DOT_Predicted_Bad = 1;
            }
        }

        /*
         * Determine the integration order of the current step.
         *
         * Special case for start-up of time integration procedure
         *           First time step = Do a predictor step as we
         *                             have recently added an initial
         *                             ydot input option. And, setting ydot=0
         *                             is equivalent to not doing a
         *                             predictor step.
         *           Second step     = If 2nd order method, do a first order
         *                             step for this time-step, only.
         *
         *           If 2nd order method with a constant time step, the
         *           first and second steps are 1/10 the specified step, and
         *           the third step is 8/10 the specified step.  This reduces
         *           the error asociated with using lower order
         *           integration on the first two steps. (RCS 11-6-97)
         *
         * If the previous time step failed for one reason or another,
         * do a linear step. It's more robust.
         */
        if (m_time_step_num == 1) {
            m_order = 1; /* Backward Euler          */
        } else if (m_time_step_num == 2) {
            m_order = 1; /* Forward/Backward Euler  */
        } else if (step_failed) {
            m_order = 1; /* Forward/Backward Euler  */
        } else if (m_time_step_num > 2) {
            m_order = 1; /* Specified
             Predictor/Corrector
             - not implemented */
        }

        /*
         * Print out an initial statement about the step.
         */
        if (m_print_flag > 1) {
            print_time_step1(m_order, m_time_step_num, time_n, delta_t_n, delta_t_nm1, step_failed, m_failure_counter);
        }

        /*
         * Calculate the predicted solution, m_y_pred_n, for the current
         * time step.
         */
        calc_y_pred(m_order);

        /*
         * Cropping of the predictor could go here ....
         *
         *     cropNorm = 0.0;
         *     if (Cur_Realm->Realm_Nonlinear.Constraint_Backtracking_Flag == Constraint_Backtrack_Enable) {
         *       cropNorm = cropPredictor(mesh, x_pred_n, abs_time_error, m_reltol);
         *     }
         */

        /*
         * Use the predicted value as the initial guess for the corrector
         * loop, for
         * every step other than the first step.
         */
        if (m_order > 0) {
            mdpUtil::mdp_copy_dbl_1(& (*m_y_n)[0], & ( (*m_y_pred_n)[0]), m_NumLcEqns);
        }

        /*
         * Calculate the new time derivative, ydot_n, that is consistent
         * with the
         * initial guess for the corrected solution vector.
         *
         */
        calc_ydot(m_order, *m_y_n, *m_ydot_n);
        /*
         * Calculate CJ, the coefficient for the jacobian corresponding to the
         * derivative of the residual wrt to the acceleration vector.
         */
        if (m_order < 2) {
            CJ = 1.0 / delta_t_n;
	} else {
	    CJ = 2.0 / delta_t_n;
        }
        /*
         * Check to see if the predicted solution is ok
         */
	step_failed = false;
        int ierror = check_predicted_soln(*m_y_n, *m_ydot_n, CJ, time_n);
        if (ierror < 0) {
          step_failed = true;
        }

        if (! step_failed)  {

	    /*
	     * Calculate a new Solution Error Weighting vector
	     */
	    setSolnWeights();
	    
	    /*
	     * We need to seed the nonlinear solver with the predicted solution and with the solution
	     * from the previous time step.
	     */
	    
	    if (m_time_step_num > 2) {
		m_nonlin->setMaxNewtIts(50);
	    } else {
		m_nonlin->setMaxNewtIts(150);
	    }
	    /*
	     *  Pass down the abs tolerances from the time stepper to the nonlinear solver
	     */
	    m_nonlin->setTolerances(m_reltol, m_NumLcEqns, & ( (*m_abstol)[0]));
	    m_nonlin->setPredicted_soln(*m_y_pred_n);
	    m_nonlin->setPreviousTimeStep(delta_t_n, *m_y_nm1, *m_ydot_nm1);
	    
	    const Epetra_Vector_Owned &abs_dd = m_func->atolDeltaDamping();
	    m_nonlin->setTolerances_deltaDamping(m_reltol, m_NumLcOwnedEqns, & (abs_dd[0]));
	    /*
	     * Solve the system of equations at the current time step.
	     * Note - x_corr_n and x_dot_n are considered to be updated,
	     * on return from this solution.
	     */
	    int num_linear_solves;
	    int numbacktracks;
	    m1d::Solve_Type_Enum ss = TimeDependentAccurate_Solve;
	    int ierror = m_nonlin->solve_nonlinear_problem(ss, m_y_n, m_ydot_n, CJ, time_n, num_newt_its, num_linear_solves,
							   numbacktracks);

#ifdef DEBUG_HKM_NOT
	    if ( m_time_step_attempts == 4) {
		ierror = -2;
	    }
	    if ( m_time_step_attempts == 7) {
		ierror = -2;
	    }
#endif
	    
	    /*
	     * Set the appropriate flags if a convergence failure is detected.
	     */
	    if (ierror < 0) { /* Step failed */
		convFailure = true;
		step_failed = true;
		rslt = "fail";
		m_numTotalConvFails++;
		m_failure_counter += 3;
		if (!mypid_ && m_print_flag > 1) {
		    printf("\tStep is Rejected, nonlinear problem didn't converge,"
			   "ierror = %d\n", ierror);
		}
	    } else { /* Step succeeded */
		convFailure = false;
		step_failed = false;
		rslt = "done";
		
		/*
		 *  Apply a filter to a new successful step
		 */
		normFilter = filterNewStep(time_n, *m_y_n, *m_ydot_n);
		if (normFilter > 1.0) {
		    convFailure = true;
		    step_failed = true;
		    rslt = "filt";
		    if (!mypid_ && m_print_flag > 1) {
			printf("\tStep is Rejected, too large filter adjustment = %g\n", normFilter);
		    }
		} else if (normFilter > 0.0) {
		    if (normFilter > 0.3) {
			if (!mypid_ && m_print_flag > 1) {
			    printf("\tStep was filtered, norm = %g, next "
				   "time step adjusted\n", normFilter);
			}
		    } else {
			if (!mypid_ && m_print_flag > 1) {
			    printf("\tStep was filtered, norm = %g\n", normFilter);
			}
		    }
		}
	    }
        }

        /*
         * Calculate the time step truncation error for the current step.
         */
        if (!step_failed) {
            time_error_factor = time_error_norm();
        } else {
            time_error_factor = 1000.;
        }
       

        /*
         * Dynamic time step control- delta_t_n, delta_t_nm1 are set here.
         */
        if (step_failed) {
            /*
             * For convergence failures, decrease the step-size by a factor of
             *  4 and try again. We store the new attempt in  delta_t_np1.
             */
            //delta_t_n_old = delta_t_n;
            delta_t_np1 = 0.25 * delta_t_n;
	    // delta_t_n = delta_t_np1;
        } else if (m_method == BEulerVarStep) {

            /*
             * If we are doing a predictor/corrector method, and we are
             * past a certain number of time steps given by the input file
             * then either correct the deltaT for the next time step based on the current time step truncation error
             * and the current integration order or chuck the step
             */
            if ((m_order > 0) && (m_time_step_num > m_numInitialConstantDeltaTSteps)) {
                /*
                 * Determine the next deltaT. If we are chucking the step, then we return a
                 * negative value of the suggested next step size here.
                 */
                delta_t_np1 = time_step_control(m_order, time_error_factor);
                if (normFilter > 0.1) {
                    if (delta_t_np1 > delta_t_n)
                        delta_t_np1 = delta_t_n;
                }
                /*
                 * Check for Current time step failing due to violation of time step
                 * truncation bounds. If negative chuck the step are redo it.
                 */
                if (delta_t_np1 < 0.0) {
                    m_numTotalTruncFails++;
                    step_failed = true;
                    delta_t_np1 = -delta_t_np1;
                    m_failure_counter += 2;
                    comment += "TIME TRUNC FAILURE";
                    rslt = "TRNC";
                }

                /*
                 * Prevent churning of the time step by not increasing the
                 * time step,
                 * if the recent "History" of the time step behavior is still bad
                 */
                else if (m_failure_counter > 0) {
                    delta_t_np1 = MIN(delta_t_np1, delta_t_n);
                }
            } else {
                delta_t_np1 = delta_t_n;
                if (!mypid_ && m_print_flag > 2) {
                    if ( (m_time_step_num <= m_numInitialConstantDeltaTSteps) && !mypid_) {
                        printf("\t  Time step error controlled ignored on %d step until step %d\n", m_time_step_num,
                                m_numInitialConstantDeltaTSteps);
                    }
                }
            }

            /* Decrease time step if a lot of Newton Iterations are
             * taken.
             * The idea being if more or less Newton iteration are taken
             * than the target number of iterations, then adjust the time
             * step downwards so that the target number of iterations or lower
             * is achieved. This
             * should prevent step failure by too many Newton iterations because
             * the time step becomes too large.  CCO
             * hkm -> put in num_new_its min of 3 because the time step
             *        was being altered even when num_newt_its == 1
             */
            int max_Newton_steps = 10000;
            int target_num_iter = 5;
            if (num_newt_its > 3000 && !step_failed) {
                if (max_Newton_steps != target_num_iter) {
                    double iter_diff = num_newt_its - target_num_iter;
                    double iter_adjust_zone = max_Newton_steps - target_num_iter;
                    double target_time_step = delta_t_n
                            * (1.0 - iter_diff * fabs(iter_diff) / ( (2.0 * iter_adjust_zone * iter_adjust_zone)));
                    target_time_step = MAX(0.5*delta_t_n, target_time_step);
                    if (!mypid_ && target_time_step < delta_t_np1) {
                        printf("\tNext time step will be decreased from %g to %g"
                                " because of new its restraint\n", delta_t_np1, target_time_step);
                        delta_t_np1 = target_time_step;
                    }
                }
            }

            /*
             * Apply any constraints given by the problem itself
             */
            double delta_t_c = m_func->delta_t_constraint(time_n, *m_y_n, *m_ydot_n);
            if (delta_t_c > 0.0) {
                if (!mypid_ && delta_t_np1 > delta_t_c) {
                    printf("\tNext time step will be decreased from %g to %g"
                            " because of the problem constraint\n", delta_t_np1, delta_t_c);
                    delta_t_np1 = delta_t_c;
                }
            }
        }

        /*
         * The final loop in the time stepping algorithm depends on whether the
         * current step was a success or not.
         */
        if (step_failed) {
            /*
             * Increment the counter indicating the number of consecutive
             * failures
             */
            numTSFailures++;
            /*
             * Print out a statement about the failure of the time step.
             */
            if (m_print_flag > 1) {
                print_time_fail(convFailure, m_time_step_num, time_n, delta_t_n, delta_t_np1, time_error_factor);
            } else if (m_print_flag == 1) {
                print_lvl1_summary(m_time_step_num, time_n, rslt, delta_t_n, num_newt_its, aztec_its, bktr_stps, time_error_factor,
                        comment.c_str());
            }

            /*
             * Change time step counters back to the previous step before
             * the failed
             * time step occurred.
             */
            time_n -= delta_t_n;
            delta_t_n = delta_t_np1;
            // delta_t_nm1 = delta_t_nm2;

            /*
             * Decide whether to bail on the whole loop
             */
            if (numTSFailures > 35) {
                giveUp = true;
            }
            /*
             * Give up on the time stepping if the time step is below a minimum. Now, the minimum
             * is set at the absolute minimum of the machine precision.
             */
            if (delta_t_np1 < 1.0E-300) {
                giveUp = true;
            }

        }

        /*
         * Do processing for a successful step.
         */
        else {

            /*
             * Decrement the number of consecutive failure counter.
             */
            m_failure_counter = MAX(0, m_failure_counter-1);

            /*
             * Print out final results of a successful time step.
             */
            if (m_print_flag > 1) {
                print_time_step2(m_time_step_num, m_order, time_n, time_error_factor, delta_t_n, delta_t_np1);
            } else if (m_print_flag == 1) {
                print_lvl1_summary(m_time_step_num, time_n, "    ", delta_t_n, num_newt_its, aztec_its, bktr_stps,
                        time_error_factor, comment.c_str());
            }

            /*
             * Output information at the end of every successful time step, if
             * requested.
             *
             * fill in
             */

        }
    } while (step_failed && !giveUp);

    /*
     *  If the time step failed we reverse the counters so that this time step routine
     *  never happened.
     */
    if (step_failed) {
        /*
         * Replace old solution vector and old time derivative solution vector.
         * into the current time step soln vectors, m_y_n, and m_ydot_n.
         */
        mdpUtil::mdp_copy_dbl_1(& (*m_y_n)[0], & ( (*m_y_nm1)[0]), m_NumLcEqns);
        mdpUtil::mdp_copy_dbl_1(& (*m_ydot_n)[0], & ( (*m_ydot_nm1)[0]), m_NumLcEqns);

        delta_t_np1 = delta_t_n;
        delta_t_n = delta_t_nm1;
        delta_t_nm1 = delta_t_nm2;
        time_n = time_nm1;
        time_nm1 = time_nm2;
        m_time_step_num--;

        if (time_n == 0.0) {
            return -1234.0;
        }
        return -time_n;
    }
    /*
     * Send back the overall result of the time step.
     */
    return time_n;
}
//====================================================================================================================

int BEulerInt::check_predicted_soln(m1d::Epetra_Vector_Ghosted & y_n, m1d::Epetra_Vector_Ghosted & ydot_n,
				     double CJ, double time_n)
{
  return 0;
}
//====================================================================================================================
// Solve for the consistent initial conditions and consistent initial time derivatives.
/*
 *  A special nonlinear problem is solved to find the consistent initial time derivatives and DAE
 *  values that cause the residual system to be solved at t = time_n.
 *
 *  This algorithm is described in the notes and utilizes the m_isAlgebraic[] field variables.
 *
 *  Given values for the non-algabraic unknowns this routine seeks to calculate a consistent value of
 *  the algebraic unknowns and ydot for the algabraic unknowns to solve the residual equations.
 *  On a successful solution, the answer is storred in the vector variables m_y_n and m_ydot_n.
 *  Note this answer depends on delta_t_n weakly. The solution depends on a integrated source and the
 *  integrated source depends on delta_t_n.
 *
 *  Note that it can be called at any time in the calculation. However, it should always be called at the
 *  start of the calculation.
 *
 *  @return 0  Always returns 0
 */
int BEulerInt::calcConsistentInitialDerivs()
{

    /*
     *  Set the max newton iterations high, as this is a do or die calculation
     */
    m_nonlin->setMaxNewtIts(150);

    m_func->setAtolVector_DAEInit(1.0E-3, *m_y_n, *m_ydot_n);

    m_func->setAtolDeltaDamping_DAEInit(1.0, *m_y_n, *m_ydot_n);

    /*
     *  Change the absolute error tolerances to that of the DAE init tolerances. Note, the time derivatives will
     *  probably have a different scaling based on the expected time response of the system.
     *
     *  -> go get atolDAE vector
     */
    const Epetra_Vector_Owned & atolDAEInitRef = m_func->atolVector_DAEInit();
    /*
     *  Tell nonlinear solver to use atolDAE
     */
    m_nonlin->setTolerances(m_reltol, m_NumLcEqns, & (atolDAEInitRef[0]));

    // not needed ???
    m_nonlin->setPredicted_soln(*m_y_n);

    // not needed ??? -> needed to specify delta_t_n
    m_nonlin->setPreviousTimeStep(delta_t_n, *m_y_nm1, *m_ydot_nm1);

    const Epetra_Vector_Owned &abs_dd = m_func->atolDeltaDamping();
    m_nonlin->setTolerances_deltaDamping(m_reltol, m_NumLcOwnedEqns, & (abs_dd[0]));
    /*
     * Solve the system of equations at the current time step.
     * Note - x_corr_n and x_dot_n are considered to be updated,
     * on return from this solution.
     */
    int num_linear_solves;
    int numbacktracks;
    int num_newt_its;
    m1d::Solve_Type_Enum ss = DAESystemInitial_Solve;
    double CJ = 1.0 / delta_t_n;

    int ierror = m_nonlin->solve_nonlinear_problem(ss, m_y_n, m_ydot_n, CJ, time_n, num_newt_its, num_linear_solves, numbacktracks);
    /*
     *  If we have experienced an error in the DAE initial solve calculation, we currently have no recourse but
     *  to end the calculation in failure. We may change this behavior in the future, I don't know.
     */
    if (ierror < 0) {
        throw CanteraError("BEulerInt::calcConsistentInitialDerivs", "Nonlinear solver failed to converge");
    }

    /*
     *   Always call writeSolution to write out the initial conditions
     */
    if (m_print_flag > 0) {
        m_func->writeSolution(0, true, time_n, delta_t_n, 0.0, *m_y_n, m_ydot_n,
                              TimeDependentAccurate_Solve, delta_t_np1);
    }
    if (m_print_flag > 0) {
        std::string snn = "Solution Time Derivative";
        m_func->showSolutionVector(snn, time_n, delta_t_n, *m_ydot_n_owned);
    }

    /*
     * Call a different user routine at the end of each step,
     * that will probably print to a file.
     */
    m_func->user_out(0, 0.0, 0.0, 0, *m_y_n, m_ydot_n);

    /*
     *  Here we show that we have in fact solved the residual equation for delta_t well. There may be some variation
     *  in the residual due to a change in the equation set. However, barring that, the residual should be near zero.
     */
    if (m_print_flag > 5) {
        double rdelta_t = 1.0 / delta_t_n;
        m_func->residEval(m_resid, true, m_y_n, m_ydot_n, 0.0, rdelta_t, Base_ResidEval, DAESystemInitial_Solve);
        double rnorm = res_error_norm(*m_resid, "DAE_Init_Residual", 10);
        printf("rnorm (DAE) = %g\n", rnorm);
    }

    return 0;
}
//====================================================================================================================
}// End of beuler namespace
//====================================================================================================================
