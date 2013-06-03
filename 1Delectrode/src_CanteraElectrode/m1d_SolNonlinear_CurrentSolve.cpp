/**
 *
 *  @file m1d_SolNonlinear.cpp
 *
 *  Damped Newton solver for 1D multi-domain problems
 */

/*
 *  $Author: hkmoffa $
 *  $Date: 2013-02-27 15:18:26 -0700 (Wed, 27 Feb 2013) $
 *  $Revision: 540 $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_defs.h"

#include "m1d_Comm.h"
#include "m1d_SolNonlinear_CurrentSolve.h"
#include "cantera/base/clockWC.h"
#include "cantera/base/mdp_allo.h"

#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_DataAccess.h"

#include "m1d_LocalNodeIndices.h"

#include "m1d_ProblemStatementCell.h"
#include "m1d_BatteryResidEval.h"

#include "m1d_DomainLayout.h"

// place where the current function is defined
#include "FuncElectrodeCurrent.h"

#include <stdio.h>
#include <math.h>

using namespace std;

#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y)    (( (x) < (y) ) ? (x) : (y))
#endif

namespace m1d {

//=====================================================================================================================
SolNonlinear_CurrentSolve::SolNonlinear_CurrentSolve() :
        SolGlobalNonlinear(),
        m_solverConstantVoltage(0),
        m_func(0),
        methodForSoln_(0),
        currentNeeded_(0.0),
        currentActual_(0.0),
        cathodeVoltageOld_(0.0),
        timeOld_(0.0),
        timeLast_(0.0),
        timeStep_(0.0),
        BC_Now_Voltage_(false),
        BC_Type_currentBC_(1),
        Value_currentBC_(0.0),
        BCFuncPtr_currentBC_(0),
        TimeDepFuncPtr_currentBC_(0),
        BC_Type_voltageBC_(0),
        Value_voltageBC_(0.0),
        BCFuncPtr_voltageBC_(0),
        TimeDepFuncPtr_voltageBC_(0),
        cathodeVoltageBest_(0.0),
        cathodeVoltageDot_(0.0),
        cathodeVoltageMin_(-10.),
        cathodeVoltageMax_(10.),
        cathodeVoltageDelta_(0.02),
        printLvl_(0),
        rtol_(1.0E-3)
{
    m_solverConstantVoltage = new SolNonlinear();
}
//=====================================================================================================================
SolNonlinear_CurrentSolve::~SolNonlinear_CurrentSolve()
{
    delete m_solverConstantVoltage;
    m_solverConstantVoltage = 0;
}
//=====================================================================================================================
void SolNonlinear_CurrentSolve::get_res(const double time_curr, const double rdelta_t, const Epetra_Vector_Ghosted *solnBase_ptr,
                                        const Epetra_Vector_Ghosted *solnDotBase_ptr)
{

}
//=====================================================================================================================
// Set the absolute tolerances for the solution variables
/*
 *   Set the absolute tolerances used in the calculation
 *
 *  @param reltol   relative tolerance used in the nonlinear solver
 *  @param n        Length of abstol. Should be equal to m_NumLcEqns
 *  @param abstol   Vector of length n that contains the tolerances to be used for the solution variables
 */
void SolNonlinear_CurrentSolve::setTolerances(double reltol, int n, const double * const abstol)
{
    rtol_ = reltol;
    m_solverConstantVoltage->setTolerances(0.5 * reltol, n, abstol);
}
//=====================================================================================================================
void SolNonlinear_CurrentSolve::setTolerances_deltaDamping(double reltol_dd, int n, const double * const abstol_dd)
{
    m_solverConstantVoltage->setTolerances_deltaDamping(reltol_dd, n, abstol_dd);
}
//=====================================================================================================================
// Setup the problem for solution
/*
 *   Here we change the underlying solution to a constant voltage representation, in which we will create an
 *   outer loop to converge on a desired current.
 *   In this routine, we store the address of the jacobian and the residual function.
 */
void SolNonlinear_CurrentSolve::setup_problem(Solve_Type_Enum solnType, Epetra_Vector_Ghosted* y_init,
                                              Epetra_Vector_Ghosted* ydot_init, double time_curr, ProblemResidEval &cv_problem,
                                              EpetraJac& jac)
{
    /*
     *  Store the address of the constant voltage or constant current problem
     */
    m_func = &cv_problem;
    BatteryResidEval *batResid = dynamic_cast<BatteryResidEval *>(m_func);

    /*
     *  Call the underlying nonlinear solver for the constant voltage problem.
     */
    m_solverConstantVoltage->setup_problem(solnType, y_init, ydot_init, time_curr, cv_problem, jac);

    /*
     * Go get the current boundary conditions. Store the results in the class variables
     *           BC_Type_currentBC_, value_currentBC_, BC_TimeDep_, TimeDep_
     */
    int BC_Type;
    double value;
    BoundaryCondition *bcFuncPtr;
    TimeDepFunctionPtr tFuncPtr;
    batResid->reportCathodeVoltageBC(time_curr, BC_Type, value, bcFuncPtr, tFuncPtr);

    if (BC_Type == 0 || BC_Type == 2 || BC_Type == 4 || BC_Type == 6 || BC_Type == 8) {
        BC_Now_Voltage_ = true;
        BC_Type_voltageBC_ = BC_Type;
        Value_voltageBC_ = value;
        BCFuncPtr_voltageBC_ = bcFuncPtr;
        TimeDepFuncPtr_voltageBC_ = tFuncPtr;
    } else {
        // 1 , 3, 5, 7, 9
        BC_Now_Voltage_ = false;
        BC_Type_currentBC_ = BC_Type;
        Value_currentBC_ = value;
        BCFuncPtr_currentBC_ = bcFuncPtr;
        TimeDepFuncPtr_currentBC_ = tFuncPtr;
    }

    /*
     *  Save some old time steps
     */
    timeOld_ = time_curr;
    timeLast_ = time_curr;

    /*
     *  Now we have to massage the problem to take out the constant current boundary condition and put in the
     *  constant voltage boundary condition
     */

    if (methodForSoln_ == 0) {
        if (BC_Now_Voltage_ == false) {
            transform_cc_to_cv(y_init, ydot_init, time_curr);
        }
    }

    // printLvl_ = 9;
}
//====================================================================================================================
// Change the problem to a constant voltage boundary condition problem
void SolNonlinear_CurrentSolve::transform_cc_to_cv(Epetra_Vector_Ghosted* soln, Epetra_Vector_Ghosted* solnDot, double time_curr)
{

    VarType v1(Voltage, 2, "CathodeVoltage");

    BatteryResidEval *batResid = dynamic_cast<BatteryResidEval *>(m_func);

    // Go get the boundary condition that is being implemented now
    int BC_Type;
    double value;
    BoundaryCondition *bcFuncPtr;
    TimeDepFunctionPtr tFuncPtr;
    batResid->reportCathodeVoltageBC(time_curr, BC_Type, value, bcFuncPtr, tFuncPtr);
    if (BC_Type == 1 || BC_Type == 3 || BC_Type == 5 || BC_Type == 7 || BC_Type == 9) {
        BC_Type_currentBC_ = BC_Type;
        Value_currentBC_ = value;
        BCFuncPtr_currentBC_ = bcFuncPtr;
        TimeDepFuncPtr_currentBC_ = tFuncPtr;
    } else {
        // 0, 2, 4, 6, 8
        printf("WARNING: SolNonlinear_CurrentSolve::transform_cc_to_cv already cv. CONFUSED!\n");
        exit(-1);
    }

    // Find the current that's needed at the current time step
    switch (BC_Type_currentBC_) {
    case 1:
        currentNeeded_ = Value_currentBC_;
        break;
    case 3:
        currentNeeded_ = Value_currentBC_ * TimeDepFuncPtr_currentBC_(time_curr);
        break;
    case 5:
    case 7:
    case 9:
        currentNeeded_ = BCFuncPtr_currentBC_->value(time_curr);
        break;
    default:
        m1d_Error("", "confused");
        break;
    }

    // Go get the current value of the cathode voltage and use that as the initial value of the voltage boundary condition
    cathodeVoltageOld_ = batResid->reportCathodeVoltage();

    // Go get the current value of the cathode current
    currentActual_ = batResid->reportCathodeCurrent();

    DomainLayout &DL = * (batResid->DL_ptr_);
    // want last surface, but be careful when we go to double tap batteries
    SurDomain1D *d_ptr = DL.SurDomain1D_List.back();
    double vCheck = d_ptr->extractSolnValue(soln, v1);
    cathodeVoltageBest_ = vCheck;
    //  if (fabs(cathodeVoltage - vCheck) > 1.0E-8) {
    // printf("we have prob\n");
    //cathodeVoltageBest_ = vCheck;
    // }

    // Store the voltage derivative -> not sure if this is used.
    cathodeVoltageDot_ = d_ptr->extractSolnValue(solnDot, v1);

    // Change the problem to a Dirichlet condition problem. We initially specify the voltage as being equal
    // to the current voltage. However, we will iterate on the voltage.
    batResid->changeCathodeVoltageBC(0, cathodeVoltageBest_);

    // Set max and min values of the voltage. These values will be used in the root-finder.
    cathodeVoltageMin_ = cathodeVoltageBest_ - 0.2;
    cathodeVoltageMax_ = cathodeVoltageBest_ + 0.2;
    cathodeVoltageDelta_ = 0.007;

    // Set the toggle for the boundary condition
    BC_Now_Voltage_ = true;
}
//====================================================================================================================
// Change the problem to a constant current boundary condition problem
void SolNonlinear_CurrentSolve::transform_cv_to_cc(Epetra_Vector_Ghosted * soln, Epetra_Vector_Ghosted * solnDot, double time_curr)
{

    VarType v1(Voltage, 2, "CathodeVoltage");

    BatteryResidEval *batResid = dynamic_cast<BatteryResidEval *>(m_func);

    // Go get the boundary condition that is being implemented now
    int BC_Type;
    double value;
    BoundaryCondition *bcFuncPtr;
    TimeDepFunctionPtr tFuncPtr;
    batResid->reportCathodeVoltageBC(time_curr, BC_Type, value, bcFuncPtr, tFuncPtr);
    /*
     *  First determine if the boundary condition is a Dirichlet condition
     *  If it is, then store the current boundary condition before it is overwritten.
     */
    if (BC_Type == 0 || BC_Type == 2 || BC_Type == 4 || BC_Type == 6 || BC_Type == 8) {
        //BC_Type_voltageBC_ = BC_Type;
        //Value_voltageBC_ = value;
        //BCFuncPtr_voltageBC_ = bcFuncPtr;
        //TimeDepFuncPtr_voltageBC_ = tFuncPtr;
    } else {
        // 1 , 3, 5, 7, 9
        printf("WARNING: SolNonlinear_CurrentSolve::transform_cv_to_cc already cc.\n");
        return;
    }

    if (BC_Type_currentBC_ != 0) {
        printf("we have prob\n");
        exit(-1);
    }

    // go get the current value of the cathode current
    currentActual_ = batResid->reportCathodeCurrent();

    // Go get the current value of the cathode voltage and use that as the initial value of the voltage boundary condition
    double cathodeVoltage = batResid->reportCathodeVoltage();

    DomainLayout &DL = * (batResid->DL_ptr_);

    // We want the last surface, but be careful when we go to double tap batteries
    SurDomain1D *d_ptr = DL.SurDomain1D_List.back();
    double vCheck = d_ptr->extractSolnValue(soln, v1);
    if (fabs(cathodeVoltage - vCheck) > 1.0E-8) {
        printf("we have prob, cv_to_cc\n");
        exit(-1);
    }

    cathodeVoltageOld_ = cathodeVoltage;
    cathodeVoltageBest_ = cathodeVoltage;
    cathodeVoltageDot_ = d_ptr->extractSolnValue(solnDot, v1);

    // Change the problem to a flux boundary condition
    batResid->changeCathodeVoltageBC(BC_Type_currentBC_, Value_currentBC_, BCFuncPtr_currentBC_, TimeDepFuncPtr_currentBC_);

    cathodeVoltageMin_ = cathodeVoltage - 0.1;
    cathodeVoltageMax_ = cathodeVoltage + 0.1;
    cathodeVoltageDelta_ = 0.007;
    BC_Now_Voltage_ = false;
}
//=====================================================================================================================
// Set the value of the maximum # of newton iterations
/*
 *  @param maxNewtIts   Maximum number of newton iterations
 *                      The default value of this is 50 iterations
 */
void SolNonlinear_CurrentSolve::setMaxNewtIts(const int maxNewtIts)
{
    m_solverConstantVoltage->setMaxNewtIts(maxNewtIts);
}
//=====================================================================================================================
void SolNonlinear_CurrentSolve::setNonLinOptions(int min_newt_its, bool matrixConditioning, bool colScaling, bool rowScaling,
                                                 int colScaleUpdateFrequency)
{
    m_solverConstantVoltage->setNonLinOptions(min_newt_its, matrixConditioning, colScaling, rowScaling, colScaleUpdateFrequency);
}
//=====================================================================================================================
//  Set the level of printing that occurs during the nonlinear solve
/*
 *   0 -> absolutely nothing is printed for a single time step.
 *   1 -> One line summary per time step
 *   2 -> short description, points of interest
 *   3 -> More printed per time step -> major algorithm issues are displayed
 *   4 -> Additional time step error control information is printed out
 *        One line summary of the nonlinear solve
 *   5 -> Summaries of the nonlinear solve iterates are printed out
 *   6 -> Algorithm information on the nonlinear iterates are printed out
 *   7 -> Additional info on the nonlinear iterates are printed out
 *   8 -> Additional info on the linear solve is printed out.
 *   9 -> Info on a per iterate of the linear solve is printed out.
 */
void SolNonlinear_CurrentSolve::setPrintFlag(int print_flag)
{
    m_print_flag = print_flag;
    if (m_solverConstantVoltage) {
        if (print_flag < 2) {
            m_solverConstantVoltage->setPrintFlag(0);
        } else {
            m_solverConstantVoltage->setPrintFlag(printLvl_);
        }
    }
}
//=====================================================================================================================
/*
 *
 */
void SolNonlinear_CurrentSolve::setPredicted_soln(const Epetra_Vector &y_pred)
{
    m_solverConstantVoltage->setPredicted_soln(y_pred);
}
//=====================================================================================================================
// Set the values for the previous time step
/*
 *   We set the values for the previous time step here. These are used in the nonlinear
 *   solve because they affect the calculation of ydot.
 *
 * @param timeStep  Time step between current time and previous time
 * @param y_nm1     Value of the solution vector at the previous time step
 * @param ydot_nm1  Value of the solution vector derivative at the previous time step
 */
void SolNonlinear_CurrentSolve::setPreviousTimeStep(const double timeStep, const Epetra_Vector& y_nm1,
                                                    const Epetra_Vector& ydot_nm1)
{
    m_solverConstantVoltage->setPreviousTimeStep(timeStep, y_nm1, ydot_nm1);

    timeStep_ = timeStep;
    // We may have to extract the old cathode voltage here.
}
//=====================================================================================================================
/*
 * solve_nonlinear_problem():
 *
 * Find the solution to F(X) = 0 by damped Newton iteration.  On
 * entry, x0 contains an initial estimate of the solution.  On
 * successful return, x1 contains the converged solution.
 *
 * SolnType = TRANSIENT -> we will assume we are relaxing a transient
 *        equation system for now. Will make it more general later,
 *        if an application comes up.
 */
int SolNonlinear_CurrentSolve::solve_nonlinear_problem(Solve_Type_Enum solnType, Epetra_Vector_Ghosted *y_comm,
                                                       Epetra_Vector_Ghosted *ydot_comm, double CJ, double time_curr,
                                                       int &num_newt_its_comm, int &num_linear_solves, int &num_backtracks)
{

    if (time_curr > timeOld_) {
        cathodeVoltageOld_ = cathodeVoltageBest_;
        timeOld_ = timeLast_;
    }

    double deltaT = 1.0 / CJ;

    if (methodForSoln_ == 1) {
        transform_cc_to_cv(y_comm, ydot_comm, time_curr);
    }
    if (BC_Now_Voltage_ != true) {
        m1d_Error("SolNonlinear_CurrentSolve::solve_nonlinear_problem", "confused");
    }

    BatteryResidEval *batResid = dynamic_cast<BatteryResidEval *>(m_func);

    double cathodeVoltage = batResid->reportCathodeVoltage();
    cathodeVoltageMin_ = cathodeVoltage - 0.75;
    cathodeVoltageMax_ = cathodeVoltage + 0.75;

    // double currentObtained = currentNeeded;

    CurrentFunc ec(this, m_solverConstantVoltage, m_func, y_comm, ydot_comm, CJ, time_curr, deltaT);

    /*
     *  Do a first order prediction of the cathode voltage.
     */
    double delta = deltaT * cathodeVoltageDot_;
    if (delta > cathodeVoltageDelta_) {
        delta = cathodeVoltageDelta_;
    } else if (delta < -cathodeVoltageDelta_) {
        delta = -cathodeVoltageDelta_;
    }
    double xBest_ = cathodeVoltageOld_ + delta;
    cathodeVoltageBest_ = cathodeVoltageOld_;

    // Define the root finder
    RootFind rf(&ec);
    rf.setPrintLvl(m_print_flag);
    ec.printLvl_ = m_print_flag;

    rf.setTol(rtol_, 1.0E-10, 0.2 * rtol_, 1.0E-10);
    rf.setFuncIsGenerallyDecreasing(true);
    rf.setDeltaX(cathodeVoltageDelta_);

    int timeRegion = batResid->m_currentTimeRegion;

    // Find the current that's needed at the current time step.
    // We need to do this because the current may change with time.
    switch (BC_Type_currentBC_) {
    case 1:
        currentNeeded_ = Value_currentBC_;
        break;
    case 3:
        currentNeeded_ = Value_currentBC_ * TimeDepFuncPtr_currentBC_(time_curr);
        break;
    case 5:
    case 7:
    case 9:
        currentNeeded_ = BCFuncPtr_currentBC_->value(time_curr, timeRegion);
        break;
    default:
        m1d_Error("", "confused");
        break;
    }

    /*
     *  Call the root finder. The actual current found is returned in currentObtained
     */
    double currentObtained = currentNeeded_;
    int status = rf.solve(cathodeVoltageMin_, cathodeVoltageMax_, 100, currentObtained, &xBest_);

    if (status == 0) {
        cathodeVoltageBest_ = xBest_;
        cathodeVoltageDot_ = (cathodeVoltageBest_ - cathodeVoltageOld_) / deltaT;
        if (printLvl_ > 1) {
            printf("CurrentSolve(): Volts (%g amps) = %g\n", currentObtained, xBest_);
        }
    } else {
        if (printLvl_) {
            printf("CurrentSolve(): bad status = %d Volts (%g amps) = %g. "
                    "Returning from nonlinear solver with an error condition\n", status, currentObtained, xBest_);
        }
        return -1;
    }

    if (methodForSoln_ == 1) {
        transform_cv_to_cc(y_comm, ydot_comm, time_curr);
    }

    return 1;
}
//=====================================================================================================================
}

