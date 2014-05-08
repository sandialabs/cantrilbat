/**
 *  @file m1d_SolGlobalNonlinear.h
 */

/*
 *  $Id: m1d_SolNonlinear_CurrentSolve.h 534 2013-02-22 21:33:41Z hkmoffa $
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */
#ifndef M1D_SOLGLOBALNONLINEAR_CURRENTSOLVE_H
#define M1D_SOLGLOBALNONLINEAR_CURRENTSOLVE_H

#include "m1d_SolGlobalNonlinear.h"
#include "m1d_SurDomain1D.h"

namespace m1d {

class BoundaryCondition;
class BatteryResidEval;
class SolNonlinear;

/**
 *  Wrapper class for 'beuler' integrator
 *  We derive the class from the class Integrator
 *  Newton iterator for multi-domain, one-dimensional problems.
 *  Used by class OneDim.
 */
class SolNonlinear_CurrentSolve: public SolGlobalNonlinear
{

public:
    /**
     * The default constructor doesn't take an argument.
     */
    SolNonlinear_CurrentSolve();

    //! Destructor
    virtual
    ~SolNonlinear_CurrentSolve();

    //! Setup the problem for solution.
    /*!
     *   Here we change the underlying solution to a constant voltage representation, in which we will create an
     *   outer loop to converge on a desired current.
     *   In this routine, we store the address of the jacobian and the residual function.
     */
    virtual void
    setup_problem(Solve_Type_Enum solnType, Epetra_Vector_Ghosted* y_init, Epetra_Vector_Ghosted* ydot_init, double time_curr,
                  ProblemResidEval &problem, EpetraJac& jac);

    //! Set the value of the maximum # of newton iterations
    /*!
     *  @param maxNewtIts   Maximum number of newton iterations
     *                      The default value of this is 50 iterations
     */
    virtual void setMaxNewtIts(const int maxNewtIts);

    //! Set nonlinear options in the underlying solver.
    /*!
     *  This is a pass-through function
     */
    virtual void
    setNonLinOptions(int min_newt_its, bool matrixConditioning, bool colScaling, bool rowScaling, int colScaleUpdateFrequency);

    //!  Set the level of printing that occurs during the nonlinear solve
    /*!
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
    void setPrintFlag(int print_flag);

    virtual void
    setPredicted_soln(const Epetra_Vector &y_pred);

    //! Set the absolute tolerances for the solution variables
    /*!
     *   Set the absolute tolerances used in the calculation
     *
     *  @param reltol   relative tolerance used in the nonlinear solver
     *  @param n        Length of abstol. Should be equal to m_NumLcEqns
     *  @param abstol   Vector of length n that contains the tolerances to be used for the solution variables
     */
    virtual void setTolerances(double reltol, int n, const double * const abstol);

    //! Set the absolute tolerances for the solution variables
    /*!
     *   Set the absolute tolerances used in the calculation
     *
     *  @param reltol   relative tolerance used in the nonlinear solver
     *  @param n        Length of abstol. Should be equal to m_NumLcEqns
     *  @param abstol   Vector of length n that contains the tolerances to be used for the solution variables
     */
    virtual void setTolerances_deltaDamping(double reltol_dd, int n, const double * const abstol_dd);

    //! get the residual
    /*!
     *
     * @param time_curr
     * @param rdelta_t
     * @param solnBase_ptr
     * @param solnDotBase_ptr
     */
    virtual void
    get_res(const double time_curr, const double rdelta_t, const Epetra_Vector_Ghosted *solnBase_ptr,
            const Epetra_Vector_Ghosted *solnDotBase_ptr);

public:

    //! Set the values for the previous time step
    /*!
     *   We set the values for the previous time step here. These are used in the nonlinear
     *   solve because they affect the calculation of ydot.
     *
     * @param timeStep  Time step between current time and previous time
     * @param y_nm1     Value of the solution vector at the previous time step
     * @param ydot_nm1  Value of the solution vector derivative at the previous time step
     */
    virtual void
    setPreviousTimeStep(const double timeStep, const Epetra_Vector& y_nm1, const Epetra_Vector& ydot_nm1);

    //! Main routine to launch a nonlinear solve at the current conditions
    //! whether it's a steady state solve or a time-dependent run.
    virtual int
    solve_nonlinear_problem(Solve_Type_Enum solnType, Epetra_Vector_Ghosted* y_comm, Epetra_Vector_Ghosted* ydot_comm, double CJ,
                            double time_curr, int &num_newt_its, int &num_linear_solves, int &num_backtracks);

    void
    transform_cc_to_cv(Epetra_Vector_Ghosted* y_init, Epetra_Vector_Ghosted* ydot_init, double time_curr);
    void
    transform_cv_to_cc(Epetra_Vector_Ghosted* soln, Epetra_Vector_Ghosted* solnDot, double time_curr);

    //---------------------------------------------------------------------------------------------------------
    // Member 
    //--------------------------------------------------------------------------------------------------------

    //! Pointer to the normal solver using a damped Newton's method
    m1d::SolNonlinear *m_solverConstantVoltage;

    //! Residual problem for the constant voltage or constant current battery
    //ProblemResidEval *m_func;
    BatteryResidEval *m_func;

    //! Flag for how the solver works
    /*!
     *  0   The solver changes the underlying problem to a constant voltage and then
     *      always does the root solver on top.
     *  1   The solver changes the underlying problem each nonlinear iteration. then
     *      before returning changes the problem back again.
     *  2   Solver first tries to solve the problem using constant current. If it has
     *      problems, then it solves it using the root finder capability using an
     *      underlying constant voltage representation.
     */
    int methodForSoln_;

    //! This is the actual value of the current needed at the current time
    /*!
     *  units: amps / m2
     */
    doublereal currentNeeded_;

    //! 
    doublereal currentActual_;

    //! Value of the cathode voltage at the last time step
    doublereal cathodeVoltageOld_;

    double timeOld_;

    //! Value of the last time from the last time step
    double timeLast_;

    double timeStep_;

    //! toggle indicating how the underlying problem is currently formulated.
    bool BC_Now_Voltage_;

    //! This is the type of constant current boundary condition that this object will solve
    /*!
     *   right now this is type 1, 3, 5, 7, or 9
     */
    int BC_Type_currentBC_;

    //! This is the value of the double used in the specification of the Cathode's current boundary condition
    /*!
     *    For some boundary conditions this isn't sufficient to describe the complexity of the boundary condition
     *    units amps/m2
     */
    double Value_currentBC_;

    //! Shallow pointer to the boundary condition function needed for the current bc
    /*!
     * We don't own this
     */
    m1d::BoundaryCondition* BCFuncPtr_currentBC_;

    //! Shallow pointer to the boundary condition time function needed for the current bc
    /*!
     * We don't own this
     */
    m1d::TimeDepFunctionPtr TimeDepFuncPtr_currentBC_;

    //! This is the type of constant voltage boundary condition that this object will solve
    /*!
     *   right now this is type 0, 2, 4, 6, or 8
     */
    int BC_Type_voltageBC_;


    //! This is the value of the double used in the specification of the Cathode's voltage boundary condition
    /*!
     *    For some boundary conditions this isn't sufficient to describe the complexity of the boundary condition
     *    units amps/m2
     */
    double Value_voltageBC_;

    //! Shallow pointer to the boundary condition function needed for the voltage bc
    /*!
     * We don't own this
     */
    m1d::BoundaryCondition* BCFuncPtr_voltageBC_;

    //! Shallow pointer to the boundary condition time function needed for the voltage bc
    /*!
     * We don't own this
     */
    m1d::TimeDepFunctionPtr TimeDepFuncPtr_voltageBC_;

    //! Value of the cathode voltage at the current time step
    doublereal cathodeVoltageBest_;

    //! Value of the derivative of the cathode voltage
    doublereal cathodeVoltageDot_;

    //! Boundary condition function
    //BoundaryCondition * BC_TimeDep_;

    //! Time dependency function
    //TimeDepFunctionPtr TimeDep_;

    //! Minimum Value of the cathode voltage to be attempted at the current step
    doublereal cathodeVoltageMin_;

    //! Maxium Value of the cathode voltage to be attempted at the current step
    doublereal cathodeVoltageMax_;

    //! Maximum Delta value of the cathode voltage to be taken during the root solver
    doublereal cathodeVoltageDelta_;

    //! Print level
    int printLvl_;

    doublereal rtol_;

};

} // namespace

#endif
