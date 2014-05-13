/**
 *  @file BEulerInt.h
 */

/*
 *  $Id: BEulerInt.h 560 2013-03-06 23:50:04Z hkmoffa $
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */
#ifndef CT_BEULERINT_H
#define CT_BEULERINT_H

#include "cantera/numerics/Integrator.h"
#include "cantera/numerics/SquareMatrix.h"

#include "m1d_EpetraJac.h"
#include "m1d_ProblemResidEval.h"
#include "m1d_SolGlobalNonlinear.h"
#include "m1d_SolNonlinear.h"

#include "Epetra_Comm.h"
#include "Epetra_IntVector.h"

#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y)    (( (x) < (y) ) ? (x) : (y))
#endif

namespace beuler {

#define OPT_SIZE 10

#define BE_SUCCESS 0
#define BE_FAILURE 1

#define BE_STEADY 0
#define BE_TRANSIENT 1

//! Macro to check pointers before deleting their contents.
/*!
 *   This macro will set the pointer to zero, indicating that the contents are
 *   now inaccessible
 */

enum BEulerMethodType
{
    BEulerFixedStep = 0, BEulerVarStep
};

/**
 * Exception class thrown when a BEuler error is encountered.
 */
class BEulerErr: public Cantera::CanteraError
{
public:
    BEulerErr(std::string msg);
};

#define BEULER_JAC_ANAL 2
#define BEULER_JAC_NUM  1

//!  Wrapper class for 'beuler' integrator We derive the class from the class Integrator
/*!
 *    This class does a basis backwards euler predictor corrector iteration.
 *
 *   It has additional hooks.
 *
 *   Treatment of Algebraic Constraints
 *   ----------------------------------------------------
 *
 *
 *   Evaluation of Solution Norms
 *  ------------------------------------------------------
 *
 *
 *   Tolerances - What they are and how to implement them.
 *   ------------------------------------------------------
 *
 *   The absolute tolerances are important inputs for the nonlinear solver and the
 *   time integrator.  These absolute tolerances are provided by the ProblemResidualEval object.
 *
 *
 *
 *
 */
class BEulerInt: public Cantera::Integrator
{

public:

    //!The default constructor doesn't take an argument.
    /*!
     *  Default settings: epetra jacobian, no user-supplied
     *  Jacobian function, Newton iteration.
     *  */
    BEulerInt();

    //! Copy constructor
    /*!
     *
     * @param r  Object to be copied
     */
    BEulerInt(const BEulerInt &r);

    //! Virtual destructor
    virtual
    ~BEulerInt();

    //! Assignment operator
    /*!
     * @param r  Object to be copied
     * @return  Returns the current object
     */
    BEulerInt &
    operator=(const BEulerInt &r);

    /*************************************** Member functions ***********************************************/
    /******************************        BASIC SETUP FUNCTIONS        *********************************/

    virtual void
    setTolerancesEpetra(double reltol, const m1d::Epetra_Vector_Owned &abstol);

    //! Set the algebraic constraints in the system
    /*!
     *  @param algebraicConstraint
     */
    virtual void setAlgebraicConstraints(const Epetra_IntVector &algebraicConstraint);

    virtual void setProblemType(int probtype);

    virtual void setTimeRegionBoundaries(const std::vector<double>& timeRegionBoundaries);

    //! Set or reset the initial time of the integration along with the
    virtual void setTimeRegion(double t0, int iregion);

    //!  Specify the nonlinear solver to be used in the problem
    /*!
     *   Note, a default nonlinear solver is used, so this step is not necessary necessarily
     *   Note, the BEulerInt object owns the nonlinear solver, so nothing needs to be freed.
     */
    void specifyNonLinearSolver(m1d::SolGlobalNonlinear *nonlin);

    //!  Setup the problem by specifying the object that generates the problem to be solved,
    //!  mallocing internal objects and memory
    /*!
     *
     * @param func       func is the object that supplies the problem.
     */
    virtual void initializePRE(m1d::ProblemResidEval& func);

    //!  Setup the problem by specifying the object that generates the problem to be solved,
    //!  mallocing internal objects and memory
    /*!
     *
     * @param func       func is the object that supplies the problem.
     */
    virtual void reinitializePRE(m1d::ProblemResidEval& func);

    //! Default initialization routines that are now not allowed
    /*!
     * Initialize the integrator for a new problem. Call after
     * all options have been set.
     * @param t0   initial time
     * @param func RHS evaluator object for system of equations.
     */
    virtual void initialize(double t0, Cantera::FuncEval& func)
    {
        throw m1d::m1d_Error("initilize", "error");
    }

    //! Default initialization routines that are now not allowed
    /*!
     * Initialize the integrator for a new problem. Call after
     * all options have been set.
     * @param t0   initial time
     * @param func RHS evaluator object for system of equations.
     */
    virtual void reinitialize(double t0, Cantera::FuncEval& func)
    {
        throw m1d::m1d_Error("initilize", "error");
    }

    //! Call the problem residual to establish initial conditions
    /*!
     *  @param t0        Initial time for the simulation
     *  @param delta_t   Initial time step for the simulation 
     *                   (this is used if no other source for this is specified
     *
     */
    void determineInitialConditions(double t0, double delta_t);

    /*************************************** Member functions ***********************************************/
    /******************************       BASIC INTEGRATION FUNCTIONS                    *********************************/

    virtual double integratePRE(double tout);

    //!    This routine advances the calculations one time step
    /*!
     * It uses a predictor  corrector approach. We use an implicit algorithm here.
     *
     * @param t_max   Final time that the simulation should be advanced to. The simulation will
     *                pick a delta_t and and t_final. t_final will not be greater than t_max.
     *
     * @return  returns the time that the calculation was advanced to
     */
    virtual doublereal step(double t_max);

    //! Recalculate the solution weights based on the current values
    virtual void setSolnWeights();

    /*************************************** Member functions ***********************************************/
    /******************************     SOLUTION ACCESS FUNCTIONS                 *********************************/

    //!    Return a changeable reference to the solution vector
    /*!
     *  @param returns a ghosted Epetra_Vector reference to the solution vector
     */
    m1d::Epetra_Vector_Ghosted & solnVector()
    {
        return (*m_y_n);
    }

    //!    Return a changeable reference to the solution dot vector
    /*!
     *  @param returns a ghosted Epetra_Vector reference to the solution vector
     */
    m1d::Epetra_Vector_Ghosted & solnDotVector()
    {
        return (*m_ydot_n);
    }

    //! Return the number of global equations in the equation system
    /*!
     *
     * @return returns the number of global equations in the equation system
     */
    int nEquations() const
    {
        return m_NumGlEqns;
    }

    //! Return the total number of function evaluations
    virtual int
    nEvals() const;

    //! find the time region
    /*!
     * @param   val
     *  @param start   If the input is a starting point. If so, we return the greater region
     *  @return returns the region
     */
    int findTimeRegion(double val, bool start);

    /*************************************** Member functions ***********************************************/
    /******************************     PARAMETER SPECIFICATION FUNCTIONS                *********************************/

    virtual void
    setMethodBEMT(BEulerMethodType t);

    //! Set the maximum time step size allowed in the calculation
    /*!
     * @param hmax value of the maximum time step size
     */
    virtual void
    setMaxStep(doublereal hmax);

    //! Set the minimum time step size allowed in the calculation
    /*!
     * @param hmin value of the maximum time step size
     */
    virtual void
    setMinStep(doublereal hmin);

    virtual void
    setMaxNumTimeSteps(int maxTimeSteps);

    //! Set the number of constant initial delta T steps
    /*!
     *  Frequently it is advantageous to carry out a number of time steps using a constant
     *  delta T, before the algorithm switches over to a variable step method whose time
     *  step control is bassed on truncation error control.
     *  This input controls the number of steps to take at the beginning of the calculation.
     *  The default is 2.
     *
     *  @param numSteps   Number of steps to keep constant
     */
    virtual void
    setNumInitialConstantDeltaTSteps(int numSteps);

    virtual void
    setPrintSolnOptions(int printSolnStepInterval, int printSolnNumberToTout, int printSolnFirstSteps = 0, bool dumpJacobians =
                                false);
    void
    setNonLinOptions(int min_newt_its = 0, bool matrixConditioning = false, bool colScaling = false, bool rowScaling = true);

    //!  Check to see that the predicted solution satisfies proper requirements.
    /*!
     *       @return Returns a negative number if the step is inappropriate. Then the stepsize is reduced
     *               and the method is checked again.
     */
    virtual int check_predicted_soln(m1d::Epetra_Vector_Ghosted & y_n, m1d::Epetra_Vector_Ghosted & ydot_n,
				     double CJ, double time_n);

    //!  Set the level of printing that occurs for each time step
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
     *
     * @param print_flag   Flag for printout of the time-stepper
     *
     * @param nonlinPrintflag    Flag for the nonlinear solver
     *                           Default is to control nonlinear printout from the time stepper.
     */
    virtual void setPrintFlag(const int print_flag, const int nonlinPrintFlag = -1);

    virtual void setColumnScales();

    //!  Calculate the L2 Weighted Norm of a delta in the solution
    /*!
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
    virtual double
    soln_error_norm(const Epetra_Vector &delta_y, const int printLargest = 10, const int typeYsoln = 1, const double dampFactor =
                            1.0) const;

    //!    L2 Weighted Norm of the residual
    /*!
     *   The vector m_residWt[i]'s are always used to weight the solution errors in
     *   the calculation.
     *
     *   The third argument has a default of false. However,
     *   if true, then a table of the largest values is printed
     *   out to standard output.
     *
     *   @param resid       Residual vector
     *   @param title      Printed out on the title line
     *   @param printLargest if True a table is printed of the largest contributors.
     */
    virtual double
    res_error_norm(const Epetra_Vector &resid, const char *title = 0, const int printLargest = 0) const;

    //! Set the initial value of the time step
    /*!
     *  @param delta_t
     */
    virtual void
    setInitialTimeStep(double delta_t);

    //!  Function called by BEuler to evaluate the Jacobian matrix and the
    //!  current residual at the current time step.
    /*!
     *  @param J = Jacobian matrix to be filled in
     *  @param f = Right hand side. This routine returns the current
     *             value of the rhs (output), so that it does
     *             not have to be computed again.
     *  @param time_curr current time
     *  @param CJ   coefficient for the time derivative term
     *  @param y_curr    current solution
     *  @param y_curr_dot  Current solution derivative
     *  @param num_newt_its  number of newton iterations
     *  @param jacType        0  Normal jacobian
     *                        1  DAE initialization problem
     */
    void
    beuler_jac(m1d::EpetraJac &J, m1d::Epetra_Vector_Owned * const f, double time_curr, double CJ,
               m1d::Epetra_Vector_Ghosted * const y_curr, m1d::Epetra_Vector_Ghosted * const y_curr_dot, int num_newt_its,
               int jacType = 0);

    /*************************************** Member functions ***********************************************/
    /******************************     INTERNAL PROCEDURE CALLS                  *********************************/

protected:
    void
    internalMalloc();

    // Function to calculate the predicted solution vector, m_y_pred_n for the (n+1)th time step.
    /*
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
     *    y_n[]       - solution vector at time n
     *    y_dot_n[]   - acceleration vector from the predictor at time n
     *    y_dot_nm1[] - acceleration vector from the predictor at time n - 1
     *
     * on output:
     *
     *    m_y_pred_n[]    - predicted solution vector at time n + 1
     */
    void
    calc_y_pred(const int order);

    //! Bound the y_pred step.
    /*!
     *  This will involve lowering the step size
     */
    void
    boundStep_y_pred();

    /**
     * Internal function to calculate the time derivative at the
     * new step
     */
    void
    calc_ydot(const int order, m1d::Epetra_Vector_Ghosted & y_curr, m1d::Epetra_Vector_Ghosted & ydot_curr);

    //! This function calculates the time step error norm using an L2 formula
    /*!
     *    on input:
     *
     *      m_y_n        -  Actual value of the solution
     *      m_y_pred_n   -  Predicted value of the solution
     *      m_ewt        - Weight vector
     *
     *   Returns the l2 norm of (m_y_n[i] - m_y_pred_n[i]) / m_ewt[i]
     */
    double time_error_norm() const;

    //! Time step control function for the selection of the time step size
    /*!
     * The algorithm is based on a desired accuracy of time integration and on an estimate of the relative
     * error of the time integration process. This routine can be called for a first order - forward Euler/backward
     * Euler predictor/ corrector and for a  second order Adams- Bashforth/Trapezoidal Rule predictor/corrector.
     * See Nachos documentation Sand86-1816 and Gresho, Lee, Sani LLNL report UCRL -
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
    double
    time_step_control(int m_order, double time_error_factor) const;

    //! computeResidWts():
    /*!
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
    void
    computeResidWts();

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
    virtual int calcConsistentInitialDerivs();

    /*
     * Filter a new step
     */
    double
    filterNewStep(double timeCurrent, m1d::Epetra_Vector_Ghosted &y_curr, m1d::Epetra_Vector_Ghosted &ydot_current);

    /*************************************** Member functions ***********************************************/
    /******************************          PRINTING FUNCTIONS                     *********************************/

    //! Returns the next time of a printout
    /*!
     * The variable m_printSolnNumberToTout causes the program to print out the solution at even intervals
     * from the beginning of the simulation to the end. This program will return the time for the next printout.
     *   @param   time_current  Current time
     *   @return  returns the next print out time
     */
    double
    getPrintTime(const double time_current) const;

    void
    print_time_step1(int order, int n_time_step, double time, double delta_t_n, double delta_t_nm1, bool step_failed,
                     int num_failures) const;

    void
    print_time_step2(int time_step_num, int order, double time, double time_error_factor, double delta_t_n,
                     double delta_t_np1) const;

    void
    print_time_fail(bool convFailure, int time_step_num, double time, double delta_t_n, double delta_t_np1,
                    double time_error_factor) const;

    void
    print_final(double time, int step_failed, int time_step_num, int num_newt_its, int total_linear_solves, int numConvFails,
                int numTruncFails, int nfe, int nJacEval) const;

    void
    print_lvl1_Header(int nTimes) const;

    void
    print_lvl3_Header() const;

    void
    print_lvl1_summary(int time_step_num, double time, const char *rslt, double delta_t_n, int newt_its, int aztec_its,
                       int bktr_stps, double time_error_factor, const char *comment) const;

public:
    void
    print_solnDelta_norm_contrib(const double * const soln0, const char * const s0, const double * const soln1,
                                 const char * const s1, const char * const title, const double * const y0, const double * const y1,
                                 double damp, int num_entries);

protected:
    //! Print a solution-like vector in a form that is easy to interpret
    /*!
     *  Prints a vector that has the same layout as the solution vector.
     *  The variable names and node positions are printed along with the values.
     *
     *  This routine uses the member function showSolutionVector() of the
     *  ProblemResidEval class, which does the layout printing. It can be used to
     *  print to a file as well.
     *
     * @param header   Header that is printed out to
     * @param v        Epetra vector -> Note only the owned elements are needed.
     *                     (mp behavior hasn't been checked)
     */
    void print_SVector(std::string header, const Epetra_MultiVector &v) const;

    //! Print a solution-like int vector in a form that is easy to interpret
    /*!
     *  Prints an int vector that has the same layout as the solution vector.
     *  The variable names and node positions are printed along with the values.
     *
     *  This routine uses the member function showSolutionIntVector() of the
     *  ProblemResidEval class, which does the layout printing. It can be used to
     *  print to a file as well.
     *
     * @param header   Header that is printed out to
     * @param v        Epetra vector -> Note only the owned elements are needed.
     *                     (mp behavior hasn't been checked)
     */
    void print_IntSVector(std::string header, const Epetra_IntVector &v) const;

    /*************************************** Member data ***********************************************/

    //! Pointer to the nonlinear solver
    m1d::SolGlobalNonlinear *m_nonlin;

      
    const Epetra_Comm *Comm_ptr_;

    //! My processor number
    int mypid_;

    /*********************
     *   Size of the Problem
     *********************/
    //! Number of global equations in the equation system
    int m_NumGlEqns;

    //! Number of equations on the processor including ghost equations
    int m_NumLcEqns;

    //! Number of equations on the processor not including ghost equations
    int m_NumLcOwnedEqns;

    /*****************************************************************************
     *   METHOD OPTIONS FOR THE TIME STEPPER
     ****************************************************************************/

    /**
     * MethodType is used to specify how the time step is to be
     * chosen. Currently, there are two choices, one is a fixed
     * step method while the other is based on a predictor-corrector
     * algorithm and a time-step truncation error tolerance.
     */
    BEulerMethodType m_method;

public:
    //! Do a calculation to find the consistent initial conditions
    int m_doSpecialStartCalc;

protected:
    //!  m_jacFormMethod determines how a matrix is formed.
    int m_jacFormMethod;

    /**
     * m_rowScaling is a boolean. If true then row sum scaling
     * of the Jacobian matrix is carried out when solving the
     * linear systems.
     */
    bool m_rowScaling;

    /**
     * m_colScaling is a boolean. If true, then column scaling
     * is performed on each solution of the linear system.
     */
    bool m_colScaling;

    /**
     * m_matrixConditioning is a boolean. If true, then the
     * Jacobian and every rhs is multiplied by the inverse
     * of a matrix that is suppose to reduce the condition
     * number of the matrix. This is done before row scaling.
     */
    bool m_matrixConditioning;

    /**
     *  Relative time truncation error tolerances
     */
    double m_reltol;

public:
    //! Map of the ghosted unknowns of the problem
    Epetra_BlockMap *m_ghostedMap;

    //! Map of the owned unknowns of the problem
    Epetra_BlockMap *m_ownedMap;
protected:
    /**
     *  Vector of absolute time truncation error tolerance
     *  when not uniform for all variables.
     */
    m1d::Epetra_Vector_Owned *m_abstol;
    /**
     * Error Weights. This is a surprisingly important quantity.
     */
    m1d::Epetra_Vector_Owned *m_ewt;

    //! Current integration order for the time stepping algorithm
    int m_order;

    //! Maximum integration order of the time-stepping algorithm
    int m_maxord;

    //!  Maximum permissible time step
    /*!
     *   No time step is allowed to be greater than this value
     */
    double delta_t_max;

    //!  Minimum permissible time step
    /*!
     *   No time step is allowed to be smaller than this value
     */
    double delta_t_min;

    /****************************************************************************
     *    METHOD OPTIONS FOR THE NONLINEAR SOLVER
     ****************************************************************************/

    //! Minimum Number of Newton Iterations per nonlinear step
    /*!
     * default = 0
     */
    int m_min_newt_its;

    /**
     * Time step number
     */
    int m_time_step_num;
    int m_time_step_attempts;
    /**
     * Max time steps allowed
     */
    int m_max_time_step_attempts;

    /******************************************************************************
     *    METHOD OPTIONS FOR THE LINEAR SOLVER
     ********************************************************************************/
    //! Dump Jacobians to disk
    /*!
     *     default false
     */
    bool m_dumpJacobians;

    //! Toggle solution damping
    bool doResidSolnDamping_;

    //! Toggle delta damping
    bool doDeltaDamping_;

    /**
     * Number of initial time steps to take where the
     * time truncation error tolerances are not checked. Instead
     * the delta T is uniform
     */
    int m_numInitialConstantDeltaTSteps;

    /**
     * Failure Counter -> keeps track of the number
     * of consequetive failures
     */
    int m_failure_counter;

    /*******************************************************************************
     *  METHOD OPTIONS FOR PRINTING ALGORITHM BEHAVIOR AND SOLUTION PROGRESS
     ********************************************************************************/

    //! Determines the level of printing for each time step.
    /*!
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
    int m_print_flag;

    //!Step Interval at which to print out the solution
    /*!
     * default = 1;
     * If set to zero, there is no printout
     */
    int m_printSolnStepInterval;
    /**
     * Number of evenly spaced printouts of the solution
     * If zero, there is no printout from this option
     *  default 1
     * If set to zero there is no printout.
     */
    int m_printSolnNumberToTout;

    /*******************************************************************************
     *      COUNTERS FOR SOLUTION PROGRESS
     ********************************************************************************/

    /**
     * Number of initial steps that the solution is
     * printed out.
     *     default = 0
     */
    int m_printSolnFirstSteps;

    /**********************************************************************************
     * INTERNAL SOLUTION DATA VALUES
     ***********************************************************************************/

    m1d::Epetra_Vector_Ghosted *m_y_n;
    m1d::Epetra_Vector_Ghosted *m_y_nm1;
    m1d::Epetra_Vector_Ghosted *m_y_pred_n;
    m1d::Epetra_Vector_Ghosted *m_ydot_n;
    m1d::Epetra_Vector_Ghosted *m_ydot_nm1;

    m1d::Epetra_Vector_Owned *m_y_n_owned;
    m1d::Epetra_Vector_Owned *m_ydot_n_owned;

    /*****************************************************************************************
     *        INTERNAL RESIDUAL VALUES
     ******************************************************************************************/

    //! Value of the residual
    /*!
     *  Vector over the owned nodes.
     */
    m1d::Epetra_Vector_Owned *m_resid;

    /*****************************************************************************************
     *        INTERNAL WEIGHTS FOR TAKING SOLUTION NORMS
     ******************************************************************************************/

    //! Boolean indicating whether the DOF is an algebraic constraint on the system.
    Epetra_IntVector *m_isAlgebraic;

    /*****************************************************************************************
     *        INTERNAL BOUNDARY INFO FOR SOLUTIONS
     *****************************************************************************************/

    /*****************************************************************************************
     *        INTERNAL TIME VARIABLES
     ******************************************************************************************/

    //! Initial time at the start of the integration
    double m_t0;
    /**
     * Final time
     */
    double m_time_final;
    /**
     *
     */
    double time_n;
    double time_nm1;
    double time_nm2;
    double delta_t_n;
    double delta_t_nm1;
    double delta_t_nm2;
    double delta_t_np1;

    //! value of the residual weights
    m1d::Epetra_Vector_Owned *m_residWts;
    m1d::Epetra_Vector_Owned *m_wksp;

    //! Main hook into the problem
    /*!
     *   Amongst other things this function calculates the residual.
     */
    m1d::ProblemResidEval *m_func;

    m1d::Epetra_Vector_Owned *m_rowScales;
    m1d::Epetra_Vector_Owned *m_colScales;

    //!  Pointer to the Jacobian representing the time dependent problem
    m1d::EpetraJac *tdjac_ptr;

    /**
     * Number of function evaluations
     */
    int m_nfe;
    /**
     * Number of Jacobian Evaluations and
     * factorization steps (they are the same)
     */
    int m_nJacEval;
    /**
     * Number of total newton iterations
     */
    int m_numTotalNewtIts;
    /**
     * Total number of linear iterations
     */
    int m_numTotalLinearSolves;
    /**
     * Total number of convergence failures.
     */
    int m_numTotalConvFails;
    /**
     * Total Number of time truncation error failures
     */
    int m_numTotalTruncFails;
    /*
     *
     */
    int num_failures;

    int m_DAEDOT_Predicted_Bad;
    int m_DOT_Predicted_Bad;
    int m_DAE_Predicted_Bad;

    Epetra_IntVector *m_isArithmeticScaled;

    //!  Number of time regions defined in the problem
    /*!
     *  A time region is a region of time where the boundary conditions are specified. Between
     *  time regions there can be step discontinuities in the boundary conditions.
     */
    int m_numTimeRegions;

    //! Current time region defined in the problem
    /*!
     *  This number varies from zero to  m_numTimeRegions.
     *  The basic idea is that step jumps in boundary conditions can occur only at time region
     *  boundaries.
     */
    int m_currentTimeRegion;

    //! Value of the time region boundaries
    /*!
     *    The code operates in two modes. If this is size(0), the m_currentTimeRegion is used
     *    to evaluate what time region the integrator is in. No information is needed as to what
     *    region or what times pertain to each region is needed in the code.
     *    If size() > 0, then we have specified the regions ahead of time. All integrations are
     *    done with these times in mind. And no step is allowed to take place across time region
     *    boundaries.
     */
    std::vector<double> m_timeRegionBoundaries;

    /***********************************************************************************************
     *   INTERNAL COMMUNICATORS
     **********************************************************************************************/
};

} // namespace

#endif // CT_BEULER
