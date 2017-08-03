/**
 *  @file m1d_SolGlobalNonlinear.h
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef M1D_SOLGLOBALNONLINEAR_H
#define M1D_SOLGLOBALNONLINEAR_H

#include "m1d_exception.h"
#include "m1d_ProblemResidEval.h"

#include "Epetra_Vector.h"
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! This is the nonlinear solver base class
/*!
 *   This is the main nonlinear solver for steady state and time-dependent problems.
 *   It employs a two-way convergence criteria: one for the residual and one for the solution update vector
 *   A base class is needed because the equation sometimes need a root solver based on a constant voltage solve to solve for
 *   a given current condition.
 */
class SolGlobalNonlinear
{
public:
    //! The default constructor doesn't take an argument.
    SolGlobalNonlinear();

    //! destructor
    virtual ~SolGlobalNonlinear();

    //! Setup the problem for solution.
    /*!
     *   Find the solution to F(X) = 0 by damped Newton iteration.  On entry, x0 contains an initial estimate of the solution. 
     *   on successful return, x1 contains the converged solution.
     *
     *  @param[in]           solveType           Enum of type Solve_Type. Describes the type of the problem to be solved.
     *                                               SteadyState_Solve 
     *                                               TimeDependentAccurate_Solve
     *                                               TimeDependentInitial
     *                                               TimeDependentRelax_Solve,
     *                                               DAESystemInitial_Solve
     *                                           This parameters needs to be passed down to residual calculation, and in some 
     *                                           cases, this effects the nonlinear solve itself.
     *  @param[in]           y_init              Pointer to a ghosted Epetra_Vector. This contains the initial solution at
     *                                           the current time.
     *  @param[in]           ydot_init           Pointer to a ghosted Epetra_Vector. This may contain the initial ydot value
     *                                           at the current time. It may be nullptr for steady state problems.
     *  @param[in]           time_curr           Current time
     *  @param[in]           problem             Reference to the  ProblemResidEval object
     *  @param[in]           jac                 Reference to teh EpetraJac object
     */
    virtual void 
    setup_problem(Zuzax::Solve_Type solveType, const Epetra_Vector_Ghosted* const y_init,
                  const Epetra_Vector_Ghosted* const ydot_init, double time_curr, ProblemResidEval& problem,
                  EpetraJac& jac);

    //! Apply hard bounds on the step size due to bounds constraints
    /*!
     *  (virtual from SolGlobalNonlinear)
     *
     *  We apply this before other terms, by decreasing the size of the original step size.
     *  There are two different types of bounds constraints that are implemented here:
     *
     *        highlow bounds constraints: The values of the unknowns can be kept within a strict bounds. Note, the
     *                                    equations have to have this as a principle. i.e., there must be a 
     *                                    maximum principle.
     *        delta bounds constraints:   The unknowns shouldn't change by an overall amount on any one step.
     *                                    This stems from the trust region concept. The nonlinear jacobian isn't 
     *                                    predictive beyond a certain trust region.
     *
     *  @param[in]           y_old               Reference to Epetra ghosted Value of the solution at the previous nonlinear step.
     *  @param[in,out]       step                Reference to the Epetra owned-only step change of the predicted solution change.
     *                                           On return step is scaled by the fbound factor.
     *
     *  @param[out]          fbound              Factor that the step size has to be reduced by.
     *
     *  @return                                  Returns the following values:
     *                                             1 :  Returns 1 if there is no reduction in the step size
     *                                             0 :  Returns 0 if there is a reduction in the step size 
     *                                                  due to bounds constraints
     *                                            -3 :  fbounds has become too small. Signal that the newton algorithm has failed.
     */
    virtual int doHardBounds(const Epetra_Vector_Ghosted& y_old, Epetra_Vector_Owned& step, double& fbound);

    //! Compute factor to keep all components in bounds.
    /*!
     *  (virtual from SolGlobalNonlinear)
     *
     *  Return the factor by which the undamped Newton step 'step0'
     *  must be multiplied in order to keep all solution components in
     *  all domains not-appreciably changing so much that the jacobian isn't representative.
     *
     *  Delta bounds: The idea behind these is that the Jacobian couldn't possibly be representative, if the
     *               variable is changed by a lot. (true for nonlinear systems, false for linear systems)
     *
     *    For variables which have a strict minimum of zero: 
     *      Maximum increase in variable in any one newton iteration: 
     *          a)   factor of 2 when above the value of the fabs(change) is above m_ewt_deltaDamping[i]
     *          b)   Equal to the change if the fabs(change) is below m_ewt_deltaDamping[i].
     *
     *      Maximum decrease in variable in any one newton iteration:
     *          a)   factor of 5 when above the value of the fabs(change) is above m_ewt_deltaDamping[i]
     *          b)   Equal to the change if the fabs(change) is below m_ewt_deltaDamping[i].
     *
     *    For arithmetically scaled variables, the maximum increase or decrease in an iteration is given by the value of 
     *    m_ewt_deltaDamping[i].
     *
     *  @param[in]           y                   Ghosted Epetra_Vector reference for the current solution 
     *  @param[in]           step0               Owned Epetra_Vector reference for the current step in the solution unknowns
     *
     *  @return                                  Returns the factor which the step should be reduced by.
     */
    virtual double deltaBoundStep(const Epetra_Vector_Ghosted& y, const Epetra_Vector_Owned& step0);

    //! Compute factor to keep all components in bounds.
    /*!
     *  (virtual from SolGlobalNonlinear)
     *
     *  Return the factor by which the undamped Newton step 'step0'
     *  must be multiplied in order to keep all solution components in
     *  all domains not-appreciably changing so much that the jacobian isn't representative.
     *
     *  Return the factor by which the undamped Newton step 'step0'
     *  must be multiplied in order to keep all solution components in
     *  all domains between their specified lower and upper bounds.
     * 
     *  This routine is meant to be used with cropping. In other words,
     *  each component is allowed to go slightly out of bounds. However, cropping may be used to enforce a strict limit.
     *
     *  Currently the bounds are hard coded into this routine:
     *
     *  Minimum value for all variables: solnLowBound[i] - 0.01 * m_ewt[i]
     *  Maximum value = none. solnHighBound[i] + 0.01 * m_ewt[i]
     *
     *  On any single step, each solution component is not allowed to go more than 9/10 of the way to the boundary, from where 
     *  it currently is.
     *
     *  If the component is already out of bounds, bounds checking is no longer carried out on that component.
     *
     *  @param[in]           y                   Ghosted Epetra_Vector reference for the current solution 
     *  @param[in]           step0               Owned Epetra_Vector reference for the current step in the solution unknowns
     *  @param[in]           loglevel            Level of printing. If greater than 3, a line is printed out about the
     *                                           damping caused by this routine.
     *
     *  @return                                  Returns the factor which the step should be reduced by.
     */
    virtual double highLowBoundStep(const Epetra_Vector_Ghosted& y, const Epetra_Vector_Owned& step0, int loglevel);

    //! Set options.
    /*!
     *  (virtual from SolGlobalNonlinear)
     *  @param[in]           maxJacAge           Maximum age of the jacobian
     */
    virtual void setOptions(int maxJacAge = 5)
    {
    }

    //! Set the absolute tolerances for the solution variables
    /*!
     *  (virtual from SolGlobalNonlinear)
     *  Set the absolute tolerances used in the calculation
     *
     *  @param[in]           reltol              relative tolerance used in the nonlinear solver
     *  @param[in]           n                   Length of abstol. Should be equal to m_NumLcEqns
     *  @param[in]           abstol              Vector of length n that contains the tolerances to be used
     *                                           for the solution variables
     */
    virtual void setTolerances(double reltol, int n, const double* const abstol);

    //! Set the absolute tolerances for the solution variables for delta damping
    /*!
     *  (virtual from SolGlobalNonlinear)
     *  Set the absolute tolerances used in the delta damping algorithm. Essentially we shouldn't control 
     *  step sizes if the step is below the absolute tolerance
     *
     *  @param[in]           reltol_dd           relative tolerance used in the delta damping algorithm.
     *  @param[in]           n                   Length of abstol. Should be equal to m_NumLcEqns
     *  @param[in]           abstol_dd           Vector of length n that contains the tolerances to be used 
     *                                           for the solution variables
     */
    virtual void setTolerances_deltaDamping(double reltol_dd, int n, const double* const abstol_dd);

    //! Set the value of the maximum # of newton iterations
    /*!
     *  (virtual from SolGlobalNonlinear)
     *  @param[in]           maxNewtIts          Maximum number of newton iterations.
     *                                           The default value of this is 50 iterations
     */
    virtual void setMaxNewtIts(const int maxNewtIts);

    //! Set the problem type
    /*!
     *  (virtual from SolGlobalNonlinear)
     *  @param[in]           probtype            problem type
     */
    virtual void setProblemType(int probtype);

    //! Set the solution weights based on the atol and rtol values
    /*!
     *  (virtual from SolGlobalNonlinear)
     */
    virtual void setDefaultSolnWeights();

    //! Set the boolean for turning on and off row scaling
    /*!
     *  (virtual from SolGlobalNonlinear)
     *  Row scaling is on by default. It's nearly always a good choice to have it on.
     *
     *  @param[in]           onoff               True if you want row scaling
     */
    virtual void setRowScaling(const bool onoff);

    //! Toggle that turns on and off column scaling
    /*!
     *  (virtual from SolGlobalNonlinear)
     *
     *  Columns scaling is turned on by default
     *
     *  @param[in]           onoff               toggle
     *  @param[in]           colScaleUpdateFrequency column scale update frequency
     *                                                  0: never
     *                                                  1: once at the start
     *                                                  2: after every jac update
     */
    virtual void setColScaling(const bool onoff, const int colScaleUpdateFrequency);

    //! Set the toggles for solution damping
    /*!
     *  (virtual from SolGlobalNonlinear)
     *
     *  Residual solution damping means that line-step updates will not be accepted if the residual norm is larger
     *  than the previous residual norm. We shrink the step size until we find an update with a smaller
     *  residual norm, or if the weigted residual norm is less than one, then we accept the step as well.
     *
     *  Delta damping means that a solution component will not be allowed to be changed beyond a certain factor
     *  or range in any one step. The default is on.
     *
     *  High/Low damping means that all solution components have high and low bounds that they can't cross.
     *  
     *
     *  @param[in]           residSolnDamping    Toggle damping due to the value of the residual
     *                                              (defaults to on) 
     *  @param[in]           deltaDamping        Toggle damping due to the delta damping criteria
     *                                              (defaults to on) 
     *  @param[in]           highLowDamping      Toggle damping due to max and min bounds on each variable.
     *                                              (defaults to on) 
     */
    virtual void setDampingToggles(const bool residSolnDamping, const bool deltaDamping, const bool highLowDamping);

    //! Set the vectors for lower and upper boundaries.
    /*!
     *  (virtual from SolGlobalNonlinear)
     *
     *  @param[in]           lowBounds           Epetra owned vector reference containing the Low bounds 
     *                                           for all solution componets
     *  @param[in]           highBounds          Epetra owned vector reference containing the high bounds 
     *                                           for all solution components
     */
    virtual void setSolutionBounds(const Epetra_Vector_Owned& lowBounds, const Epetra_Vector_Owned& highBounds);

    //! Sets the print solution options
    /*!
     *  (virtual from SolGlobalNonlinear)
     *
     *  @param[in]           dumpJacobians       Boolean indicating whether jacobians are to be dumped out
     */
    virtual void setPrintSolnOptions(bool dumpJacobians);

    //! Set some nonlinear solver options
    /*!
     *  (virtual from SolGlobalNonlinear)
     *
     *  @param[in]           min_newt_its        Set the minimum number of newton iterations to carry out at every step.
     *                                           A value of 2 is useful to monitor convergence.
     *                                           (defaults to 0)
     *  @param[in]           matrixConditioning  If true, carry out matrix conditioning process before trying to find the
     *                                           inverse of the matrix (Not implemented).
     *                                           (defaults to false)
     *  @param[in]           colScaling          Implement column scaling to the matrix system. The unknowns are then scaled
     *                                           inversely to the matrix. Scaling is carried out so that the deltas at
     *                                           convergence of the matrix are scaled to a value of one. This means that
     *                                           the column scales are equal to the weighting matrix.
     *                                           Column scaling is carried out before row scaling.
     *                                           (default = false)
     *  @param[in]           rowScaling          Implement row scaling to the matrix system. The max row element in the
     *                                           jacobian is caled to one by multiplying all terms in the row and rhs
     *                                           by a scale factor. Note, this is effective in practise in reducing
     *                                           the condition number of the matrix.
     *                                           (default = true).
     *  @param[in]           colScalingUpdateFrequency frequency of column scaling
     *                                             0  Column scales are never updated by this routine. They are
     *                                                an input to this routine.
     *                                             1  Column scales are updated once at the start of the sol_nonlinear_problem
     *                                                call. A call to setDefaultColScales() is made just after the
     *                                                solution weights are evaluated (default).
     *                                             2  Column scales are updated after each jacobian evaluation.
     *                                                A call to setDefaultColScales() is made during the scaleMatrix routine.
     *                                               (Defaults to 1)
     */
    virtual void 
    setNonLinOptions(int min_newt_its = 0, bool matrixConditioning = false, bool colScaling = false, 
                     bool rowScaling = true, int colScalingUpdateFrequency = 1);

    //!  Set the level of printing that occurs during the nonlinear solve
    /*!
     *  Printing is done to stdout.
     *
     *  @param[in]           print_flag          Set up the level of printing according to the table below
     *
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

    //! Supply a predicted solution which will be used to help set up the solution weights
    /*!
     *  (virtual from SolGlobalNonlinear)
     *
     *  This routine is used to set up the scales for the solution error weights. 
     * 
     *  @param[in]           y_pred              Value of the predicted solution with ghost values
     */
    virtual void 
    setPredicted_soln(const Epetra_Vector_Ghosted& y_pred);

    //! L2 Weighted Norm of a delta in the solution
    /*!
     *  The vector m_ewt[i]'s are always used to weight the solution errors in the calculation.
     *
     *  The second argument has a default of false. However,
     *  if true, then a table of the largest values is printed  out to standard output.
     *
     *  @param[in]           delta_soln          Norm of a delta of the solution vector
     *  @param[in]           printLargest        if True a table is printed of the largest contributors.
     *  @param[in]           title               Printed out on the title line
     *  @param[in]           typeYsoln           Parameter indicating whether m_y_curr[] is currently
     *                                           evaluated before the delta or after the delta has been implemented.
     *  @param[in]           dampFactor          Damping factor that will be applied to delta_y before creating a new ysoln
     *
     *  @return                                  L2 weighted norm of the solution update vector
     */
    virtual double
    soln_error_norm(const Epetra_Vector& delta_soln, const bool printLargest = false, const char* title = 0,
                    const int typeYsoln = 1, const double dampFactor = 1.0) const;

    //!  L2 Weighted Norm of the residual
    /*!
     *   The vector m_residWt[i]'s are always used to weight the solution errors in the calculation.
     *
     *   The second argument has a default of false. However, if true, then a table of the largest values is printed
     *   out to standard output.
     *
     *   @param[in]          resid               Owned Epetra vector reference containing the computed residual
     *   @param[in]          title               Printed out on the title line
     *   @param[in]          printLargest        if True a table is printed of the largest contributors.
     *
     *   @return                                 L2 weighted norm of the residual vector
     */
    virtual double
    res_error_norm(const Epetra_Vector_Owned& resid, const char* title = 0, const int printLargest = 0) const;

   //! Evaluate the residual at the current conditions
    /*!
     *  The evaluated residual is stored in the internal vector, m_resid.
     *
     *  @param[in]           time_curr           Current value of the time
     *  @param[in]           rdelta_t            one over the current deltaT value
     *  @param[in]           solnBase_ptr        ptr to Ghosted Epetra vector containing the current value of the solution
     *  @param[in]           solnDotBase_ptr     ptro to Ghosted Epetra vector containing the current value of the solution dot
     */
    virtual void
    get_res(const double time_curr, const double rdelta_t,
            const Epetra_Vector_Ghosted* const solnBase_ptr, const Epetra_Vector_Ghosted* const solnDotBase_ptr);

public:

    //! Main routine to launch a nonlinear solve at the current conditions
    //! whether it's a steady state solve or a time-dependent run.
    /*!
     *   Find the solution to F(X) = 0 by damped Newton iteration.  On
     *   entry, y_comm contains an initial estimate of the solution.  On successful return, y_comm contains the converged solution.
     *
     *  @param[in]           solveType           = TRANSIENT -> we will assume we are relaxing a transient equation system for now.
     *                                                          Will make it more general later,if an application comes up.
     *                                           = DAEINIT Solve the initial conditions problem for a DAE system.

     *  @param[in,out]       y_comm              ptr to Ghosted Epetra vector containing the current value of the solution
     *  @param[in,out]       ydot_comm           ptr to Ghosted Epetra vector containing the current value of the solution dot
     *  @param[in]           CJ                  Proportionality factor between the dependence of ydot on the solution y.
     *                                           (Comes from the BDF formulas).
     *  @param[in]           time_curr           Current value of the time
     *  @param[out]          num_newt_its        Number of newton iterations it took to solve the system
     *  @param[out]          num_linear_solves   Number of linear solution steps it took during the linear solves
     *  @param[out]          num_backtracks      Number of back tracks along line searches taken to solve the system.
     *
     *  @return                                  Return code:
     *                                             1 :  Successful step was taken: Next step was less than previous step.
     *                                                  s1 is calculated
     *                                             2 :  Successful step: Next step's norm is less than 0.8
     *                                             3 :  Success:  The final residual is less than 1.0
     *                                                  A predicted deltaSoln is not produced however. s1 is estimated.
     *                                             4 :  Success:  The final residual is less than the residual
     *                                                  from the previous step.
     *                                                  A predicted deltaSoln is not produced however. s1 is estimated.
     *                                             1 :  Nonlinear problem solved successfully.
     *
     *                                            -1 :  We ran out of the allowed newton iterations. Nothing happened
     *                                                  that was a show stopper. We just ran out of iterations
     *                                                  Return the current solution as the best case.
     *                                            -2 :
     *                                            -3 :  We ran into a hard bounds problem. This is a show stopper, because
     *                                                  it is a hard bounds. This usually indicates that the either the nonlinear
     *                                                  solution path veers into prohibited territory before settling down into
     *                                                  valid territory, or, it may mean that there is sufficient round off error
     *                                                  from ill-conditioning that the step goes into prohibited territory.
     */
    virtual int
    solve_nonlinear_problem(Zuzax::Solve_Type solveType, Epetra_Vector_Ghosted* const y_comm,
                            Epetra_Vector_Ghosted* const ydot_comm, double CJ,
                            double time_curr, int& num_newt_its, int& num_linear_solves, int& num_backtracks);

    //!  Scale the matrix.
    /*!
     * @param delta_soln
     * @param y_curr
     * @param ydot_curr
     * @param time_curr
     * @param rdelta_t
     * @param loglevel
     */
    virtual void
    scaleMatrix(Epetra_Vector_Owned& delta_soln,
                const Epetra_Vector_Ghosted& y_curr,
                const Epetra_Vector_Ghosted& ydot_curr,
                const double time_curr,
                const double rdelta_t,
                int loglevel);

    //! Compute the undamped Newton step.
    /*!
     * @param delta_soln
     * @param y_curr
     * @param ydot_curr
     * @param time_curr
     * @param rdelta_t
     * @param loglevel
     */
    virtual void
    doNewtonSolve(Epetra_Vector_Owned& delta_soln,
                  const Epetra_Vector& y_curr,
                  const Epetra_Vector& ydot_curr,
                  const double time_curr,
                  const double rdelta_t,
                  int loglevel);



    //! Set the column scales used in the program
    /*!
     *
     * @param colScales
     */
    virtual void
    setColumnScaleVector(const Epetra_Vector_Owned& colScales);

    //! Set the column scaling vector at the current time
    /*!
     * The column scales are set equal to the current value of the weighting vector
     */
    virtual void
    setDefaultColumnScaleVector();

    //! return the column scales used in the program
    /*!
     *
     * @param colScales
     * @return returns true if column scaling is being used
     */
    virtual bool
    getColumnScaleVector(Epetra_Vector_Owned& colScales) const;

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

    //! Compute the Residual Weights
    /*!
     *  The residual weights are defined here to be equal to the inverse of the row scaling factors used to
     *  row scale the matrix, after column scaling is used. They are multiplied by 10-3 because the column
     *  weights are also multiplied by that same quantity.
     *
     *  The basic idea is that a change in the solution vector on the order of the convergence tolerance
     *  multiplied by  [RJC] which is of order one after row scaling should give you the relative weight
     *  of the row. Values of the residual for that row can then be normalized by the value of this weight.
     *  When the tolerance in delta x is achieved, the tolerance in the residual is also achieved.
     */
    virtual void
    computeResidWts();

    //! Get the Residual Weights
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
     *  this routine reports the weights.
     *
     * @param residWts  Returns the residual weights in a vector of length
     */
    virtual void
    getResidWts(Epetra_Vector_Owned& residWts);

    //! Calculate a consistent ydot
    /*!
     * Function to calculate the acceleration vector ydot for the first or
     * second order predictor/corrector time integrator.  This routine can be
     * called by a first order - forward Euler / backward Euler predictor /
     * corrector or for a second order Adams - Bashforth / Trapezoidal Rule
     * predictor / corrector.  See Nachos documentation Sand86-1816 and Gresho,
     * Lee, Sani LLNL report UCRL - 83282 for more information.
     *
     *  other variables:
     *
     *      delta_t_n   - Magnitude of the current time step at time n  (i.e., = t_n - t_n-1)
     *      y_nm1[]     - Solution vector at time n-1
     *      ydot_nm1[]  - Acceleration vector at time n-1
     *
     *  output:
     *
     *      ydot_curr[]   - Current acceleration vector at time n
     *
     * Note we use the current attribute to denote the possibility that
     * y_curr[] may not be equal to m_y_n[] during the nonlinear solve because we may be using a look-ahead scheme.
     *
     *  @param[in]           order               Order of the time stepping method
     *                                              = 1 -> first order forward Euler/backward Euler
     *                                                     predictor/corrector
     *                                              = 2 -> second order Adams-Bashforth/Trapezoidal Rule
     *                                                     predictor/corrector
     *  @param[in]           y_curr              Current value of the solution vector at time n
     *  @param[out]          ydot_curr           Calculated value of the solution dot vector at time n
     */
    virtual void
    calc_ydot(int order, const Epetra_Vector& y_curr, Epetra_Vector& ydot_curr);


    /************************************************ Member data ***********************************/
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
    int m_print_flag;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------

#endif
