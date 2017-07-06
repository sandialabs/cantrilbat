/**
 *  @file m1d_SolNonlinear.h
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef M1D_SOLNONLINEAR_H
#define M1D_SOLNONLINEAR_H

#include "m1d_SolGlobalNonlinear.h"
#include "m1d_exception.h"

#include "m1d_EpetraJac.h"
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//!  Exception class thrown when an error in the nonlinear solver is encountered.
/*!
 *
 */
class SolNonlinearErr : public m1d_Error
{
public:

    //! Constructor 
    /*!
     *  @param[in]           msg                 Message to be written
     */
    SolNonlinearErr(std::string msg);
};
//==================================================================================================================================

//! Analytical jacobian
#define BEULER_JAC_ANAL 2
//! Numerical jacobian
#define BEULER_JAC_NUM  1

//==================================================================================================================================
//!  This is the nonlinear solver.
/*!
 *   This is the main nonlinear solver for steady state and time-dependent problems.
 *   It employs a two-way convergence criteria: one for the residual and one for the solution update vector
 *
 */
class SolNonlinear : public SolGlobalNonlinear
{
public:

    //!  The default constructor doesn't take an argument.
    SolNonlinear();

    //! Destructor
    virtual ~SolNonlinear();

    //! Setup the problem for solution.
    /*!
     *  (virtual from SolGlobalNonlinear)
     *
     *   Find the solution to F(X) = 0 by damped Newton iteration.  On entry, x0 contains an initial estimate of the solution. 
     *   on successful return, x1 contains the converged solution.
     *
     *  @param[in]           solveType           Enum of type Solve_Type_Enum. Describes the type of the problem to be solved.
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
    setup_problem(Solve_Type_Enum solveType, const Epetra_Vector_Ghosted* const y_init,
                  const Epetra_Vector_Ghosted* const ydot_init, double time_curr, ProblemResidEval& problem, 
                  EpetraJac& jac) override;

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
    virtual int
    doHardBounds(const Epetra_Vector_Ghosted& y_old, Epetra_Vector_Owned& step, double& fbound) override;

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
    virtual double
    deltaBoundStep(const Epetra_Vector_Ghosted& y, const Epetra_Vector_Owned& step0) override;

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
     * Currently the bounds are hard coded into this routine:
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
    virtual double
    highLowBoundStep(const Epetra_Vector_Ghosted& y, const Epetra_Vector_Owned& step0, int loglevel) override;

    //! Set options.
    /*!
     *  (virtual from SolGlobalNonlinear)
     *  @param[in]           maxJacAge           Maximum age of the jacobian
     */
    virtual void setOptions(int maxJacAge = 5) override
    {
        m_maxAge = maxJacAge;
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
    virtual void setTolerances(double reltol, int n, const double* const abstol) override;

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
    virtual void setTolerances_deltaDamping(double reltol_dd, int n, const double* const abstol_dd) override;

    //! Set the value of the maximum # of newton iterations
    /*!
     *  (virtual from SolGlobalNonlinear)
     *  @param[in]           maxNewtIts          Maximum number of newton iterations.
     *                                           The default value of this is 50 iterations
     */
    virtual void setMaxNewtIts(const int maxNewtIts) override;

    //! Set the problem type
    /*!
     *  (virtual from SolGlobalNonlinear)
     *  @param[in]           probtype            problem type
     */
    virtual void setProblemType(int probtype) override;

    //! Set the solution weights based on the atol and rtol values
    /*!
     *  (virtual from SolGlobalNonlinear)
     */
    virtual void setDefaultSolnWeights() override;

    //! Set the boolean for turning on and off row scaling
    /*!
     *  (virtual from SolGlobalNonlinear)
     *  Row scaling is on by default. It's nearly always a good choice to have it on.
     *
     *  @param[in]           onoff               True if you want row scaling
     */
    virtual void setRowScaling(const bool onoff) override;

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
    virtual void setColScaling(const bool onoff, const int colScaleUpdateFrequency) override;

    //! Set the toggles for solution damping
    /*!
     *  (virtual from SolGlobalNonlinear)
     *
     *  @param[in]           residSolnDamping    Toggle damping due to the value of the residual
     *  @param[in]           deltaDamping        Toggle damping due to the delta damping criteria
     *  @param[in]           highLowDamping      Toggle damping due to max and min bounds on each variable.
     */
    virtual void 
    setDampingToggles(const bool residSolnDamping, const bool deltaDamping, const bool highLowDamping) override;

    //! Set the vectors for lower and upper boundaries.
    /*!
     *  (virtual from SolGlobalNonlinear)
     *
     *  @param[in]           lowBounds           Epetra owned vector reference containing the Low bounds 
     *                                           for all solution componets
     *  @param[in]           highBounds          Epetra owned vector reference containing the high bounds 
     *                                           for all solution components
     */
    virtual void
    setSolutionBounds(const Epetra_Vector_Owned& lowBounds, const Epetra_Vector_Owned& highBounds) override;

    //! Print the largest contributors to the solution update on a relative term
    /*!
     *  A table of the largest contributors to the solution update is printed out, according to the weights.
     *  This is useful in debugging to understand what is not converging.
     *
     *  @param[in]           solnDelta0          Delta in the current solution
     *  @param[in]           s0                  Title for the solnDelta0 contribution usually "DeltaSoln"
     *  @param[in]           solnDelta1          Next soln delta step after the current delta is taken
     *  @param[in]           s1                  title for the solnDelta1 contribution usually "DeltaSolnTrialTest"
     *  @param[in]           title               Title of the table , usually
     *                                                  "dampNewt: Important Entries for Weighted Soln Updates"
     *  @param[in]           y0                  Old soln
     *  @param[in]           y1                  New solution after solnDelta0 is applied
     *  @param[in]           damp                current value of the damping factor
     *  @param[in]           num_entries         Number of entries in the table to be printed out
     */
    void
    print_solnDelta_norm_contrib(const Epetra_Vector& solnDelta0, const char* const s0,
                                 const Epetra_Vector& soln1, const char* const s1,
                                 const char* const title, const Epetra_Vector& y0, const Epetra_Vector& y1,
                                 double damp, int num_entries);

    //!  Update the solution vector using the step change that was just computed.
    /*!
     *    We update the solution vector and the solution dot vector (this is the time derivative vector),
     *    putting the answer into the fixed location,
     *       m_y_new[]   and m_ydot_new[]
     *    given the previous solution vector,
     *          y0[]     and ydot0_ptr[]
     *    and the update step vector with a damping factor
     *          ff           step_1[]
     *
     *   @param  y0         INPUT     Input solution vector     - Ghosted Epectra_Vector reference
     *   @param  ydot0_ptr  INPUT     Input solution dot vector - Ghosted Epectra_Vector ptr
     *   @param  ff         INPUT     Damping factor - double
     *   @param  step_1     INPUT     Input step vector         - Ghosted Epectra_Vector reference
     *
     *    OUTPUT
     * ----------------
     *    This routine changes
     *      m_y_new         OUTPUT     New solution vector     - Ghosted Epectra_Vector ptr
     *      m_ydot_new      OUTPUT     New solution dot vector - Ghosted Epectra_Vector ptr
     *
     *   DISCUSSION
     * ----------------
     *
     *    This routine will update the locally owned values. Then, it will call updateGhostEqns()
     *    for both the m_y_new and  m_ydot_new vectors to update the ghost unknowns on neighboring processors.
     *
     *    For the  DAESystemInitial_Solve problem this routine will scatter the step vector unknowns
     *    into the correct locations in the m_y_new and *m_ydot_new vectors according to whether the
     *    the value of (*m_isAlgebraic)[j] is equal to one or not. For DAE unknowns, theoretically
     *    it doesn't matter what the value of the time derivative is, since it doesn't enter into
     *    the equation set. However, here we set (*m_ydot_new)[j] = 0.0 for DAE unknowns to avoid
     *    doing nothing with the entry.
     *
     *  @param y0    Base solution vector for the current time. The delta is applied to this variable.
     *  @param ydot0_ptr   Base Solution dot vector for the current time. The delta is applied
     *                     to this after division of the time step coefficient.
     *  @param ff    Damping coefficient to use. The step is reduced by this factor.
     *  @param step0 Step to be taken. If we are doing the DAE initial problem, this is a
     *               mixture of the base solution step and the ydot solution step, depending
     *               upon the value of (*m_isAlgebraic)[j].
     */
    void updateSoln(const Epetra_Vector_Ghosted& y0, const Epetra_Vector_Ghosted* ydot0_ptr,
                    double ff, const Epetra_Vector_Ghosted&  step0);

    //! Sets some printing options
    /*!
     *  (virtual from SolGlobalNonlinear)
     *
     *  @param[in]           dumpJacobians       Dump jacobians to disk.
     */
    virtual void
    setPrintSolnOptions(bool dumpJacobians) override;

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
    setNonLinOptions(int min_newt_its = 0, bool matrixConditioning = false, bool colScaling = false, bool rowScaling =
                         true, int colScalingUpdateFrequency = 1) override;

    //!  Supply a predicted solution which will be used to help set up the solution weights
    /*!
     *  (virtual from SolGlobalNonlinear)
     *
     *  This routine is used to set up the scales for the solution error weights. 
     * 
     *  @param[in]           y_pred              Value of the predicted solution with ghost values
     */
    virtual void
    setPredicted_soln(const Epetra_Vector_Ghosted& y_pred) override;

    //!    L2 Weighted Norm of a delta in the solution
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
     *   @param title      Printed out on the title line
     *   @param typeYsoln  Parameter indicating whether m_y_curr[] is currently
     *                     evaluated before the delta or after the delta has been implemented.
     *   @param dampFactor Damping factor that will be applied to delta_y before creating a new ysoln
     */
    virtual double
    soln_error_norm(const Epetra_Vector_Owned& delta_y,
                    const bool printLargest = false,
                    const char* title = 0,
                    const int typeYsoln = 1,
                    const double dampFactor = 1.0) const;

    //!    L2 Weighted Norm of the residual
    /*!
     *   The vector m_residWt[i]'s are always used to weight the solution errors in
     *   the calculation.
     *
     *   The second argument has a default of false. However,
     *   if true, then a table of the largest values is printed
     *   out to standard output.
     *
     *   @param delta_y  Norm of a delta of the solution vector
     *   @param title      Printed out on the title line
     *   @param printLargest if True a table is printed of the largest contributors.
     */
    virtual double
    res_error_norm(const Epetra_Vector_Owned& resid, const char* title = 0, const int printLargest = 0) const;

    //!  Function called by SolNonlinear to evaluate the Jacobian matrix and the
    /*!
     *
     * @param jac
     * @param doTimeDependentResid
     * @param time_curr
     * @param rdelta_t
     * @param solnBase_ptr
     * @param solnDotBase_ptr
     * @return
     */
    int
    get_jac(EpetraJac& jac,
            Epetra_Vector_Owned* res,
            const bool doTimeDependentResid,
            double time_curr,
            double rdelta_t,
            const Epetra_Vector_Ghosted* solnBase_ptr,
            const Epetra_Vector_Ghosted* solnDotBase_ptr);

    //! get the residual
    /*!
     *
     * @param time_curr
     * @param rdelta_t
     * @param solnBase_ptr
     * @param solnDotBase_ptr
     */
    virtual void
    get_res(const double time_curr,
            const double rdelta_t,
            const Epetra_Vector_Ghosted* solnBase_ptr,
            const Epetra_Vector_Ghosted* solnDotBase_ptr);

public:

    //! Main routine to launch a nonlinear solve at the current conditions
    //! whether it's a steady state solve or a time-dependent run.
    virtual int
    solve_nonlinear_problem(Solve_Type_Enum solveType,
                            Epetra_Vector_Ghosted* y_comm,
                            Epetra_Vector_Ghosted* ydot_comm,
                            double CJ,
                            double time_curr,
                            int& num_newt_its,
                            int& num_linear_solves,
                            int& num_backtracks);

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
    doNewtonSolve(Epetra_Vector_Owned& delta_soln, const Epetra_Vector& y_curr, const Epetra_Vector& ydot_curr,
                  const double time_curr, const double rdelta_t, int loglevel);

    //! Attempt to find a damping step using a line search algorithm that leads to a better solution
    /*!
     *  On entry, the member variable m_stp  contains an undamped Newton step for the solution (y0, ydot0_ptr). 
     *  This method attempts to  find a damping coefficient such that the next undamped step would have
     *  a norm smaller than that of step0. If successful, the new solution after taking the
     *  damped step is returned in y1, and the undamped step at y1 is returned in step1.
     *
     * On entry, the member variable m_stp must contain the Newton step for the
     * solution y0. The step size at this point may have already gone through a HighLow bounds check and a 
     * delta solution bounds check. Therefore, it may have already been reduced. The factor for the reduction is kept
     * in the variable m_fbound.
     * 
     *  This method attempts to find a damping coefficient for the current value of m_stp
     *  such that the next undamped step would have a norm smaller than
     *  that of step0. If successful, the new solution after taking the
     *  damped step is returned in m_y_new and m_ydot_new, and the undamped next step from m_y_new[] is storred 
     *  in m_step_2[].
     *
     *  However, before it tries a new newton step, it will accept a step if the new residual norm is less than the first
     *  residual norm or if the new residual norm is less than one. This option leads to return codes of 3 and 4. In this
     *  case the weighted solution norm for the next step, s1, is estimated and not calculated.
     *
     *  @param[in]           time_curr           Current time
     *  @param[in]           y0                  Current value of the solution vector
     *  @param[in]           ydot0_ptr           Current pointer to the solution Time derivative
     *  @param[in]           s1                  Overall damping factor for the step that produces a viable result.
     *  @param[in]           loglevel            Loglevel controls the amount of printing
     *  @param[out]          num_backtracks      Returns the number of backtrack steps
     *
     *  @return                                  1 Successful step was taken: Next step was less than previous step.
     *                                                                        s1 is calculated
     *                                           2 Successful step: Next step's norm is less than 0.8
     *                                           3 Success:  The final residual is less than 1.0.
     *                                                       A predicted deltaSoln is not produced however. s1 is estimated.
     *                                           4 Success:  The final residual is less than the residual from the previous step.
     *                                                       A predicted deltaSoln is not produced however. s1 is estimated.
     *                                           0 Uncertain Success: s1 is about the same as s0
     *                                          -2 Unsuccessful step.
     *
     *  Return member data:
     *
     *       m_y_new[]         New solution vector to be used in the next nonlin step 
     *       m_ydot_new[]      New solutionDot vector to be used in the next nonlin step
     *       m_step_2[]        New raw step vector computed from m_y_new[] using a Newton's method.
     */
    int dampStep(double time_curr, const Epetra_Vector& y0, const Epetra_Vector* ydot0_ptr, double& s1,
                 int& loglevel, int& num_backtracks);


    //! Attempt to find a damping step that leads to a better solution
    /*!
     *  On entry, the member variable m_stp  contains an undamped Newton step for the solution (y0, ydot0_ptr). 
     *  This method attempts to  find a damping coefficient such that the next undamped step would have
     *  a norm smaller than that of step0. If successful, the new solution after taking the
     *  damped step is returned in y1, and the undamped step at y1 is returned in step1.
     *
     *  @param[in]           time_curr           Current time
     *  @param[in]           y0                  Current value of the solution vector
     *  @param[in]           ydot0_ptr           Current pointer to the solution Time derivative
     *  @param[in]           s1                  Overall damping factor for the step that produces a viable result.
     *  @param[in]           loglevel            Loglevel controls the amount of printing
     *  @param[out]          num_backtracks      Returns the number of backtrack steps
     *
     *  @return                                  1 Successful step was taken: Next step was less than previous step.
     *                                                                        s1 is calculated
     *                                           2 Successful step: Next step's norm is less than 0.8
     *                                           3 Success:  The final residual is less than 1.0
     *                                                       A predicted deltaSoln is not produced however. s1 is estimated.
     *                                           4 Success:  The final residual is less than the residual
     *                                                       from the previous step.
     *                                                       A predicted deltaSoln is not produced however. s1 is estimated.
     *                                           0 Uncertain Success: s1 is about the same as s0
     *                                          -2 Unsuccessful step.
     */
    int dampStep_alt(double time_curr, const Epetra_Vector& y0, const Epetra_Vector* ydot0_ptr, double& s1, int& loglevel,
                     int& num_backtracks);


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
     *  We set the values for the previous time step here. These are used in the nonlinear
     *  solve, because they affect the calculation of ydot.
     *
     *  @param[in]           timeStep            Time step between current time and previous time
     *  @param[in]           y_nm1               Value of the solution vector at the previous time step
     *  @param[in]           ydot_nm1            Value of the solution vector derivative at the previous time step
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
     *  row scale the matrix, after column scaling is used. This routine returns the residual weights.
     *
     *  @param[in]           residWts            Returns the residual weights in a vector of length
     */
    virtual void getResidWts(Epetra_Vector_Owned& residWts);

    //! Print a solution-like vector in a form that is easy to interpret
    /*!
     *  Prints a vector that has the same layout as the solution vector.
     *  The variable names and node positions are printed along with the values.
     *
     *  This routine uses the member function showSolutionVector() of the
     *  ProblemResidEval class, which does the layout printing. It can be used to print to a file as well.
     *
     *  @param[in]           header              Header that is printed out to
     *  @param[in]           v                   Epetra vector -> Note only the owned elements are needed.
     *                                           (mp behavior hasn't been checked)
     */
    void print_SVector(std::string header, const Epetra_MultiVector& v) const;


    //! Print a solution-like int vector in a form that is easy to interpret
    /*!
     *  Prints an int vector that has the same layout as the solution vector.
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
    void print_IntSVector(std::string header, const Epetra_IntVector& v) const;


    //! Calculate a consistent ydot
    /*!
     * Function to calculate the acceleration vector ydot for the first or
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
    virtual void
    calc_ydot(int order, const Epetra_Vector& y_curr, Epetra_Vector& ydot_curr);

private:
    //! Check to see if the nonlinear problem has converged
    /*!
     *  @param[in]           dampCode            Code returned from the line solve damping algorithm
     *
     *  @param[in]           s1                  Predicted value of the normalized step size in the next nonlinear iteration
     *
     * @return integer is returned. If positive, then the problem has converged
     *           1 Successful step was taken: Next step's norm is less than 1.0.
     *                                        The final residual norm is less than 1.0.
     *           2 Successful step: Next step's norm is less than 0.8.
     *                              This step's norm is less than 1.0.
     *                              The residual norm can be anything.
     *           3 Success:  The final residual is less than 1.0
     *                        The predicted deltaSoln is below 1.0.
     *           0 Uncertain Success: s1 is about the same as s0
     */
    int
    convergenceCheck(int m, double s1);

protected:
    /************************************************ Member data ***********************************/
    /*********************
     * METHOD FLAGS
     *********************/

    //! Nonlinear problem type
    /*!
     *  0 SteadyState
     *  1 TimeDependentAccurate
     *  2 TimeDependentRelax
     *  3 DAE
     */
    Solve_Type_Enum solnType_;

    //! m_jacFormMethod determines how a matrix is formed.
    int m_jacFormMethod;

    //!  If true then row sum scaling of the Jacobian matrix is carried out when solving the linear systems.
    bool m_rowScaling;
    
    //!  If true, then column scaling is performed on each solution of the linear system.
    bool m_colScaling;

    //! Frequency with which the column scales are updated
    /*!
     *   0  Column scales are never updated by this routine. They are
     *      an input to this routine.
     *   1  Column scales are updated once at the start of the sol_nonlinear_problem
     *      call. A call to setDefaultColScales() is made just after the
     *      solution weights are evaluated (default).
     *   2  Column scales are updated after each jacobian evaluation.
     *      A call to setDefaultColScales() is made during the scaleMatrix
     *      routine.
     */
    int colScaleUpdateFrequency_;

    //!  Boolean that turn son  matrix conditioning
    /*!
     *  If true, then the Jacobian and every rhs is multiplied by the inverse
     *  of a matrix that is suppose to rfeduce the condition number of the matrix. This is done before row scaling.
     */
    bool m_matrixConditioning;

    //! Relative nonlinear solver truncation error tolerances
    double m_reltol;
    
    //! Vector of absolute nonlinear solver truncation error tolerance when not uniform for all variables.
    /*!
     *  only the owned unknowns
     */
    Epetra_Vector_Owned* m_abstol;

    //! Solution error weighting vector
    /*!
     *  pointer to an Epetra_vector for owned unknowns
     */
    Epetra_Vector_Owned* m_ewt;

    //! Solution delta damping weighting vector
    /*!
     *  pointer to an Epetra_vector for owned unknowns
     */
    Epetra_Vector_Owned* m_ewt_deltaDamping;

    //! Absolute tolerance criteria for delta damping
    /*!
     *  pointer to an Epetra_vector for owned unknowns
     */
    Epetra_Vector_Owned* m_absTol_deltaDamping;

    //! Pointer to the Epetra_Comm object
    const Epetra_Comm* Comm_ptr_;

    //! Current integration order
    int m_order;

    //! Failure Counter -> keeps track of the number * of consequetive failures
    int m_failure_counter;

    /**
     * Minimum Number of Newton Iterations per nonlinear step
     * default = 0
     */
    int m_min_newt_its;

    //! Maximum number of newton iterations
    int maxNewtIts_;

    //! Age of the current jacobian
    int m_jacAge;

    //! Number of linear iterations used to solve the linear system
    int m_curr_linearIts;

    //! Norm of the solution for the linear system
    double m_curr_normLin;

    /****************************************************************************
     *    METHOD OPTIONS FOR THE NONLINEAR SOLVER
     ****************************************************************************/

    //! Toggle solution damping
    bool doResidSolnDamping_;

    //! Toggle delta damping
    bool doDeltaDamping_;

    //! Toggle High Low Bounds Damping
    bool doHighLowDamping_;

    //! Number of equations in the ode integrator on the current processor
    int m_NumLcEqns;
    
    //! Number of equations in the ode integrator on the current processor that are owned by the current processor.
    int m_NumLcOwnedEqns;

    //! Number of equations on all processors
    int m_NumGbEqns;

    /******************************************************************************
     *    METHOD OPTIONS FOR THE LINEAR SOLVER
     ********************************************************************************/

    //! Value of the current solutionDot vector with ghost node values
    Epetra_Vector_Ghosted* m_y_curr;

    //! Value of the current solutionDot vector with only owned values
    Epetra_Vector_Owned* m_y_curr_owned;

    //! Value of the proposed new solution vector with ghost node values
    Epetra_Vector_Ghosted* m_y_new;

    //! Value of the proposed new solution vector with only owned values
    Epetra_Vector_Owned* m_y_new_owned;

    //! Stored solution vector at the previous time step
    Epetra_Vector_Ghosted* m_y_nm1;

    //! "Predicted" solution vector at the current time step, really used as a scaling value for the error weights
    /*!
     *  We use this as an initial scaling vector for the solution tolerances
     *  We keep a copy of the previously converged solution around here to help set the value 
     *  of the solution weights for the next call to this nonlinear solver.
     */
    Epetra_Vector_Ghosted* m_y_pred_n;

    //! Value of the current solutionDot vector with ghost node values
    Epetra_Vector_Ghosted* m_ydot_curr;

    //! Value of the current solutionDot vector but only with owned values
    Epetra_Vector_Owned* m_ydot_curr_owned;

    //! Value of the proposed new solutionDot vector with ghost node values
    Epetra_Vector_Ghosted* m_ydot_new;

    //! Value of the proposed new solutionDot vector but only with owned values
    Epetra_Vector_Owned* m_ydot_new_owned;

    //! Value of the previous solutionDot vector from the previous time step
    /*!
     *  This is sometimes needed to update the solution and solutionDot vector for the current time step
     */
    Epetra_Vector_Ghosted* m_ydot_nm1;

    /*******************************************************************************
     *  METHOD OPTIONS FOR PRINTING ALGORITHM BEHAVIOR AND SOLUTION PROGRESS
     ********************************************************************************/

    //! Dump Jacobians to disk -  default false
    bool m_dumpJacobians;

    /*******************************************************************************
     *      COUNTERS FOR SOLUTION PROGRESS
     ********************************************************************************/

    /************************
     * TIME VARIABLES
     ************************/

    //! True if we are solving a time-dependent problem.
    //! False if we are solving a steady state problem.
    bool doTimeDependentResid_;

    //! Current value of the time
    double time_n;

    /**********************************************************************************
     * INTERNAL SOLUTION DATA VALUES
     ***********************************************************************************/

    //! Current value of the time step
    /*!
     *  This doesn't change within the routine.  This is set by setPrevoiusTimeStep
     */
    double delta_t_n;

    /*****************************************************************************************
     *        INTERNAL RESIDUAL VALUES
     ******************************************************************************************/

    //!  Residual value.
    /*!
     *   This may either be scaled or unscaled depending on the value of m_resid_scaled
     */
    Epetra_Vector_Owned* m_resid;

    //! Boolean indicating whether the vector m_resid is currently scaled or not.
    bool m_resid_scaled;

    //!   RHS of the linear problem
    Epetra_Vector_Owned* m_rhs;

    //! Vector of residual weights
    /*!
     *   These are used to establish useful and informative weighted norms of the residual vector.
     */
    Epetra_Vector_Owned* m_residWts;

    /*****************************************************************************************
     *        INTERNAL WEIGHTS FOR TAKING SOLUTION NORMS
     ******************************************************************************************/

    //! Norm of the residual at the start of each nonlinear iteration
    double m_normResid0;

    //! Norm of the residual before damping
    double m_normResidFRaw;

    //! Norm of the solution update created by the iteration in its raw, undamped form.
    double m_normSolnFRaw;

    //! Norm of the residual for a trial calculation which may or may not be used
    double m_normResidTrial;

    //! Vector of the norm
    double m_normResidPoints[15];

    /*****************************************************************************************
     *        INTERNAL BOUNDARY INFO FOR SOLUTIONS
     *****************************************************************************************/

    //! Vector of solution lower bounds
    Epetra_Vector_Owned* solnLowBound_;

    //! Vector of solution high bounds
    Epetra_Vector_Owned* solnHighBound_;

    /*****************************************************************************************
     *        INTERNAL TIME VARIABLES
     ******************************************************************************************/

    //! Main step change variable
    /*!
     *  Epectra vector for owned unknowns
     */
    Epetra_Vector_Owned* m_stp;

    //! Tentative second step change - used to see if the first step change is a good step
    Epetra_Vector_Owned* m_step_2;

    //! Pointer to the Residual Function
    ProblemResidEval* m_func;

    //! Epetra_Vector pointer  for  Vector of row scaling factors -  owned unknowns
    Epetra_Vector_Owned* m_rowScales;

    //! Epetra_Vector pointer  for  Vector of column scaling factors -  owned unknowns
    Epetra_Vector_Owned* m_colScales;

    //! Boolean Vector indicating whether its algebraic
    Epetra_IntVector* m_isAlgebraic;

    //! Integer vector
    Epetra_IntVector* m_isArithmeticScaled;

    //! Damping factor for the bounds algorithm
    double m_fbound;

    //! Current value of the damping factor
    double m_fdamp;

    //! Pointer to the jacobian representing the time dependent problem
    EpetraJac* tdjac_ptr;

    //! Number of function evaluations
    int m_nfe;
    
    //! Number of Jacobian Evaluations and factorization steps (they are the same)
    int m_nJacEval;

    //! Number of newton its in the current solve
    int m_num_newt_its;

    //!  Number of total newton iterations
    int m_numTotalNewtIts;

    //! Total number of linear iterations
    int m_numTotalLinearSolves;

    //! Freqency of residual updates
    int m_frequencyResidWtUpdatesTD;

    //! Maximum value of the jacobian age (unused atm)
    int m_maxAge;

    /***********************************************************************************************
     *   INTERNAL COMMUNICATORS
     **********************************************************************************************/
    //! My processor number
    int mypid_;

public:
    //! Static variable that can be set to dump out Jacobians.
    /*!
     *   Dump the jacobian if true and print level is >= 7
     */
    static bool s_print_NumJac;
};
//==================================================================================================================================
} // namespace
//----------------------------------------------------------------------------------------------------------------------------------


#endif
