/**
 *  @file m1d_SolnNonlinear.h
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

namespace m1d
{
//==================================================================================================================================
/**
 * Exception class thrown when a BEuler error is encountered.
 */
class SolNonlinearErr : public m1d_Error {
public:
  SolNonlinearErr(std::string msg);
};

#define BEULER_JAC_ANAL 2
#define BEULER_JAC_NUM  1

/**
 *  Wrapper class for 'beuler' integrator
 *  We derive the class from the class Integrator
 *  Newton iterator for multi-domain, one-dimensional problems.
 *  Used by class OneDim.
 */
class SolNonlinear : public SolGlobalNonlinear {

public:
  /**
   * The default constructor doesn't take an argument.
   */
  SolNonlinear();

  //! Destructor
  virtual
  ~SolNonlinear();

  //! Setup the problem for solution.
  virtual  void
  setup_problem(Solve_Type_Enum solveType,
                Epetra_Vector_Ghosted* y_init,
                Epetra_Vector_Ghosted* ydot_init,
                double time_curr,
                ProblemResidEval &problem,
                EpetraJac& jac);

  //! Apply hard bounds on the step size
  /*!
   *  We apply this before other terms, by decreasing the size of the original
   *  step size.
   *
   *  @param fbound Factor that the step size had to be reduced by
   */
  virtual int
  doHardBounds(const Epetra_Vector_Ghosted &y_old, Epetra_Vector_Owned &step, double &fbound);

  //! Compute factor to keep all components in bounds.
   virtual double
  deltaBoundStep(const Epetra_Vector_Ghosted& y, const Epetra_Vector_Owned& step0);

  virtual double
  highLowBoundStep(const Epetra_Vector_Ghosted& y, const Epetra_Vector_Owned& step0, int loglevel);

  /// Set options.
  virtual  void
  setOptions(int maxJacAge = 5)
  {
    m_maxAge = maxJacAge;
  }

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

  //! Set the value of the maximum # of newton iterations
  /*!
   *  @param maxNewtIts   Maximum number of newton iterations
   *                      The default value of this is 50 iterations
   */
  void setMaxNewtIts(const int maxNewtIts);

  virtual void
  setProblemType(int probtype);

  virtual void
  setDefaultSolnWeights();

  virtual void
  setRowScaling(bool onoff);

  //! Toggle that turns on and off column scaling
  /*!
   * Columns scaling is turned on by default
   * @param onoff  toggle
   * @param colScaleUpdateFrequency column scale update frequency
   *
   *        0 never
   *        1 once at the start
   *        2 after every jac update
   */
   virtual void setColScaling(bool onoff, int colScaleUpdateFrequency);

  //! Set the toggles for solution damping
  /*!
   *
   * @param residSolnDamping
   * @param deltaDamping
   * @param highLowDamping
   */
   virtual void setDampingToggles(const bool residSolnDamping, const bool deltaDamping, const bool highLowDamping);

  //! Set the vectors for lower and upper boundaries.
  /*!
   *
   * @param lowBounds
   * @param highBounds
   */
  virtual void
  setSolutionBounds(const Epetra_Vector_Owned &lowBounds, const Epetra_Vector_Owned &highBounds);

  //! Print the largest contributors to the solution update
  //! on a relative term
  /*!
   *
   * @param soln0
   * @param s0
   * @param soln1
   * @param s1
   * @param title
   * @param y0
   * @param y1
   * @param damp
   * @param num_entries
   */
  void
  print_solnDelta_norm_contrib(const Epetra_Vector &soln0,
                               const char * const s0,
                               const Epetra_Vector & soln1,
                               const char * const s1,
                               const char * const title,
                               const Epetra_Vector &y0,
                               const Epetra_Vector &y1,
                               double damp,
                               int num_entries);

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
   *   @param dumpJacobians Dump jacobians to disk.
   */
  virtual void
  setPrintSolnOptions(bool dumpJacobians);
  
  //! Set some nonlinear solver options
  /*!
   *    @param min_newt_its          Set the minimum number of newton iterations to carry out at every step.
   *                                 A value of 2 is useful to monitor convergence.
   *    @param matrixConditioning    If true, carry out matrix conditioning process before trying to find the
   *                                 inverse of the matrix (Not implemented).
   *    @param colScaling            Implement column scaling to the matrix system. The unknowns are then scaled
   *                                 inversely to the matrix. Scaling is carried out so that the deltas at 
   *                                 convergence of the matrix are scaled to a value of one. This means that
   *                                 the column scales are equal to the weighting matrix. 
   *                                 Column scaling is carried out before row scaling.
   *                                 (default = off)
   *    @param rowScaling            Implement row scaling to the matrix system. The max row element in the
   *                                 jacobian is caled to one by multiplying all terms in the row and rhs
   *                                 by a scale factor. Note, this is effective in practise in reducing
   *                                 the condition number of the matrix.
   *                                 (default = on).
   */
  void
  setNonLinOptions(int min_newt_its = 0, bool matrixConditioning = false, bool colScaling = false, bool rowScaling =
      true, int colScalingUpdateFrequency = 1);


  //!
  virtual void
  setPredicted_soln(const Epetra_Vector &y_pred);

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
  soln_error_norm(const Epetra_Vector_Owned &delta_y,
                  const bool printLargest = false,
                  const char *title = 0,
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
  res_error_norm(const Epetra_Vector_Owned &resid, const char *title = 0, const int printLargest = 0) const;

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
  get_jac(EpetraJac &jac,
          Epetra_Vector_Owned *res,
          const bool doTimeDependentResid,
          double time_curr,
          double rdelta_t,
          const Epetra_Vector_Ghosted *solnBase_ptr,
          const Epetra_Vector_Ghosted *solnDotBase_ptr);

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
          const Epetra_Vector_Ghosted *solnBase_ptr,
          const Epetra_Vector_Ghosted *solnDotBase_ptr);

public:

  //! Main routine to launch a nonlinear solve at the current conditions
  //! whether it's a steady state solve or a time-dependent run.
  virtual int
  solve_nonlinear_problem(Solve_Type_Enum solveType,
                          Epetra_Vector_Ghosted* y_comm,
                          Epetra_Vector_Ghosted* ydot_comm,
                          double CJ,
                          double time_curr,
                          int &num_newt_its,
                          int &num_linear_solves,
                          int &num_backtracks);

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
  scaleMatrix(Epetra_Vector_Owned &delta_soln,
              const Epetra_Vector_Ghosted &y_curr,
              const Epetra_Vector_Ghosted &ydot_curr,
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
  doNewtonSolve(Epetra_Vector_Owned &delta_soln, const Epetra_Vector &y_curr, const Epetra_Vector &ydot_curr, 
                const double time_curr, const double rdelta_t, int loglevel);

  //! Attempt to find a damping step
  /*!
   * On entry, step0 must contain an undamped Newton step for the
   * solution x0. This method attempts to find a damping coefficient
   * such that the next undamped step would have a norm smaller than
   * that of step0. If successful, the new solution after taking the
   * damped step is returned in y1, and the undamped step at y1 is
   * returned in step1.
   *
   * @param time_curr
   * @param y0
   * @param ydot0
   * @param step0
   * @param y1
   * @param ydot1
   * @param step1
   * @param s1
   * @param r
   * @param jac
   * @param loglevel
   * @param num_backtracks
   * @return
   *           1 Successful step was taken: Next step was less than previous step.
   *                                        s1 is calculated
   *           2 Successful step: Next step's norm is less than 0.8
   *           3 Success:  The final residual is less than 1.0
   *                        A predicted deltaSoln is not produced however. s1 is estimated.
   *           4 Success:  The final residual is less than the residual
   *                       from the previous step.
   *                        A predicted deltaSoln is not produced however. s1 is estimated.
   *           0 Uncertain Success: s1 is about the same as s0
   *          -2 Unsuccessful step.
   */
  int
  dampStep(double time_curr,
           const Epetra_Vector& y0,
           const Epetra_Vector* ydot0_ptr,
           double &s1,
           int& loglevel,
           int& num_backtracks);

 int dampStep_alt(double time_curr, const Epetra_Vector& y0, const Epetra_Vector* ydot0_ptr, double &s1, int& loglevel, int& num_backtracks);


  //! Set the column scales used in the program
  /*!
   *
   * @param colScales
   */
  virtual void
  setColumnScaleVector(const Epetra_Vector_Owned &colScales);

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
  getColumnScaleVector(Epetra_Vector_Owned & colScales) const;

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
   *  row scale the matrix, after column scaling is used. This routine returns the residual weights.
   *
   * @param residWts  Returns the residual weights in a vector of length 
   */
  virtual void getResidWts(Epetra_Vector_Owned &residWts);

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
   *  This routine uses the member function showSolutionVector() of the
   *  ProblemResidEval class, which does the layout printing. It can be used to
   *  print to a file as well.
   *
   * @param header   Header that is printed out to 
   * @param v        Epetra vector -> Note only the owned elements are needed.
   *                     (mp behavior hasn't been checked)
   */
  void print_IntSVector(std::string header, const Epetra_IntVector &v) const;

 
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
  calc_ydot(int order, const Epetra_Vector &y_curr, Epetra_Vector &ydot_curr);

private:
  //! Check to see if the nonlinear problem has converged
  /*!
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

  /**
   * m_jacFormMethod determines how a matrix is formed.
   */
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

  /**
   * m_matrixConditioning is a boolean. If true, then the
   * Jacobian and every rhs is multiplied by the inverse
   * of a matrix that is suppose to rfeduce the condition
   * number of the matrix. This is done before row scaling.
   */
  bool m_matrixConditioning;

  /**
   *  Relative nonlinear solver truncation error tolerances
   */
  double m_reltol;
  /**
   *  Vector of absolute nonlinear solver truncation error tolerance
   *  when not uniform for all variables.
   */
  Epetra_Vector_Owned *m_abstol;
  /**
   * Error Weights. This is a surprisingly important quantity.
   */
  Epetra_Vector_Owned *m_ewt;

  Epetra_Vector_Owned *m_ewt_deltaDamping;
  Epetra_Vector_Owned *m_absTol_deltaDamping;

  const Epetra_Comm *Comm_ptr_;

  /**
   * Current integration order
   */
  int m_order;

  /**
   * Failure Counter -> keeps track of the number
   * of consequetive failures
   */
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

  int m_curr_linearIts;

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

  /**
   * Number of equations in the ode integrator on the
   * current processor
   */
  int m_NumLcEqns;

  /**
   * Number of equations in the ode integrator on the
   * current processor that are owned by the current processor.
   */
  int m_NumLcOwnedEqns;

  //! Number of equations on all processors
  int m_NumGbEqns;

  /******************************************************************************
   *    METHOD OPTIONS FOR THE LINEAR SOLVER
   ********************************************************************************/

  Epetra_Vector_Ghosted *m_y_curr;
  Epetra_Vector_Owned   *m_y_curr_owned;

  Epetra_Vector_Ghosted *m_y_new;
  Epetra_Vector_Owned   *m_y_new_owned;

  //! Stored solution vector
  //Epetra_Vector_Ghosted *m_y_n;

  //! Stored solution vector at the previous time step
  Epetra_Vector_Ghosted *m_y_nm1;

  //! Predicted solution vector at the current time step
  Epetra_Vector_Ghosted *m_y_pred_n;

  Epetra_Vector_Ghosted *m_ydot_curr;
  Epetra_Vector_Owned   *m_ydot_curr_owned;

  Epetra_Vector_Ghosted *m_ydot_new;
  Epetra_Vector_Owned   *m_ydot_new_owned;

  Epetra_Vector_Ghosted *m_ydot_nm1;

  /*******************************************************************************
   *  METHOD OPTIONS FOR PRINTING ALGORITHM BEHAVIOR AND SOLUTION PROGRESS
   ********************************************************************************/

  /**
   * Dump Jacobians to disk
   *     default false
   */
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

  double time_n;
  //double timeStep;

  /**********************************************************************************
   * INTERNAL SOLUTION DATA VALUES
   ***********************************************************************************/

  //! Current value of the time step
  /*!
   *  This doesn't change within the routine
   *  This is set by setPrevoiusTimeStep
   */
  double delta_t_n;

  /*****************************************************************************************
   *        INTERNAL RESIDUAL VALUES
   ******************************************************************************************/

  //!  Residual value.
  /*!
   *   This may either be scaled or unscaled depending on the value of m_resid_scaled
   */
  Epetra_Vector_Owned *m_resid;

  //! boolean indicating whether the vector m_resid is scaled or not.
  /*!
   *
   */
  bool m_resid_scaled;

  //!   RHS of the linear problem
  Epetra_Vector_Owned *m_rhs;

  //! Vector of residual weights
  /*!
   *   These are used to establish useful and informative weighted norms of the residual vector.
   */
  Epetra_Vector_Owned *m_residWts;

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
   *
   */
  Epetra_Vector_Owned *m_stp;

  //! Tentative second step change - used to see if the first step change is a good step
  Epetra_Vector_Owned *m_step_2;

  //! Pointer to the Residual Function
  ProblemResidEval *m_func;

  Epetra_Vector_Owned *m_rowScales;
  Epetra_Vector_Owned *m_colScales;

  //! Boolean Vector indicating whether its algebraic
  Epetra_IntVector *m_isAlgebraic;

  //! Integer vector 
  Epetra_IntVector *m_isArithmeticScaled;

  double m_fbound;
  //! Current value of the damping factor
  double m_fdamp;
  

  /**
   * Pointer to the jacobian representing the
   * time dependent problem
   */
  EpetraJac *tdjac_ptr;

  /**
   * Number of function evaluations
   */
  int m_nfe;
  /**
   * Number of Jacobian Evaluations and
   * factorization steps (they are the same)
   */
  int m_nJacEval;

  //! Number of newton its in the current solve
  int m_num_newt_its;

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

  //! Freqency of residual updates 
  int m_frequencyResidWtUpdatesTD;

  int m_maxAge;
  int m_nv, m_np, m_n;
  doublereal m_elapsed;

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

} // namespace


#endif
