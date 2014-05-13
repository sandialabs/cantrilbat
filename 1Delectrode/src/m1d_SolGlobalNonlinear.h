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

namespace m1d
{


enum BEulerMethodType {
  BEulerFixedStep = 0, 
  BEulerVarStep
};

/**
 * Exception class thrown when a BEuler error is encountered.
 */
class SolGlobalNonlinearErr : public m1d_Error {
public:
  SolGlobalNonlinearErr(std::string msg);
};

#define BEULER_JAC_ANAL 2
#define BEULER_JAC_NUM  1

/**
 *  Wrapper class for 'beuler' integrator
 *  We derive the class from the class Integrator
 *  Newton iterator for multi-domain, one-dimensional problems.
 *  Used by class OneDim.
 */
class SolGlobalNonlinear {

public:
  /**
   * The default constructor doesn't take an argument.
   */
  SolGlobalNonlinear();

  //! destructor
  virtual
  ~SolGlobalNonlinear();

  //! Setup the problem for solution.
  virtual void
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
  virtual void
  setOptions(int maxJacAge = 5)
  {
   
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
  virtual void setMaxNewtIts(const int maxNewtIts);

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
   *        0 never
   *        1 once at the start
   *        2 after every jac update
   */
  virtual void
  setColScaling(bool onoff, int colScaleUpdateFrequency);

  //! Set the toggles for solution damping
  /*!
   *
   * @param residSolnDamping
   * @param deltaDamping
   * @param highLowDamping
   */
  virtual void
  setDampingToggles(const bool residSolnDamping, const bool deltaDamping, const bool highLowDamping);

  //! Set the vectors for lower and upper boundaries.
  /*!
   *
   * @param lowBounds
   * @param highBounds
   */
  virtual void
  setSolutionBounds(const Epetra_Vector_Owned &lowBounds, const Epetra_Vector_Owned &highBounds);

 
  virtual void
  setPrintSolnOptions(bool dumpJacobians);

  virtual void
  setNonLinOptions(int min_newt_its = 0, bool matrixConditioning = false, bool colScaling = false, bool rowScaling =
		   true, int colScalingUpdateFrequency = 1);

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
  virtual void
  setPrintFlag(int print_flag);

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
  soln_error_norm(const Epetra_Vector &delta_y,
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
  doNewtonSolve(Epetra_Vector_Owned &delta_soln,
                const Epetra_Vector &y_curr,
                const Epetra_Vector &ydot_curr,
                const double time_curr,
                const double rdelta_t,
                int loglevel);



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
  getResidWts(Epetra_Vector_Owned &residWts);



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

} // namespace


#endif
