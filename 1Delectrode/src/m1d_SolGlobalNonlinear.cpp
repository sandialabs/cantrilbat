/**
 *
 *  @file m1d_SolNonlinear.cpp
 *
 *  Damped Newton solver for 1D multi-domain problems
 */

/*
 *  $Author: hkmoffa $
 *  $Date: 2013-01-07 15:32:48 -0700 (Mon, 07 Jan 2013) $
 *  $Revision: 504 $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_defs.h"

#include "m1d_Comm.h"
#include "m1d_SolNonlinear.h"
#include "cantera/base/clockWC.h"
#include "mdp_allo.h"

#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_DataAccess.h"

#include "m1d_LocalNodeIndices.h"
#include <stdio.h>
#include <math.h>


using namespace std;

#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y)    (( (x) < (y) ) ? (x) : (y))
#endif

namespace m1d
{

//-----------------------------------------------------------
//                 Constants
//-----------------------------------------------------------

const double DampFactor = 4;
const int NDAMP = 10;

//-----------------------------------------------------------
//                 Static Functions
//-----------------------------------------------------------



  class errBC : public m1d_Error {
  public:
    errBC(std::string procedure) :
      m1d_Error(procedure, "Base Class SolGlobalNonlinear called")
    {
    }
  };


//=====================================================================================================================
SolGlobalNonlinear::SolGlobalNonlinear() :
 m_print_flag(3)
{
}
//=====================================================================================================================
SolGlobalNonlinear::~SolGlobalNonlinear() {

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
 */
double
SolGlobalNonlinear::soln_error_norm(const Epetra_Vector &delta_y,
                              const bool printLargest,
                              const char *title,
                              const int typeYsoln,
                              const double dampFactor) const
{
  throw errBC("soln_error_norm()");
  return 0.0;
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
double
SolGlobalNonlinear::res_error_norm(const Epetra_Vector_Owned &resid, const char *title, const int printLargest) const
{
 
  throw errBC("res_error_norm()");
  return 0.0;
}
//=====================================================================================================================
void
SolGlobalNonlinear::get_res(const double time_curr,
                      const double rdelta_t,
                      const Epetra_Vector_Ghosted *solnBase_ptr,
                      const Epetra_Vector_Ghosted *solnDotBase_ptr)
{
 throw errBC("scaleMatrix()");
}

//=====================================================================================================================
void
SolGlobalNonlinear::scaleMatrix(Epetra_Vector_Owned &delta_soln,
                          const Epetra_Vector_Ghosted &y_curr,
                          const Epetra_Vector_Ghosted &ydot_curr,
                          const double time_curr,
                          const double rdelta_t,
                          int loglevel)
{
  throw errBC("scaleMatrix()");
}
//=====================================================================================================================
//  Compute the undamped Newton step.
/*
 * The residual function is
 * evaluated at the current time, time_curr, at the current values of the
 * solution vector, y_curr, and the solution time derivative, ydot_curr,
 * but the Jacobian is not recomputed.
 */
void
SolGlobalNonlinear::doNewtonSolve(Epetra_Vector_Owned &delta_soln,
                            const Epetra_Vector_Ghosted &y_curr,
                            const Epetra_Vector_Ghosted &ydot_curr,
                            const double time_curr,
                            const double rdelta_t,
                            int loglevel)
{
  throw errBC("doNewtonSolve()");
}
//=====================================================================================================================

int
SolGlobalNonlinear::doHardBounds(const Epetra_Vector_Ghosted &y_old, Epetra_Vector_Owned &step, double &fbound)
{
  throw errBC("doHardBounds()");
  return 1;
}
//=====================================================================================================================
double
SolGlobalNonlinear::deltaBoundStep(const Epetra_Vector_Ghosted& y, const Epetra_Vector_Owned& step0)
{
  throw errBC("deltaBoundStep()");
  return 0.0;
}
//=====================================================================================================================
/*
 *
 * boundStep():
 *
 * Return the factor by which the undamped Newton step 'step0'
 * must be multiplied in order to keep all solution components in
 * all domains between their specified lower and upper bounds.
 * *
 * This routine is meant to be used with cropping. In other words,
 * each component is allowed to go slightly out of bounds. However,
 * cropping is used to enforce a strict limit.
 *
 * Currently the bounds are hard coded into this routine:
 *
 *  Minimum value for all variables: solnLowBound[i] - 0.01 * m_ewt[i]
 *  Maximum value = none. solnHighBound[i] + 0.01 * m_ewt[i]
 *
 * Thus, this means that all solution components are expected
 * to be numerical greater than zero in the limit of time step
 * truncation errors going to zero.
 *
 * Delta bounds: The idea behind these is that the Jacobian
 *               couldn't possibly be representative, if the
 *               variable is changed by a lot. (true for
 *               nonlinear systems, false for linear systems)
 *  Maximum increase in variable in any one newton iteration:
 *   factor of 2
 *  Maximum decrease in variable in any one newton iteration:
 *   factor of 5
 */
double
SolGlobalNonlinear::highLowBoundStep(const Epetra_Vector_Ghosted& y, const Epetra_Vector_Owned& step0, int loglevel)
{
  throw errBC("highLowBoundStep()");
  return 1.0;
}

//=====================================================================================================================
/*
 * setColumnScales():
 *
 * Set the column scaling vector at the current time
 */
void
SolGlobalNonlinear::setDefaultColumnScaleVector()
{
  throw errBC("setDefaultColumnScaleVector()");
}
//=====================================================================================================================
// Return the column scales
/*
 *   Note, if there are no column scaling, then 1's are returned in the vector.
 *
 * @param colScales
 */
bool
SolGlobalNonlinear::getColumnScaleVector(Epetra_Vector_Owned & colScales) const
{
 throw  errBC("getColumnScaleVector()");
 return true;
}
//=====================================================================================================================
// Set the column scales
/*
 *   Note, if there are no column scaling, then 1's are returned in the vector.
 *
 * @param colScales
 */
void
SolGlobalNonlinear::setColumnScaleVector(const Epetra_Vector_Owned & colScales)
{
 throw errBC("setColumnScaleVector()");
}

//=====================================================================================================================
// Setup the problem for solution.
void
SolGlobalNonlinear::setup_problem(Solve_Type_Enum solveType,
				  Epetra_Vector_Ghosted* y_init,
				  Epetra_Vector_Ghosted* ydot_init,
				  double time_curr,
				  ProblemResidEval &problem,
				  EpetraJac& jac)
{
  throw errBC("setup_problem()");
}
//=====================================================================================================================
/*
 *
 */
void
SolGlobalNonlinear::setPredicted_soln(const Epetra_Vector &y_pred)
{
  throw errBC("setPredicted_soln()");
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
int
SolGlobalNonlinear::solve_nonlinear_problem(Solve_Type_Enum solveType,
					    Epetra_Vector_Ghosted *y_comm,
					    Epetra_Vector_Ghosted *ydot_comm,
					    double CJ,
					    double time_curr,
					    int &num_newt_its_comm,
					    int &num_linear_solves,
					    int &num_backtracks)
{
  throw errBC("solve_nonlinear_problem()");
  return 1;
}
//=====================================================================================================================
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
void
SolGlobalNonlinear::calc_ydot(int order, const Epetra_Vector &y_curr, Epetra_Vector &ydot_curr)
{
  throw errBC("calc_ydot()");
}
//=====================================================================================================================
void
SolGlobalNonlinear::setTolerances(double reltol, int n, const double * const abstol)
{
  throw errBC("setTolerances()");
}
//=====================================================================================================================
void
SolGlobalNonlinear::setTolerances_deltaDamping(double reltol_dd, int n, const double * const abstol_dd)
{
  throw errBC("setTolerances_deltaDamping()");
}

//=====================================================================================================================
  // Set the value of the maximum # of newton iterations
  /*
   *  @param maxNewtIts   Maximum number of newton iterations
   */
  void SolGlobalNonlinear:: setMaxNewtIts(const int maxNewtIts) {
    throw errBC("setMaxnewtIts()");
  }
//=====================================================================================================================
void
SolGlobalNonlinear::setProblemType(int jacFormMethod)
{
  throw errBC("setProblemType()");
}
//=====================================================================================================================
void
SolGlobalNonlinear::setDefaultSolnWeights()
{
  throw errBC("setDefaultSolnWeights()");
}
//=====================================================================================================================
void
SolGlobalNonlinear::setPrintFlag(int print_flag)
{
  m_print_flag = print_flag;
}
//=====================================================================================================================
/*
 *
 * @param dumpJacobians Dump jacobians to disk.
 *
 *                   default = false
 */
void
SolGlobalNonlinear::setPrintSolnOptions(bool dumpJacobians)
{
  throw errBC("setPrintSolnOptions()");
}
//=====================================================================================================================
void
SolGlobalNonlinear::setNonLinOptions(int min_newt_its, bool matrixConditioning, bool colScaling, bool rowScaling,
                               int colScaleUpdateFrequency)
{
  throw errBC("setNonLinOptions()");
}
//=====================================================================================================================
// set the previous time step
/*
 *
 * @param jac
 */
void
SolGlobalNonlinear::setPreviousTimeStep(const double timeStep_comm, const Epetra_Vector& y_nm1, const Epetra_Vector& ydot_nm1)
{
  throw errBC("setPreviousTimeStep()");
}
//=====================================================================================================================
// Compute the Residual Weights
/*
 *  The residual weights are defined here to be equal to the inverse of the row scaling factors used to
 *  row scale the matrix, after column scaling is used. They are multiplied by 10-3 because the column
 *  weights are also multiplied by that same quantity.
 *
 *  The basic idea is that a change in the solution vector on the order of the convergence tolerance
 *  multiplied by  [RJC] which is of order one after row scaling should give you the relative weight
 *  of the row. Values of the residual for that row can then be normalized by the value of this weight.
 *  When the tolerance in delta x is achieved, the tolerance in the residual is also achieved.
 */
void
SolGlobalNonlinear::computeResidWts()
{ 
  throw errBC("computeResidWts()");
 
}
//=====================================================================================================================
void
SolGlobalNonlinear::getResidWts(Epetra_Vector_Owned &residWts)
{
  throw errBC("getResidWts()");
}
//=====================================================================================================================
void
SolGlobalNonlinear::setRowScaling(bool onoff)
{
  throw errBC("setRowScaling()");
}
//=====================================================================================================================
void
SolGlobalNonlinear::setColScaling(bool onoff, int colScaleUpdateFrequency)
{
  throw errBC("setColScaling()");
}
//=====================================================================================================================
// Set the toggles for solution damping
/*
 *
 * @param residSolnDamping
 * @param deltaDamping
 * @param highLowDamping
 */
void
SolGlobalNonlinear::setDampingToggles(const bool residSolnDamping, const bool deltaDamping, const bool highLowDamping)
{
  throw errBC("setDampingToggles()");
}
//=====================================================================================================================
// Set the vectors for lower and upper boundaries.
/*
 * @param lowBounds
 * @param highBounds
 */
void
SolGlobalNonlinear::setSolutionBounds(const Epetra_Vector_Owned &lowBounds, const Epetra_Vector_Owned &highBounds)
{
  throw errBC("setSolutionBounds()");
}
//=====================================================================================================================
}

