/**
 *  @file m1d_FlatBatteryResidEval.h
 *  File that contains the description of a single m1d problem that
 *  will be solved.
 */

/*
 *  $Author: hkmoffa $
 *  $Date: 2013-03-08 16:35:51 -0700 (Fri, 08 Mar 2013) $
 *  $Revision: 564 $
 *
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef M1D_FLATBATTERYRESIDEVAL_H
#define M1D_FLATBATTERYRESIDEVAL_H

#include "m1d_ProblemResidEval.h"
#include "m1d_BoundaryCondition.h"
#include "m1d_SurDomain1D.h"
#include "m1d_BatteryResidEval.h"  


namespace m1d {

//!  Residual for 1D cell  battery evaluations
/**
 * 
 */
class FlatBatteryResidEval: public BatteryResidEval
{

public:

    //! Default constructor
    /*!
     *
     * @param atol   Absolute tolerance for calculation
     */
    FlatBatteryResidEval(double atol = 1.0e-13);

    //! Destructor
    /*!
     *
     */
    virtual ~FlatBatteryResidEval();

    //! Default copy constructor
    /*!
     *
     * @param r  Object to be copied
     */
    FlatBatteryResidEval(const FlatBatteryResidEval &r);

    //! Assignment operator
    /*!
     *
     * @param r  Object to be copied
     * @return   Returns a copy of the current problem
     */
    FlatBatteryResidEval &
    operator=(const FlatBatteryResidEval&r);



    //! Write out to a file or to standard output the current solution
    /*!
     *   These functions are affected by the print controls of the nonlinear solver
     *   and the time stepper.
     *
     *      ievent is a description of the event that caused this
     *      function to be called.
     *
     *      @param ievent  Event that's causing this routine to be called.
     *                     =  0 Initial conditions for a calculation
     *                     =  1 Completion of a successful intermediate time step.
     *                     =  2 Completion of a successful Final time or final calculation.
     *                     =  3 Completion of a successful Intermediate nonlinear step
     *                     = -1 unsuccessful time step that converged, but failed otherwise
     *                     = -2 unsuccessful nonlinear step.
     *      @param doTimeDependentResid  true if solving a time dependent problem
     *      @param time_current      Current time
     *      @param delta_t_n         Current value of delta_t
     *      @param istep             Current step number
     *      @param soln_n               Current value of the solution vector
     *      @param solnDot_n_ptr            Current value of the time deriv of the solution vector
     *      @param delta_t_np1       Suggested next delta t value (defaults to 0.0 if there isn't a
     *                               good guess for the next delta_t).
     */
    virtual void
    writeSolution(const int ievent, const bool doTimeDependentResid, const double time_current, const double delta_t_n,
                  const int istep, const Epetra_Vector_Ghosted &soln_n, const Epetra_Vector_Ghosted * const solnDot_n_ptr,
                  const Solve_Type_Enum solveType = TimeDependentAccurate_Solve, 
	          const double delta_t_np1 = 0.0);

    //! This function may be used to create output at various points in the
    //! execution of an application.
    /*!
     *   These functions are not affected by the print controls of the nonlinear solver
     *   and the time stepper.
     *
     *      ievent is a description of the event that caused this
     *      function to be called.
     *
     *      @param ievent  Event that's causing this routine to be called.
     *                     =  0 Initial conditions for a calculation
     *                     =  1 Completion of a successful intermediate time step.
     *                     =  2 Completion of a successful Final time or final calculation.
     *                     =  3 Completion of a successful Intermediate nonlinear step
     *                     =  4 Write out current and voltage to timeDep_IV.dat
     *                     = -1 unsuccessful time step that converged, but failed otherwise
     *                     = -2 unsuccessful nonlinear step.
     *
     *      @param time_current      Current time
     *      @param delta_t_n         Current value of delta_t
     *      @param istep             Current step number
     *      @param y_n               Current value of the solution vector
     *      @param ydot_n_ptr        Current value of the time deriv of the solution vector
     */
    virtual void
    user_out(const int ievent, const double time_current, const double delta_t_n, const int istep, const Epetra_Vector_Ghosted &y_n,
             const Epetra_Vector_Ghosted * const ydot_n_ptr);

    void write_IV(const int ievent, const bool doTimeDependentResid, const double time_current, const double delta_t_n, int istep,
                  const Epetra_Vector_Ghosted &y_n, const Epetra_Vector_Ghosted * const ydot_n_ptr);

    

};
// *****************************************************************
//  end of m1d namespace
}
// ******************************************************************

#endif

