/**
 *  @file m1d_BatteryResidEval.h
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

#ifndef M1D_BATTERYRESIDEVAL_H
#define M1D_BATTERYRESIDEVAL_H

#include "m1d_ProblemResidEval.h"
#include "m1d_BoundaryCondition.h"
#include "m1d_SurDomain1D.h"

namespace m1d {

//!  Residual for 1D cell  battery evaluations
/**
 * 
 */
class BatteryResidEval: public ProblemResidEval
{

public:

    //! Default constructor
    /*!
     *
     * @param atol   Absolute tolerance for calculation
     */
    BatteryResidEval(double atol = 1.0e-13);

    //! Destructor
    /*!
     *
     */
    virtual ~BatteryResidEval();

    //! Default copy constructor
    /*!
     *
     * @param r  Object to be copied
     */
    BatteryResidEval(const BatteryResidEval &r);

    //! Assignment operator
    /*!
     *
     * @param r  Object to be copied
     * @return   Returns a copy of the current problem
     */
    BatteryResidEval &
    operator=(const BatteryResidEval&r);

    //! Calculate the initial conditions
    /*!
     *   This calls the parent class initialConditions method to loop over the volume and surface domains.
     *   Then the method tried to better estimate the electrolyte potential.
     *
     * @param doTimeDependentResid    Boolean indicating whether we should
     *                                formulate the time dependent residual
     * @param soln                    Solution vector. This is the input to
     *                                the residual calculation.
     * @param solnDot                 Solution vector. This is the input to
     *                                the residual calculation.
     * @param t                       Time
     * @param delta_t                 delta_t for the initial time step
     */
    void
    initialConditions(const bool doTimeDependentResid, Epetra_Vector_Ghosted *soln, Epetra_Vector_Ghosted *solnDot, const double t,
                      const double delta_t);

    //! Make an attempt to improve the initial guess for the electrolyte voltage  based on the boundary conditions
    void
    improveInitialConditions(Epetra_Vector_Ghosted *soln);


    //! Calculate a residual vector
    /*!
     *   The basic algorithm is to loop over the volume domains.
     *   Then, we loop over the surface domains
     *
     * @param res                     residual output
     * @param doTimeDependentResid    Boolean indicating whether the time
     *                                dependent residual is requested
     * @param doTimeDependentResid    Boolean indicating whether we should
     *                                formulate the time dependent residual
     * @param soln                    Pointer to the solution vector. This is the input to the residual calculation.
     * @param solnDot                 Pointer to the solution Dot vector. This is the input to the residual calculation.
     * @param t                       current time
     * @param rdelta_t                delta t inverse
     * @param residType               Residual type
     * @param solveType               Solve type
     */
    virtual void
    residEval(Epetra_Vector* const &  res,
	      const bool doTimeDependentResid,
	      const Epetra_Vector_Ghosted *soln,
	      const Epetra_Vector_Ghosted *solnDot,
	      const double t,
	      const double rdelta_t,
	      const ResidEval_Type_Enum residType = Base_ResidEval,
	      const Solve_Type_Enum solveType = TimeDependentAccurate_Solve);

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
     */
    virtual void
    writeSolution(const int ievent, const bool doTimeDependentResid, const double time_current, const double delta_t_n,
                  const int istep, const Epetra_Vector_Ghosted &soln_n, const Epetra_Vector_Ghosted * const solnDot_n_ptr,
                  const Solve_Type_Enum solveType = TimeDependentAccurate_Solve);

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

    //! Set a solution parameter 
    /*!
     *  @param paramName   String identifying the parameter to be set
     *  @param paramVal    Single double value of the parameter to be set
     *
     *  @return returns the number of parameters set by the command, zero or a positive number.
     *          Returns a negative number if the parameter is unknown
     */
    virtual int
    setSolutionParam(std::string paramName, double paramVal);

    //! Get a solution parameter 
    /*!
     *  @param paramName   String identifying the parameter to be set
     *  @param paramVal    Vector of parameters returned
     *
     *  @return returns the number of parameters returned.
     */
    virtual int
    getSolutionParam(std::string paramName, double * const paramVal);

    //!   Report on the boundary condition applied to the cathode voltage equation
    /*!
     *
     *   @param[in]  time        Current time for evaluating time dependent BC
     *   @param[out] BC_Type     Type of the boundary condition
     *   @param[out] value       Value of the dirichlet condition or flux - default 0.0
     *   @param[out] BC_TimeDep  BoundaryCondition Pointers for time dependent BC for BC_Tppe = 3,4,5
     *                            (default 0)
     *   @param[out] TimeDep     Function pointer to a function that returns a double given a single parameter (the time).
     *                           Defaults to a NULL pointer.
     */
    void
    reportCathodeVoltageBC(double time, int &BC_Type, double &value, BoundaryCondition * &BC_TimeDep,
                           TimeDepFunctionPtr &TimeDep) const;

    //!   Change the boundary condition applied to the cathode voltage equation
    /*!
     *   @param[in] BC_Type     Type of the boundary condition
     *   @param[in] value       Value of the Dirichlet condition or flux - default 0.0
     *   @param[in] BC_TimeDep  BoundaryCondition Pointers for time dependent BC for BC_Type = 3,4,5
     *                          Defaults to a NULL pointer.
     *   @param[in] TimeDep     Function pointer to a function that returns a double given a single parameter (the time).
     *                          Defaults to a NULL pointer.
     */
    void
    changeCathodeVoltageBC(int BC_Type, double value, BoundaryCondition * BC_TimeDep = 0, TimeDepFunctionPtr TimeDep = 0);

    double reportCathodeVoltage() const;

  
    double reportCathodeCurrent() const;

    //! Get the max value of the sub grid time step number from the last residual calculation
    /*!
     *   @return   Returns the max subgrid time step number from all of the electrodes. Special steps
     *             aren't counted in this number.
     */
    int getMaxSubGridTimeSteps() const;

 protected:
    int maxSubGridTimeSteps_;

};

// *****************************************************************
//  end of m1d namespace
}
// ******************************************************************

#endif

