/**
 *  @file BEulerInt_Battery.h
 */


/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */
#ifndef CT_BEULERINT_BATTERY_H
#define CT_BEULERINT_BATTERY_H

#include "m1d_BEulerInt.h"
#include "m1d_SurDomain1D.h"
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
class BoundaryCondition;
}
//----------------------------------------------------------------------------------------------------------------------------------
namespace beuler
{
//=================================================================================================================================
//! Integrator for Battery applications
/*!
 *  the main need for this layer is a reworking of the initial condition problem for constant current boundary conditions.
 */
class BEulerInt_Battery : public beuler::BEulerInt
{

public:

    //!The default constructor doesn't take an argument.
    /*!
     *  Default settings: epetra jacobian, no user-supplied
     *  Jacobian function, Newton iteration.
     *  */
    BEulerInt_Battery();

    //! Copy constructor
    /*!
     *
     * @param r  Object to be copied
     */
    BEulerInt_Battery(const BEulerInt_Battery& r);

    //! Virtual destructor
    virtual
    ~BEulerInt_Battery();

    //! Assignment operator
    /*!
     *  @param r  Object to be copied
     *  @return  Returns the current object
     */
    BEulerInt_Battery& operator=(const BEulerInt_Battery& r);

    //! Solve for the consistent initial conditions and consistent initial time derivatives.
    /*!
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
     *  This top level routine does a trick. For problems with a constant current boundary condition. It
     *  first solves the consistent initial condition problem using a close constant voltage problem.
     *  Then, it changes the problem back to the original problem and resolves the consistent initial
     *  condition problem. The reason is that the constant current problem fails often when given bad
     *  initial guesses.
     *
     *  @return                                  0  Always returns 0
     */
    virtual int calcConsistentInitialDerivs() override;

    //! Solve for the consistent initial conditions and consistent initial time derivatives.
    /*!
     *  This is the inner calculation. It gets called by calcConsistentInitialDerivs() to do all of  the work
     *
     *  @return                                  0  Always returns 0
     */
    int calcConsistentInitialDerivs_Inner();

    //!  Check to see that the predicted solution satisfies proper requirements.
    /*!
     *  @return                                  Returns a negative number if the step is inappropriate.
     *                                           Then the stepsize is reduced and the method is checked again.
     */
    virtual int check_predicted_soln(m1d::Epetra_Vector_Ghosted& y_n, m1d::Epetra_Vector_Ghosted& ydot_n,
                                     double CJ, double time_n) override;

    /*************************************** Member data ***********************************************/

    //! This is the type of constant current boundary condition that this object will solve
    int BC_Type_current_;

    //! This is the value of the double used in the specification of the Cathode's current boundary condition
    double CurrentValueBC_;

    //! Value of the derivative of the cathode voltage
    double cathodeVoltageDot_;

    //! Boundary condition function
    m1d::BoundaryCondition* BC_TimeDep_;

    //! Time dependency function
    m1d::TimeDepFunctionPtr TimeDep_;

    //! Value of the current (amps)
    double currentNeeded_;

    //! Best current value of the cathode voltage
    double CathodeVoltageBest_;
};
//=================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

