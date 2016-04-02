/*
 * @File GFCEO_Electrode_Integrator.h coupling 
 */
/*
 * Copywrite 2015 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _GFCEO_ELECTRODE_H
#define _GFCEO_ELECTRODE_H


#include "Electrode.h"
#include "cantera/numerics/ResidJacEval.h"
//-----------------------------------------------------------------------------------------------------------------------------------
namespace Cantera
{

//===================================================================================================================================
//! This class is a derived class used to carry out fully coupled simulations
/*!
 * Complete problem statement
 *
 */
class GFCEO_Electrode : public Cantera::ResidJacEval
{
public:

    // ---------------------------------------------------------------------------------------------
    // ----------------------- BASIC SETUP ROUTINES  -----------------------------------------------
    // ---------------------------------------------------------------------------------------------

    //! Constructor
    /*!
     *  @param[in]      atol                     Default value for the absolute tolerance
     */
    GFCEO_Electrode(Electrode* ee, doublereal atol = 1.0E-13, int iOwn = 0);

    //! Destructor
    virtual ~GFCEO_Electrode();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    GFCEO_Electrode(const GFCEO_Electrode& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    GFCEO_Electrode& operator=(const GFCEO_Electrode& right);

    //! Returns a changeable reference to the electrode object
    /*!
     *  @return                         Returns a changeable reference to the electrode object
     */
    Electrode& electrode();


  //! Return the number of equations in the equation system
    /*!
     *  @return                                        Returns the number of equations in the equation system
     */
    virtual int nEquations() const;

   //! Constrain solution component k. 
    /*!
     *
     *  @param[in]        k                 index of the unknown
     *  @param[in]        flag              type of the constraint
     *
     *        Possible values for  'flag' are:
     *   - c_NONE       no constraint
     *   - c_GE_ZERO    >= 0
     *   - c_GT_ZERO    >  0
     *   - c_LE_ZERO    <= 0
     *   - c_LT_ZERO    <  0
     */
    virtual void constrain(const int k, const int flag);

    //! Value for the constraint condition for unknown k
    /*!
     *  @param[in]          k                  index of the unknown
     *
     *  @return                                returns the value of the constaint condition for soln component
     */
    int constraint(const int k) const;

    //! Initialization function that sets the sizes in the object
    virtual void initSizes();

    //! Specify that solution component k is purely algebraic 
    /*!
     *    That is, the derivative of this component does not appear in the residual function.
     *
     *  @param[in]       k                    Solution index
     */
    virtual void setAlgebraic(const int k);


    //! Evaluate the residual function at the current conditions
    /*!
     * @param t             Time                    (input)
     * @param delta_t       The current value of the time step (input)
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input, do not modify)
     * @param resid         Value of the residual that is computed (output)
     * @param evalType      Type of the residual being computed (defaults to Base_ResidEval)
     * @param id_x          Index of the variable that is being numerically differenced to find
     *                      the jacobian (defaults to -1, which indicates that no variable is being
     *                      differenced or that the residual doesn't take this issue into account)
     * @param delta_x       Value of the delta used in the numerical differencing
     *
     * @return Returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *           -0 or neg value Means an unsuccessful operation
     */
    virtual int evalResidNJ(const doublereal t, const doublereal delta_t,
                            const doublereal* const y,
                            const doublereal* const ydot,
                            doublereal* const resid,
                            const ResidEval_Type_Enum evalType = Base_ResidEval,
                            const int id_x = -1,
                            const doublereal delta_x = 0.0);

   //! Evaluate the residual function at the current conditions
    /*!
     *  This is a wrapper around the evalResidNJ function, which should be used instead because 
     *  delta_t is not supplied in this interface.
     *  delta_t is assigned to a value of -1, and evalResidNJ() is then called.
     *
     * @param[in]         t             Time                    (input)
     * @param[in]         y             Solution vector (input, do not modify)
     * @param[in]         ydot          Rate of change of solution vector. (input, do not modify)
     * @param[out]        resid         Value of the residual that is computed (output)
     *
     * @return Returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *           -0 or neg value Means an unsuccessful operation
     */
    virtual int eval(const doublereal t, const doublereal* const y, const doublereal* const ydot,
                     doublereal* const resid);

    //! Input the initial conditions for the problem
    /*!
     * Values for both the solution and the value of ydot may be provided.
     *
     * @param[in] t0             Time
     * @param[out] y             Solution vector
     * @param[out] ydot          Rate of change of solution vector.
     *
     * @return                   Returns a flag to indicate that operation is successful.
     *                              1  Means a successful operation
     *                             -0 or neg value Means an unsuccessful operation
     */
    virtual int getInitialConditions(const doublereal t0, doublereal* const y, doublereal* const ydot);


  //! Filter the solution predictions
    /*!
     * Codes might provide a predicted step change. This routine filters the predicted
     * solution vector eliminating illegal directions.
     *
     * @param t             Time                    (input)
     * @param ybase         Solution vector (input, output)
     * @param step          Proposed step in the solution that will be cropped
     *
     * @return              Return the norm of the amount of filtering
     */
    virtual doublereal filterNewStep(const doublereal t, const doublereal* const ybase,
                                     doublereal* const step);

    //! Filter the solution predictions
    /*!
     * Codes might provide a predicted solution vector. This routine filters the predicted
     * solution vector.
     *
     * @param t             Time                    (input)
     * @param y             Solution vector (input, output)
     *
     * @return              Return the norm of the amount of filtering
     */
    virtual doublereal filterSolnPrediction(const doublereal t, doublereal* const y);

    //! Set a global value of the absolute tolerance
    /*!
     *  @param atol   Value of atol
     */
    void setAtol(doublereal atol);

    //! Evaluate the time tracking equations, if any
    /*!
     * Evaluate time integrated quantities that are calculated at the
     * end of every successful time step.  This call is made once at the end of every successful
     * time step that advances the time. It's also made once at the start of the time stepping.
     *
     * @param t             Time                    (input)
     * @param delta_t       The current value of the time step (input)
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input, do not modify)
     *
     * @return Returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *           -0 or neg value Means an unsuccessful operation
     */
    virtual int  evalTimeTrackingEqns(const doublereal t, const doublereal delta_t, const doublereal* const y,
                                      const doublereal* const ydot);

    //! Evaluate any stopping criteria other than a final time limit
    /*!
     *  If we are to stop the time integration for any reason other than reaching a final time limit, tout,
     *  provide a test here. This call is made at the end of every successful time step iteration
     *
     *  @return    If true, the the time stepping is stopped. If false, then time stepping is stopped if t >= tout
     *             Defaults to false.
     *
     * @param t             Time                    (input)
     * @param delta_t       The current value of the time step (input)
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input, do not modify)
     */
    virtual bool evalStoppingCritera(const doublereal t,
                                     const doublereal delta_t,
                                     const doublereal* const y,
                                     const doublereal* const ydot);

    //! Return a vector of delta y's for calculation of the numerical Jacobian
    /*!
     *   There is a default algorithm provided.
     *
     *        delta_y[i] = atol[i] + 1.0E-6 ysoln[i]
     *        delta_y[i] = atol[i] + MAX(1.0E-6 ysoln[i] * 0.01 * solnWeights[i])
     *
     * @param t             Time                    (input)
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input, do not modify)
     * @param delta_y       Value of the delta to be used in calculating the numerical jacobian
     * @param solnWeights   Value of the solution weights that are used in determining convergence (default = 0)
     *
     * @return Returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *           -0 or neg value Means an unsuccessful operation
     */
    virtual int
    calcDeltaSolnVariables(const doublereal t,
                           const doublereal* const y,
                           const doublereal* const ydot,
                           doublereal* const delta_y,
                           const doublereal* const solnWeights = 0);

    //!  Returns a vector of column scale factors that can be used to column scale Jacobians.
    /*!
     *  Default to yScales[] = 1.0
     *
     * @param t             Time                    (input)
     * @param y             Solution vector (input, do not modify)
     * @param y_old         Old Solution vector (input, do not modify)
     * @param yScales       Value of the column scales
     */
    virtual void calcSolnScales(const doublereal t, const doublereal* const y,
                                const doublereal* const y_old, doublereal* const yScales);


 //! This function may be used to create output at various points in the execution of an application.
    /*!
     *  @param ifunc     identity of the call
     *                          0  Initial call
     *                          1  Called at the end of every successful time step
     *                         -1  Called at the end of every unsuccessful time step
     *                          2  Called at the end of every call to integrateRJE()
     *
     * @param t             Time                    (input)
     * @param delta_t       The current value of the time step (input)
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input)
     */
    virtual void user_out2(const int ifunc, const doublereal t,
                           const doublereal delta_t,
                           const doublereal* const y,
                           const doublereal* const ydot);

    //! This function may be used to create output at various points in the execution of an application.
    /*!
     *  This routine calls user_out2().
     *
     *  @param ifunc     identity of the call
     * @param t             Time                    (input)
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input)
     */
    virtual void user_out(const int ifunc, const doublereal t,
                          const doublereal* y,
                          const doublereal* ydot);

    //! Multiply the matrix by another matrix that leads to better conditioning
    /*!
     *  Provide a left sided matrix that will multiply the current jacobian, after scaling
     *  and lead to a better conditioned system.
     *  This routine is called just before the matrix is factored.
     *
     *  Original Problem:
     *        J delta_x = - Resid
     *
     *  New problem:
     *      M (J delta_x) = - M Resid
     *
     *  @param    matrix     Pointer to the current jacobian (if zero, it's already been factored)
     *  @param    nrows      offsets for the matrix
     *  @param    rhs        residual vector. This also needs to be lhs multiplied by M
     *
     * @return Returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *           -0 or neg value Means an unsuccessful operation
     */
    virtual int  matrixConditioning(doublereal* const matrix, const int nrows,
                                    doublereal* const rhs);

    //! Calculate an analytical jacobian and the residual at the current time and values.
    /*!
     *  Only called if the jacFormation method is set to analytical
     *
     * @param t             Time                    (input)
     * @param delta_t       The current value of the time step (input)
     * @param cj            Coefficient of yprime used in the evaluation of the jacobian
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input, do not modify)
     * @param J             Reference to the SquareMatrix object to be calculated (output)
     * @param resid         Value of the residual that is computed (output)
     *
     * @return Returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *           -0 or neg value Means an unsuccessful operation
     */
    virtual int evalJacobian(const doublereal t, const doublereal delta_t, doublereal cj,
                             const doublereal* const y, const doublereal* const ydot,
                             GeneralMatrix& J, doublereal* const resid);

    //! Calculate an analytical jacobian and the residual at the current time and values.
    /*!
     *  Only called if the jacFormation method is set to analytical
     *
     * @param t             Time                    (input)
     * @param delta_t       The current value of the time step (input)
     * @param cj            Coefficient of yprime used in the evaluation of the jacobian
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input, do not modify)
     * @param jacobianColPts   Pointer  to the vector of pts to columns of the SquareMatrix
     *                         object to be calculated (output)
     * @param resid         Value of the residual that is computed (output)
     *
     * @return Returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *           -0 or neg value Means an unsuccessful operation
     */
    virtual int evalJacobianDP(const doublereal t, const doublereal delta_t, doublereal cj,
                               const doublereal* const y,
                               const doublereal* const ydot,
                               doublereal* const* jacobianColPts,
                               doublereal* const resid);

     
 //!      Write out to a file or to standard output the current solution
    /*!
     *  @param[in]     ievent                  a description of the event that caused this function to be called.
     *  @param[in]     time                    Time
     *  @param[in]     deltaT                  Delta Time for the last time step
     *  @param[in]     time_step_num           time step number
     *  @param[in]     y                       Solution vector
     *  @param[out]    ydot                    Rate of change of solution vector.
     * 
     */
    virtual void writeSolution(int ievent, const double time, const double deltaT,  const int time_step_num,
                               const double* const y, const double* const ydot);




    //! Change ownership of Electrode object to this object if it is not already
    void assertOwnership();    

protected:

    //! Electrode object
    Electrode* ee_;

    //! This object owns the Electrode object and is responsible for deletion.
    bool iOwnObject_;

};

} // End namespace Cantera
//-----------------------------------------------------------------------------------------------------------------------------------
#endif

