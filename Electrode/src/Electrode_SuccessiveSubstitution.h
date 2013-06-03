/*
 * $Id: Electrode.h 533 2013-02-21 19:50:50Z vebruni $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "Electrode.h"

#ifndef _ELECTRODE_SUCCESSIVESUBSTITUTION_H
#define _ELECTRODE_SUCCESSIVESUBSTITUTION_H



namespace Cantera
{

// -----------------------------------------------------------------------------------------------------------------

class Electrode_SuccessiveSubstitution : public Electrode
{
public:

    //! Constructor
    Electrode_SuccessiveSubstitution();

    //! Destructor
    virtual ~Electrode_SuccessiveSubstitution();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    Electrode_SuccessiveSubstitution(const Electrode_SuccessiveSubstitution& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    Electrode_SuccessiveSubstitution& operator=(const Electrode_SuccessiveSubstitution& right);

    //! Duplicator function
    /*!
     *  Duplicate the current electrode object, returning a base electrode pointer
     *
     *  (Virtual function from Electrode.h)
     */
    Electrode* duplMyselfAsElectrode() const;



    Electrode_Types_Enum electrodeType() const;


    //!  Calculate the change in the state of the system when integrating from Tglobal_initial
    //!  to Tglobal_final
    /*!
     *  All information is kept internal within this routine. This may be done continuously
     *  and the solution is not updated.
     *
     *  Note the tolerance parameters refere to the nonlinear solves within the calculation
     *  They do not refer to time step parameters.
     *
     *  @param deltaT        DeltaT for the integration step.
     *  @param GlobalRtolSrcTerm Relative tolerance for the source term vector calcualted from
     *                       this routine.
     *                       Defaults to 1.0E-3
     *  @param fieldInterpolationType Type of interpolation of field variables defaults to 0
     *  @param subIntegrationType     Type of subintegration. defaults to 0
     *
     *  @return Returns the number of subcycle steps it took to complete the full step.
     */
    virtual int integrate(double deltaT, double  GlobalRtolSrcTerm = 1.0E-3,
                          Electrode_Exterior_Field_Interpolation_Scheme_Enum fieldInterpolationType = T_FINAL_CONST_FIS,
                          Subgrid_Integration_RunType_Enum subIntegrationType = BASE_TIMEINTEGRATION_SIR);


    //!  Residual calculation for the solution of the Nonlinear integration problem
    /*!
     *    Given tfinal and delta_t, and given y and ydot which are estimates of the solution
     *    and solution derivative at tfinal, this function calculates the residual equations.
     *    It is the residual function used in the nonlinear solver that relaxes the equations
     *    at each time step.
     *
     *    This is typically called from evalResidNJ(), which is called directly from the
     *    nonlinear solver. However, we expose this routine so that the residual can be queried
     *    given all of the inputs.
     *
     * @param tfinal        Time                    (input)
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
     * @return
     */
    virtual int integrateResid(const doublereal tfinal, const doublereal delta_t,
                               const doublereal* const y, const doublereal* const ydot,
                               doublereal* const resid,
                               const ResidEval_Type_Enum evalType, const int id_x, const doublereal delta_x);



    // -----------------------------------------------------------------------------------------------------------------
    // ---------------------------- SOLUTION OF NONLINEAR TIME DEPENDENT SYSTEM  ---------------------------------------
    // -----------------------------------------------------------------------------------------------------------------

    //! Return the number of equations in the Nonlinear equation system used to solve the system
    //! at every time step
    /*!
     *  This is also equal to the number of state variables in the problem
     */
    virtual int nEquations() const;


};
// -----------------------------------------------------------------------------------------------------------------


}

#endif
/*****************************************************************************/
