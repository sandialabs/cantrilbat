/**
 * @file Electrode_InfCapacity.h
 *   Definitions for an electrode object that has an infinite capacity
 *   (see \ref electrode_mgr and class \link Zuzax::Electrode_InfCapacity Electrode_InfCapacity \endlink).
 */

/*
 * Copywrite 2010 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_INFCAPACITY_H
#define _ELECTRODE_INFCAPACITY_H

#include "Electrode.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace Zuzax
{
//==================================================================================================================================
//!  Electrode with infinite capacity
/*!
 *            How the capacity of an infinite capacity electrode is calculated
 *
 *   The capacity is set by the setCapacityParams() routine. We use the optional parameter
 *   for that routine to set the variable capacitySpeciesSpMoles_. Then, the capacity is
 *   calculcated from the following formula:
 *
 *       capacity = capacitySpeciesSpMoles_ * capacitySpeciesCoeff_ * Faraday
 *
 */
class Electrode_InfCapacity : public Zuzax::Electrode
{
public:

    //! Constructor
    Electrode_InfCapacity();

    //! Destructor
    virtual ~Electrode_InfCapacity();

    //! Copy Constructor
    /*!
     * @param[in]            right               Object to be copied
     */
    Electrode_InfCapacity(const Electrode_InfCapacity& right);

    //! Assignment operator
    /*!
     *  @param[in]           right               object to be copied
     *
     *  @return                                  Returns a reference to the current object
     */
    Electrode_InfCapacity& operator=(const Electrode_InfCapacity& right);

    //! Duplicator function
    /*!
     *  Duplicate the current Electrode object, returning a base Electrode pointer
     *
     *  @return                                   Returns a duplicate of the current object as a base class pointer
     */
    virtual Electrode* duplMyselfAsElectrode() const override;

    //! Return the type of electrode
    /*!
     *  Returns the enum type of the electrode. This is used in the factory routine.
     *
     *  @return                                  Returns an enum type, called   Electrode_Types_Enum
     */
    virtual Electrode_Types_Enum electrodeType() const override;

    //!  Setup the electrode
    /*!
     *  @param[in]           ei                  ELECTRODE_KEY_INPUT pointer object
     *
     *  @return                                  Returns an int
     */
    virtual int electrode_model_create(ELECTRODE_KEY_INPUT* ei) override;

    //!  Calculate the change in the state of the system when integrating from T_initial_initial_  to t_final_final_
    /*!
     *  (virtual from Electrode)
     *  All information is kept internal within this routine. This may be done continuously
     *  and the solution is not updated.
     *
     *  Note the tolerance parameters refere to the nonlinear solves within the calculation. They do not refer to time step parameters.
     *
     *  @param[in]           deltaT              DeltaT for the integration step.
     *  @param[in]            GlobalRtolSrcTerm  Relative tolerance for the source term vector calcualted from
     *                                           this routine.    Defaults to 1.0E-3
     *
     *  @param[in]           fieldInterpolationType Type of interpolation of field variables defaults to T_FINAL_CONST_FIS,
     *  @param[in]           subIntegrationType     Type of subintegration. Defaults to BASE_TIMEINTEGRATION_SIR.
     *                                              In this integration, the program determines its own strategy for the time step.
     *
     *  @return                                  Returns the number of subcycle steps it took to complete the full step.
     *                                           Failures to complete the integration due to time truncation error issues return a -1.
     *                                           Failures due to invalid function calculation attempts return a -2.
     *                                           Failures due to invalid arguments return a -3.
     */
    virtual int integrate(double deltaT, double  GlobalRtolSrcTerm = 1.0E-3,
                          Electrode_Exterior_Field_Interpolation_Scheme_Enum fieldInterpolationType = T_FINAL_CONST_FIS,
                          Subgrid_Integration_RunType_Enum subIntegrationType = BASE_TIMEINTEGRATION_SIR) override;

    //! The internal state of the electrode must be kept for the initial and final times of an integration step.
    /*!
     *  This function advances the initial state to the final state that was calculated
     *  in the last integration step.
     *
     *  @param[in]           Tinitial            This is the new initial time. This time is compared against the "old"
     *                                           final time, to see if there is any problem. The two should be the same.
     *                                           If Tinitial is t_init_init, we redo the time step.
     *
     *  @param[in]           doAdvancementAlways Always do the advancement, no matter what. Normally, Tinitial is checked against the 
     *                                           current t_init_init value. If they are the same, then we redo the time step.
     *                                           However, if  doResetAlways is true, we advance the solution unknowns to the 
     *                                           final_final values produced in the last global step no matter what.
     *
     *  @return                                  Returns true if the time step is reset to t_init_init.
     */
    virtual bool resetStartingCondition(double Tinitial, bool doAdvancementAlways = false) override;

    //! Set the internal initial intermediate and initial global state from the internal final state
    /*!
     *  (virtual function)
     *
     *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
     *  routine as well.
     *
     * @param setInitInit   Boolean indicating whether you should set the init_init state as well
     */
    virtual void setInitStateFromFinal(bool setInitInit = false) override;

    //! Returns the integrated moles transfered for each phase in the electrode object over the time step
    /*!
     *   (virtual from Electrode.h)
     *
     *    Rewrite because there is no check on phase_mole changes.
     *
     *    @param  phaseMolesTransfered   Vector of moles transfered (length = number of total phases in the electrode object)
     *                                   units = kmol
     */
    virtual void getIntegratedPhaseMoleTransfer(double* const phaseMolesTransfered) override;

};
//====================================================================================================================
}
//======================================================================================================================
#endif
//======================================================================================================================
