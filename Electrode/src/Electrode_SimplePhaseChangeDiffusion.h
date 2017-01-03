/*
 * $Id: Electrode_SimplePhaseChangeDiffusion.h 571 2013-03-26 16:44:21Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_SIMPLEPHASECHANGEDIFFUSION_H
#define _ELECTRODE_SIMPLEPHASECHANGEDIFFUSION_H

#include "Electrode.h"


#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{


//! This class is a derived class used to model phase - change electrodes
/*!
 * Complete problem statement
 *
 */
class Electrode_SimplePhaseChangeDiffusion : public Electrode
{
public:
    //! Constructor
    Electrode_SimplePhaseChangeDiffusion();

    //! Destructor
    virtual  ~Electrode_SimplePhaseChangeDiffusion();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    Electrode_SimplePhaseChangeDiffusion(const Electrode_SimplePhaseChangeDiffusion& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    Electrode_SimplePhaseChangeDiffusion& operator=(const Electrode_SimplePhaseChangeDiffusion& right);

    //! Return the type of electrode
    /*!
     *  Returns the enum type of the electrode. This is used in the factory routine.
     *
     *  @return Returns an enum type, called   Electrode_Types_Enum
     */
    virtual Electrode_Types_Enum electrodeType() const;

    //! create the electrode model
    virtual int electrode_model_create(ELECTRODE_KEY_INPUT* ei);

    //!  Set the electrode initial conditions from the input file.
    /*!
     *   (virtual from Electrode)
     *   (This is a serial virtual function or an overload function)
     *
     *    This is one of the most important routines. It sets up the initial conditions of the electrode
     *    from the input file. The electrode itself has been set up from a call to electrode_model_create().
     *    After the call to this routine, the electrode should be internally ready to be integrated and reacted.
     *    It takes its input from an ELECTRODE_KEY_INPUT object which specifies the setup of the electrode
     *    object and the initial state of that object.
     *
     *    The routine works like an onion initialization. The parent object is initialized before the
     *    child. This means the child object first calls the parent, before it does its own initializations.
     *
     *  @param[in]           ei                  ELECTRODE_KEY_INPUT pointer object
     *
     *  @return                                  Returns zero if successful, and -1 if not successful.
     */
    virtual int setInitialConditions(ELECTRODE_KEY_INPUT* ei) override;


    void calcRate(double deltaT);

    void extractInfo(std::vector<size_t>& justBornMultiSpecies);

    virtual void updateState();

    //! Print conditions of the electrode for the current integration step to stdout
    /*!
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values
     */
    virtual void printElectrode(int pSrc = 1, bool subTimeStep = false) override;

    //! Print condition of a phase in the electrode
    /*!
     *  @param iPhase        Print the phase
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values
     */
    virtual void printElectrodePhase(size_t iPhase, int pSrc = 1, bool subTimeStep = false) override;

    // ---------------------- SURFACE AREAS -------------------------------------------------------

    //! The internal state of the electrode must be kept for the initial and
    //! final times of an integration step.
    /*
     *  This function advances the initial state to the final state that was calculated
     *  in the last integration step.
     *
     * @param Tinitial   This is the New initial time. This time is compared against the "old"
     *                   final time, to see if there is any problem.
     */
    virtual void resetStartingCondition(double Tinitial, bool doTestsAlways = false) override;

    //!  Calculate the change in the state of the system when integrating from T_initial_initial_
    //!  to t_final_final_
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
     *  @param fieldInterpolationType Type of interpolation of field variables defaults to T_FINAL_CONST_FIS,
     *  @param subIntegrationType     Type of subintegration. Defaults to BASE_TIMEINTEGRATION_SIR.
     *                                In this integration, the program determines its own strategy
     *                                for the time step.
     *
     *  @return Returns the number of subcycle steps it took to complete the full step.
     *          Failures to complete the integration due to time truncation error issues return a -1.
     *          Failures due to invalid function calculation attempts return a -2.
     *          Failures due to invalid arguments return a -3.
     */
    virtual int integrate(double deltaT, double  GlobalRtolSrcTerm = 1.0E-3,
                          Electrode_Exterior_Field_Interpolation_Scheme_Enum fieldInterpolationType = T_FINAL_CONST_FIS,
                          Subgrid_Integration_RunType_Enum subIntegrationType = BASE_TIMEINTEGRATION_SIR);


    void check_initial_CAP();
    void check_final_CAP();
    void check_final_OuterVol();

protected:

    //! we identify the phases here as being an outer phase.
    /*!
     *  The outer phase will be in direct contact with the electrolyte.
     *  the outer phase will be in contact with surrounding particles. Therefore,
     *  the voltage of the outer phase will be considered the voltage of the electrode.
     */

    //! Phase index of the outer solid phase
    int phaseIndexOuterSolidPhase_;

    //! Molar density of the outer solid phase
    /*!
     *  units are kmol m-3.
     */
    double MD_OuterSolidPhase_;

    size_t index_Li_int_;

    size_t pindex_Li7_int;


    //! Phase index of the solid phase that comprises the inner solid
    int phaseIndexInnerSolidPhase_;
    double MD_InnerSolidPhase_;


    //! we identify the phases here as being an outer phase.
    /*!
     *  The outer phase will be in direct contact with the electrolyte.
     *  the outer phase will be in contact with surrounding particles. Therefore,
     *  the voltage of the outer phase will be considered the voltage of the electrode.
     */



    int surfIndexOuterSurface_;
    int surfIndexInnerSurface_;


    //!  Value of the flux at the outer edge - kmol m-2 s-1
    double Nflux_final_;
    //!  value of the rate from a single particle - kmol s-1;
    double Nrate_final_;

    //! calculated values
    double Da_s_;
    double k_r_internal_;
    double k_f_internal_;
    double k_r_external_;
    double k_f_external_;
    double a_lip_;

    double C_external_final_;
    double C_internal_final_;
    double mf_internal_final_;
    double mf_internal_init_;

    double C_external_init_;
    double C_external_init_init_;
    double C_internal_init_;
    double C_internal_init_init_;

    double mf_external_final_;
    double mf_external_init_;

    double MD_int_;

    double Radius_internal_final_;
    double Radius_internal_final_final_;
    double Radius_internal_eff_;
    double Radius_internal_init_;
    double Radius_internal_init_init_;


    double MN_internal_init_;
    double MN_internal_final_;

    double DiffCoeff_;

    double ROP_inner_;
    double ROP_rate_inner_;

    double ROP_outer_;
    double ROP_rate_outer_;

    bool   zeroedInnerRadius_;
    double deltaTZeroed_;
    bool   NoInnerSolid_;

    int    SolidInnerKSpecies_;
    double SolidInnerKSpeciesReacStoichCoeff_;

    int    SolidOuterKSpecies_;
    double SolidOuterKSpeciesProdStoichCoeff_;

    double CAP_init_;
    double CAP_final_;

    double deltaCAPdt_;

};

}


#endif
/*****************************************************************************/
