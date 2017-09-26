/*
 * $Id: Electrode_SimpleDiffCVode.h 571 2013-03-26 16:44:21Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_SIMPLEDIFF_H
#define _ELECTRODE_SIMPLEDIFF_H


#include "Electrode_Integrator.h"

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
class Electrode_SimpleDiff : public Electrode_Integrator
{
public:
    //! Constructor
    Electrode_SimpleDiff();

    //! Destructor
    virtual  ~Electrode_SimpleDiff();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    Electrode_SimpleDiff(const Electrode_SimpleDiff& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    Electrode_SimpleDiff& operator=(const Electrode_SimpleDiff& right);

    //! create the electrode model
    int electrode_model_create(ELECTRODE_KEY_INPUT* ei);

    void calcRate(double deltaT);

    void extractInfo(std::vector<size_t>& justBornMultiSpecies);


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
     *  @param iph           Print the phase
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values
     */
    virtual void printElectrodePhase(size_t iph, int pSrc = 1, bool subTimeStep = false) override;

    // ---------------------- SURFACE AREAS -------------------------------------------------------

    //! The internal state of the electrode must be kept for the initial and
    //! final times of an integration step.
    /*
     *  This function advances the initial state to the final state that was calculated
     *  in the last integration step.
     *
     *  @param Tinitial   This is the New initial time. This time is compared against the "old"
     *                   final time, to see if there is any problem.
     *  @param[in]           doAdvancementAlways Always do the reset, no matter what. Normally, Tinitial is checked against the 
     *                                           current t_init_init value. If they are the same, then we redo the time step.
     *                                           However, if  doAdvancementAlways is true, we advance the solution unknowns to the 
     *                                           final_final values produced in the last global step no matter what.
     *
     *  @return                                  Returns true if the time step is reset to t_init_init.
     */
    virtual bool resetStartingCondition(double Tinitial, bool doAdvancementAlways = false) override;

    //! Return the number of equations in the equation system that is used to solve the ODE integration
    /*!
     *  @return                                  Returns the number of equations
     */
    int nEquations() const;



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

    int index_Li_int_;

    int pindex_Li7_int;


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

    Cantera::Integrator* m_integ;

    //! Absolute tolerance for nonlinear residual
    std::vector<double> atolBaseResid_;

    //! Relative tolerance for nonlinear residual
    std::vector<double> rtolBaseResid_;


};

}


#endif
/*****************************************************************************/
