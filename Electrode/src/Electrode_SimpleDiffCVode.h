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



#include "Electrode.h"
#include "Electrode_Integrator.h"
#include "cantera/integrators.h"
#include "cantera/numerics/ResidJacEval.h"

namespace Cantera
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
    void extractInfo(std::vector<int>& justBornMultiSpecies);


    //! Print conditions of the electrode for the current integration step to stdout
    /*!
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values
     */
    virtual void printElectrode(int pSrc = 1, bool subTimeStep = false);

    //! Print condition of a phase in the electrode
    /*!
     *  @param iPhase        Print the phase
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values
     */
    virtual void printElectrodePhase(int iPhase, int pSrc = 1, bool subTimeStep = false);

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
    virtual void  resetStartingCondition(double Tinitial, bool doTestsAlways = false);


    //! Return the number of equations in the equation system that is used to solve the ODE integration
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
    doublereal MD_OuterSolidPhase_;

    int index_Li_int_;

    int pindex_Li7_int;


    //! Phase index of the solid phase that comprises the inner solid
    int phaseIndexInnerSolidPhase_;
    doublereal MD_InnerSolidPhase_;


    //! we identify the phases here as being an outer phase.
    /*!
     *  The outer phase will be in direct contact with the electrolyte.
     *  the outer phase will be in contact with surrounding particles. Therefore,
     *  the voltage of the outer phase will be considered the voltage of the electrode.
     */



    int surfIndexOuterSurface_;
    int surfIndexInnerSurface_;


    //!  Value of the flux at the outer edge - kmol m-2 s-1
    doublereal Nflux_final_;
    //!  value of the rate from a single particle - kmol s-1;
    doublereal Nrate_final_;

    //! calculated values
    doublereal Da_s_;
    doublereal k_r_internal_;
    doublereal k_f_internal_;
    doublereal k_r_external_;
    doublereal k_f_external_;
    doublereal a_lip_;

    doublereal C_external_final_;
    doublereal C_internal_final_;
    doublereal mf_internal_final_;
    doublereal mf_internal_init_;

    doublereal C_external_init_;
    doublereal C_external_init_init_;
    doublereal C_internal_init_;
    doublereal C_internal_init_init_;

    doublereal mf_external_final_;
    doublereal mf_external_init_;

    doublereal MD_int_;

    doublereal Radius_internal_final_;
    doublereal Radius_internal_final_final_;
    doublereal Radius_internal_eff_;
    doublereal Radius_internal_init_;
    doublereal Radius_internal_init_init_;


    doublereal MN_internal_init_;
    doublereal MN_internal_final_;

    doublereal DiffCoeff_;

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
    std::vector<doublereal> atolBaseResid_;

    //! Relative tolerance for nonlinear residual
    std::vector<doublereal> rtolBaseResid_;


};

}


#endif
/*****************************************************************************/
