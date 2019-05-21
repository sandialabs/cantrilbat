/**
 *  @file Electrod_Equilibrium.h Predict the stability of a single phase that is currently zereod.
 *
 *  This routine fills in an estimate for the solution
 *  Return 1 if the phases are stable and 0 if they are not
 *
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_EQUILIBRIUM_H
#define _ELECTRODE_EQUILIBRIUM_H

#include "Electrode.h"

#include <string>
#include <vector>
//----------------------------------------------------------------------------------------------------------------------------------
namespace Zuzax
{

class MP_EquilStatic;
class FixedChemPotSSTP;
//==================================================================================================================================
//! Class which determines the stability of phases due to kinetics
/*!
 *  Note, I believe I can make this class simpler
 *
 */
class Electrode_Equilibrium {

public:
    //! Constructor
    /*!
     *
     * @param elect  Electrode object pertaining to the stability problem
     *               The electrode object assigns this object as a "friend"
     */
    Electrode_Equilibrium(Electrode* elect);

    //! Destructor
    virtual ~Electrode_Equilibrium();

    //! Copy Constructor
    /*!
     *  @param[in]           right               Object to be copied
     */
    Electrode_Equilibrium(const Electrode_Equilibrium& right);

    //! Assignment operator
    /*!
     *  @param[in]           right               Object to be copied
     *
     *  @return                                  Returns a reference to the object being copied
     */
    Electrode_Equilibrium& operator=(const Electrode_Equilibrium& right);

    //! Set up the equilibrium problem
    /*!
     *   Sets up an MP_EquilStatic object with the PhaseList information from the electrode object
     *
     *  @return                                  Returns 0
     */
    int setupEquilibriumProblem();


    //! Add a Li ion bath capability
    /*!
     *  Adds a fixed chemical potential Li phase to the MP_EquilStatic object created by setupEquilibriumProblem().
     *  Using this approach, we can establish a bath phase of Li ions at a particular voltage.
     *
     *  @return                                  Returns 0
     */
    int addLiIonPlusElectronBathPhase();

    //! Set the state of the captive ThermoPhase objects and the MultiPhase object from the current
    //! final state of this Electrode
    void update_MP();

    //! Set the state of the electrode object, from this current MultiPhase object
    void uploadMP();

    //! Predict the stability of a single phase that is currently zeroed.
    /*!
     *  This routine fills in an estimate for the solution. It assumes that there is a reference
     *  electrode
     *
     *  @param[in] iphase     Phase id in the list of phases
     *  @param[in] funcStab   value of the stability function calculated by this routine
     *
     *  @return   Return 1 if the phases are stable and 0 if they are not
     */
    int predictStabilitySinglePhase(int iphase, double& funcStab);

    //! Returns the mole fractions of a phase just calculated within the routine.
    /*!
     *
     *  @param[in]           iphase      Phase id in the list of phases within the electrode object
     *  @param[in]           moleFractions
     */
    void getPhaseMoleFractions(int iphase, double* const moleFractions);

    //! Return a complete multiphase object representing equilibrium for the electrode
    /*!
     *  @return                                  Returns a multiphase object
     */
    MP_EquilStatic* MultiPhase_Obj();

protected:

    //! This is a reference to the friend object, where we will pull most of the data from.
    Electrode* ee_;

    //! MultiPhase object that is used in the equilibrium solves
    /*!
     *    Note, this object contains pointers to ThermoPhase objects. These objects are shared
     *    with the electrode object
     *
     */
    MP_EquilStatic* m_mp;

    //! Mapping between MultiPhase objects ThermoPhase objects to the Electrode objects
    /*!
     *   A value of -1 indicates the ThermoPhase doesn't exist in the Electrode object
     *
     *   PhaseIndex_mp_[MultiPhaseIndex] = ElectrodeVolumeIndex
     */
    std::vector<int> PhaseIndex_mp_;

    //! Fixed elemental chemical potential for Lithium that depends on the voltage of the electrode
    FixedChemPotSSTP* LiFixed_;

public:
    //! Print level for input to vcs routines within this object
    /*!
     *    This is a public member so that it can be manipulated
     */
    int printLvl_;

};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

