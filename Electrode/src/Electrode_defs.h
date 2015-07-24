/*
 * $Id: Electrode_defs.h 571 2013-03-26 16:44:21Z hkmoffa $
 */

#ifndef CT_ELECTRODE_DEFS_H
#define CT_ELECTRODE_DEFS_H

namespace Cantera
{

/*!
 *    \addtogroup electrodeobj
 *    @{
 */

//!  Enum Type identifying the models for the Electrodes
/*!
 *  Note that this Enum class may be extended in other contexts.
 *  Therefore, loops over this enum should use a case statement with a default that either throws
 *  a catchable signal or returns a similar type of error message.
 */
enum Electrode_Types_Enum {
    UNKNOWN_ET = -1,

    //! Base Electrode Type with string name "BaseType"
    /*!
     *  This type is not a user type. In the future this class may turn into a virtual class.
     */
    BASE_TYPE_ET = 0,

    //! Infinite Capacity Electrode type with string name "InfCapacity"
    INF_CAPACITY_ET,

    //! Multi-plateau model employing a reaction extent mode, string name "MP_RxnExtent"
    MP_RXNEXTENT_ET,

    //! Multi-plateau model employing a full chemistry model, string name "MultiPlateau_NoDiff"
    MULTIPLATEAU_NODIFF_ET,

    //! Simple diffusion model, employing chemistry at interface, string name "SimpleDiff"
    SIMPLE_DIFF_ET,

    //! Diffusion model, based on the TALE algorithm, string named "DiffTALE"
    DIFF_TALE_ET,

    //! Two %Phase algorithm combined with steady state diffusion, string named "SimplePhaseChangeDiffusion"
    SIMPLE_PHASE_CHANGE_DIFFUSION_ET,

    //! CSTR model with reactions on the surface of the particle, no diffusion, string named "CSTR"
    CSTR_ET,

    //! CSTR mode, defined in a child class for zn metal anodes, string named "CSTR_Zn_Anode"
    CSTR_ZN_ANODE_ET,

    //! CSTR model, defined in a child class for MCMB anodes, string named, "CSTR_MCMBAnode"
    CSTR_MCMB_ANODE_ET,

    //! CSTR model, defined in a child class for LiCoO2, string named, "CSTR_LiCoO2Cathode"
    CSTR_LICO2_CATHODE_ET,

    //! Early model where successive substitution was used to relax the CSTR equations, string named "SuccessiveSubstitution"
    SUCCESSIVE_SUBSTITUTION_ET,

    //! Multi-plateau model employing a reaction extent mode, refined for FeS2 cathodes, string name "MP_RxnExtent_FeS2"
    MP_RXNEXTENT_FES2_ET,

    //! Multi-plateau model employing a reaction extent mode, refined for LiSi anodes, string name "MP_RxnExtent_LiSi"
    MP_RXNEXTENT_LISI_ET,

    //! Network of one or more radial diffusion regions with surrounding interfacial kinetics regions, string name "Radial_Diff_Regions"
    RADIAL_DIFF_REGIONS_ET

};
/*!    @}                       */


//! Description of what the capacity means in terms of direction, and the specification of 
//! whether this electrode will be used as an anode or a cathode in the calculation.
/*!
 *   \ingroup electrodeobj
 *
 *   The electrode workings are remarkedly devoid of the need to specify whether
 *   the particular model is a model for an anode or for a cathode.
 *   Basically, if the overpotential is positive i.e., phiMetal is greater than
 *   phiElectrode, the electrode will be producing electrons and therefore it will
 *   be acting as an anode. IF the overpotential is negative it will be consuming
 *   electrons and therefore the electrode will be acting as a cathode.
 *
 *   However, for the purposes of calculating the capacity and for calculating
 *   the state of charge, there needs to be a specification of the sign for
 *   whether has a state of charge of zero when it can't accept any more electrons
 *   (the anode situation) or for when it can't give any more electrons
 *   (the cathode situation). This is a sign convention on the concept of
 *   the current capacity of the electrode. This sign convention is given by
 *   the following enum variable. The CAPACITY_OTHER_ECT value is currently not used.
 */
enum Electrode_Capacity_Type_Enum {
    //!  Capacity and state of charge concepts pertinent for an anode
    CAPACITY_ANODE_ECT = 0,
    //!  Capacity and state of charge concepts pertinent for a cathode.
    CAPACITY_CATHODE_ECT,
    //!  Capacity that is neither anodic or cathodic.
    CAPACITY_OTHER_ECT
};

};

#endif

