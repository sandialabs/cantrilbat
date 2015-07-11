/*
 * $Id: Electrode_defs.h 571 2013-03-26 16:44:21Z hkmoffa $
 */

#ifndef CT_ELECTRODE_DEFS_H
#define CT_ELECTRODE_DEFS_H

namespace Cantera
{


const int cBaseType = 1;
const int cInfCapacity= 2;
const int cMP_RxnExtent = 3;
const int cMultiPlateau_NoDiff = 4;
const int cSimpleDiff = 5;
const int cSimplePhaseChangeDiffusion = 6;
const int cCSTR = 7;
const int cCSTR_MCMBAnode = 8;
const int cCSTR_LiCoO2Cathode = 9;
const int cDiffTALE = 10;

//! Description of what the capacity means in terms of direction, and the
//! specification of whether this electrode will be used as an anode or a
//! cathode in the calculation.
/*!
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
 *   the following enum variable.
 */
enum Electrode_Capacity_Type_Enum {
    //!  Capacity and state of charge concepts are pertinent for an anode
    CAPACITY_ANODE_ECT = 0,
    //!  Capacity and state of charge concepts are pertinent for a cathode.
    CAPACITY_CATHODE_ECT,
    //!  Capacity that is neither anodic or cathodic.
    CAPACITY_OTHER_ECT
};

enum Electrode_Function_Type_Enum {
    //! Anode
    ELECTRODE_ANODE_EFT = 0,
    //! Cathode
    ELECTRODE_CATHODE_EFT,
    //! Reference electrode
    ELECTRODE_REFERENCE_EFT
    
};


};

#endif

