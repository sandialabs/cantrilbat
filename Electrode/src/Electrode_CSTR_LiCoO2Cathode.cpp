/**
 *  @file Electrode_CSTR_LiCoO2Cathode.cpp
 *     Definitions for the Electrode_CSTR class, used to model 
 *     Electrode processes in particles with no transport limits hard-coded for LiCoO2 Cathodes
 *     (see \ref electrode_mgr and class \link Zuzax::Electrode_CSTR_LiCoO2Cathode Electrode_CSTR_LiCoO2Cathode\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "Electrode_CSTR_LiCoO2Cathode.h"

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
/*
 *  ELECTRODE_INPUT: constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.
 */
Electrode_CSTR_LiCoO2Cathode::Electrode_CSTR_LiCoO2Cathode() :
    Electrode_CSTR(),
    Global_LiCoO2_Model_(3),
    ig_SolidLi_(npos),
    ig_SolidV_(npos),
    ip_LiCoO2_(npos)
{
}
//==================================================================================================================================
Electrode_CSTR_LiCoO2Cathode::Electrode_CSTR_LiCoO2Cathode(const Electrode_CSTR_LiCoO2Cathode& right) :
    Electrode_CSTR_LiCoO2Cathode()
{
    operator=(right);
}
//==================================================================================================================================
Electrode_CSTR_LiCoO2Cathode& Electrode_CSTR_LiCoO2Cathode::operator=(const Electrode_CSTR_LiCoO2Cathode& right)
{
    if (this == &right) {
        return *this;
    }
    Electrode_CSTR::operator=(right);

    Global_LiCoO2_Model_ = right.Global_LiCoO2_Model_;
    ig_SolidLi_          = right.ig_SolidLi_;
    ig_SolidV_           = right.ig_SolidV_;
    ip_LiCoO2_           = right.ip_LiCoO2_;
    return *this;
}
//==================================================================================================================================
Electrode_CSTR_LiCoO2Cathode::~Electrode_CSTR_LiCoO2Cathode()
{
}
//==================================================================================================================================
Electrode_Types_Enum  Electrode_CSTR_LiCoO2Cathode::electrodeType() const
{
    return CSTR_LICO2_CATHODE_ET;
}
//==================================================================================================================================
void  Electrode_CSTR_LiCoO2Cathode::setCapacityCoeff_LiCoO2()
{
    double  capacityLeftSpeciesCoeff = 1.0;
    for (size_t iph = 0; iph < m_NumVolPhases; iph++) {
        if (iph == solnPhase_ || iph == metalPhase_) {
            continue;
        }
        size_t kStart = m_PhaseSpeciesStartIndex[iph];
        std::string pname = PhaseNames_[iph];
        bool found = false;
        for (size_t k = 0; k < NumSpeciesList_[iph]; k++) {
            size_t iGlobSpeciesIndex = kStart + k;
            std::string sss = speciesName(iGlobSpeciesIndex);
            found = false;
            capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] =0.0;
            capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = 0.0;
            if (sss == "CoO2") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex]    = capacityLeftSpeciesCoeff;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = capacityLeftSpeciesCoeff;
                found = true;
            }

            if (sss == "LiCoO2") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex]    = 0.0;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = capacityLeftSpeciesCoeff;
                found = true;
            }
            if (!found) {
                throw Electrode_Error(":setCapacityCoeff_LiCoO2()", "unknown species: " + sss + " in phase " + pname);
            }
        }
    }
}
//==================================================================================================================================
int Electrode_CSTR_LiCoO2Cathode::electrode_model_create(ELECTRODE_KEY_INPUT* ei)
{
    int flag = Electrode_CSTR::electrode_model_create(ei);
    if (flag != 0) {
        return flag;
    }
    electrodeType_ = ELECTRODETYPE_CATHODE;

    ig_SolidV_ = globalSpeciesIndex("CoO2","LiCoO2_Interstitials_cathode");
    if (ig_SolidV_ == npos) {
        ig_SolidV_ = globalSpeciesIndex("CoO2","LiCoO2_Margules_1");
        if (ig_SolidV_ == npos) {
            throw Electrode_Error("Electrode_CSTR_LiCoO2Cathode::electrode_model_create()",
                               "Expected to find species CoO2 in phase LiCoO2_Interstitials_cathode");
        } else {
            Global_LiCoO2_Model_ = 1;
        }
    } else {
        Global_LiCoO2_Model_ = 3;
    }
    ig_SolidLi_ = npos;
    if (Global_LiCoO2_Model_ == 3) {
        ig_SolidLi_ = globalSpeciesIndex("LiCoO2", "LiCoO2_Interstitials_cathode");
    } else if (Global_LiCoO2_Model_ == 1) {
        ig_SolidLi_ = globalSpeciesIndex("LiCoO2", "LiCoO2_Margules_1");
    }
    if (ig_SolidLi_ == npos) {
        throw ZuzaxError("Electrode_CSTR_LiCoO2Cathode::electrode_model_create()",
                         "Expected to find species CoO2 in phase LiCoO2_Interstitials_cathode");
    }
    if (Global_LiCoO2_Model_ == 3) {
        ip_LiCoO2_ = globalPhaseIndex("LiCoO2_Interstitials_cathode");
    } else  if (Global_LiCoO2_Model_ == 1) {
        ip_LiCoO2_ = globalPhaseIndex("LiCoO2_Margules_1");
    } else {
        ip_LiCoO2_ = globalPhaseIndex("LiCoO2_Interstitials_cathode");
    }
    if (ip_LiCoO2_ == npos) {
        throw Electrode_Error("Electrode_CSTR_MCMBAnode::electrode_model_create()",
                           "I tried to find the phase LiCoO2_Interstitials_cathode  but failed. May need to generalize the code now");
    }

    RelativeExtentRxn_RegionBoundaries_.resize(2);
    RelativeExtentRxn_RegionBoundaries_[1] = 0.94;
    RelativeExtentRxn_RegionBoundaries_[0] = 0.55;
    double li_mf = moleFraction(ig_SolidLi_);

    RelativeExtentRxn_init_ = li_mf;
    RelativeExtentRxn_final_ = RelativeExtentRxn_init_;
    RelativeExtentRxn_final_final_ = RelativeExtentRxn_init_;
    RelativeExtentRxn_init_init_ = RelativeExtentRxn_init_;
    if (RelativeExtentRxn_final_ <  RelativeExtentRxn_RegionBoundaries_[0]) {
        throw Electrode_Error("Electrode_CSTR_LiCoO2Cathode::electrode_model_create()", "Relative Extent Rxn outside of bounds");
    }
    if (RelativeExtentRxn_final_ >  RelativeExtentRxn_RegionBoundaries_[1]) {
        throw Electrode_Error("Electrode_CSTR_LiCoO2Cathode::electrode_model_create()", "Relative Extent Rxn outside of bounds");
    }
    xRegion_init_ = 0;
    xRegion_init_init_ = 0;
    xRegion_final_ = 0;
    xRegion_final_final_ = 0;

    setCapacityCoeff_LiCoO2();

    capacityInitialZeroDod_  = capacity();

    return 0;
}
//======================================================================================================================
// Calculate the relative extent of reaction from the current state of the object
/*
 *  Calculate the relative extent of reaction from the final state, spmoles_final.
 *  This is a virtual function because there is no way to do this except by knowing about
 *  the system.
 *  The relative extent of reaction is a dimensionless number that varies. It doesn't
 *  always vary between 0 and 1. Sometimes there are Li's that can be reacted or sites
 *  that can't be filled with Li....
 *
 *  @return returns the relative extent of reaction (dimensionless).
 */
double Electrode_CSTR_LiCoO2Cathode::calcRelativeExtentRxn_final() const
{
    double li_mf = moleFraction(ig_SolidLi_);
    return (li_mf);
}
//==================================================================================================================================
// Set the final state of the electrode using the relExtentRxn
/*
 *  This sets the state of the system, i.e., spmoles_final_[] for the solid phase
 *  components of the electrode using a single number.
 *
 *  It is virtual because there is no way to do this except by knowing about the system
 *
 *  It must be the case that  calcRelativeExtentRxn_final() and setStateFinal_relativeExtentRxn()
 *  are inverses of one another. Note, this means that if the state of the system has more than one rank,
 *  then the other ranks are unperturbed by the round trip.
 *
 *  @param relExtentRxn  input of the relative extent of reaction
 */
void Electrode_CSTR_LiCoO2Cathode::setState_relativeExtentRxn(double relExtentRxn)
{
    double li_mf =  relExtentRxn;
    double sum =  spMoles_final_[ig_SolidV_] +  spMoles_final_[ig_SolidLi_];
    spMoles_final_[ig_SolidLi_] =         li_mf * sum;
    spMoles_final_[ig_SolidV_]  = (1.0 - li_mf) * sum;
    /*
     *  Now Update all of the other numbers from spMoles_final_[]
     */
    updateState_Phase(ip_LiCoO2_);
}
//==================================================================================================================================
} // End of ZZCantera namespace
//----------------------------------------------------------------------------------------------------------------------------------
