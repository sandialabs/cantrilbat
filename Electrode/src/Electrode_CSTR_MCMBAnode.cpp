/**
 *  @file Electrode_CSTR_MCMBAnode.cpp
 *     Definitions for the Electrode_CSTR class, used to model 
 *     Electrode processes in particles with no transport limits hard-coded for the MCMB Anode
 *     (see \ref electrode_mgr and class \link Zuzax::Electrode_CSTR_MCMBAnode Electrode_CSTR_MCMBAnode \endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "Electrode_CSTR_MCMBAnode.h"

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
Electrode_CSTR_MCMBAnode::Electrode_CSTR_MCMBAnode() :
    Electrode_CSTR(),
    ig_SolidLi_(npos),
    ig_SolidV_(npos),
    ip_MCMB_(npos)
{
}
//==================================================================================================================================
Electrode_CSTR_MCMBAnode::Electrode_CSTR_MCMBAnode(const Electrode_CSTR_MCMBAnode& right) :
    Electrode_CSTR_MCMBAnode()
{
    operator=(right);
}
//==================================================================================================================================
Electrode_CSTR_MCMBAnode& Electrode_CSTR_MCMBAnode::operator=(const Electrode_CSTR_MCMBAnode& right)
{
    if (this == &right) {
        return *this;
    }
    Electrode_CSTR::operator=(right);
    ig_SolidLi_ = right.ig_SolidLi_;
    ig_SolidV_  = right.ig_SolidV_;
    ip_MCMB_    = right.ip_MCMB_;
    return *this;
}
//==================================================================================================================================
Electrode* Electrode_CSTR_MCMBAnode::duplMyselfAsElectrode() const
{
    Electrode_CSTR_MCMBAnode* dd = new Electrode_CSTR_MCMBAnode(*this);
    return dd;
}
//==================================================================================================================================
Electrode_CSTR_MCMBAnode::~Electrode_CSTR_MCMBAnode()
{
}
//==================================================================================================================================
Electrode_Types_Enum Electrode_CSTR_MCMBAnode::electrodeType() const
{
    return CSTR_MCMB_ANODE_ET;
}
//==================================================================================================================================
int Electrode_CSTR_MCMBAnode::electrode_model_create(ELECTRODE_KEY_INPUT* ei)
{
    int flag = Electrode_CSTR::electrode_model_create(ei);
    if (flag != 0) {
        return flag;
    }
    init_sizes();

    ig_SolidLi_ = globalSpeciesIndex("Li_C6-bulk");
    if (ig_SolidLi_ == npos) {
        throw Electrode_Error("Electrode_CSTR_MCMBAnode::electrode_model_create()",
                           "I tried to find the species Li_C6-bulk but failed. May need to generalize the code now");
    }

    ig_SolidV_  = globalSpeciesIndex("V_C6-bulk");
    if (ig_SolidV_ == npos) {
        throw Electrode_Error("Electrode_CSTR_MCMBAnode::electrode_model_create()",
                           "I tried to find the species V_C6-bulk but failed. May need to generalize the code now");
    }

    ip_MCMB_ = globalPhaseIndex("MCMB_Interstitials_anode");
    if (ip_MCMB_ == npos) {
        throw Electrode_Error("Electrode_CSTR_MCMBAnode::electrode_model_create()",
                           "I tried to find the phase MCMB_Interstitials_anode  but failed. May need to generalize the code now");
    }

    RelativeExtentRxn_RegionBoundaries_.resize(2);

    RelativeExtentRxn_RegionBoundaries_[1] = 0.94;
    RelativeExtentRxn_RegionBoundaries_[0] = 0.10;
    double li_mf = moleFraction(ig_SolidLi_);

    RelativeExtentRxn_init_ = 1.0 -  li_mf;
    RelativeExtentRxn_final_ = RelativeExtentRxn_init_;
    RelativeExtentRxn_final_final_ = RelativeExtentRxn_init_;
    RelativeExtentRxn_init_init_ = RelativeExtentRxn_init_;

    if (RelativeExtentRxn_final_ <  RelativeExtentRxn_RegionBoundaries_[0]) {
        throw Electrode_Error("Electrode_CSTR::electrode_model_create()", "Relative Extent Rxn outside of bounds");
    }
    if (RelativeExtentRxn_final_ >  RelativeExtentRxn_RegionBoundaries_[1]) {
        throw Electrode_Error("Electrode_CSTR::electrode_model_create()", "Relative Extent Rxn outside of bounds");
    }
    xRegion_init_ = 0;
    xRegion_init_init_ = 0;
    xRegion_final_ = 0;
    xRegion_final_final_ = 0;
    return 0;
}
//==================================================================================================================================
double Electrode_CSTR_MCMBAnode::calcRelativeExtentRxn_final() const
{
    double li_mf = moleFraction(ig_SolidLi_);
    return (1.0 - li_mf);
}
//==================================================================================================================================
void Electrode_CSTR_MCMBAnode::setState_relativeExtentRxn(double relExtentRxn)
{
    double li_mf =  1.0 - relExtentRxn;
    double sum =  spMoles_final_[ig_SolidV_] +  spMoles_final_[ig_SolidLi_];
    spMoles_final_[ig_SolidLi_] =         li_mf * sum;
    spMoles_final_[ig_SolidV_]  = (1.0 - li_mf) * sum;

    if (ip_MCMB_ == npos) {
        throw Electrode_Error("ddddd", "ERROR");
    }
    updateState_Phase(ip_MCMB_);
}
//==================================================================================================================================
} 
//----------------------------------------------------------------------------------------------------------------------------------
