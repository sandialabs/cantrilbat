/*
 * $Id: Electrode_CSTR_MCMBAnode.cpp 571 2013-03-26 16:44:21Z hkmoffa $
 */

#include "Electrode_CSTR_MCMBAnode.h"

namespace Cantera
{
//======================================================================================================================
/*
 *  ELECTRODE_INPUT: constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.
 */
Electrode_CSTR_MCMBAnode::Electrode_CSTR_MCMBAnode() :
    Electrode_CSTR(),
    ig_SolidLi_(-1),
    ig_SolidV_(-1),
    ip_MCMB_(-1)
{

}
//======================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
Electrode_CSTR_MCMBAnode::Electrode_CSTR_MCMBAnode(const Electrode_CSTR_MCMBAnode& right) :
    Electrode_CSTR(),
    ig_SolidLi_(-1),
    ig_SolidV_(-1),
    ip_MCMB_(-1)
{
    /*
     * Call the assignment operator.
     */
    *this = operator=(right);
}
//======================================================================================================================
// Assignment operator
/*
 *  @param right object to be copied
 */
Electrode_CSTR_MCMBAnode& Electrode_CSTR_MCMBAnode::operator=(const Electrode_CSTR_MCMBAnode& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    Electrode_CSTR::operator=(right);

    ig_SolidLi_ = right.ig_SolidLi_;
    ig_SolidV_  = right.ig_SolidV_;
    ip_MCMB_    = right.ip_MCMB_;

    /*
     * Return the reference to the current object
     */
    return *this;
}
//======================================================================================================================
/*
 *
 *  ELECTRODE_INPUT:destructor
 *
 * We need to manually free all of the arrays.
 */
Electrode_CSTR_MCMBAnode::~Electrode_CSTR_MCMBAnode()
{

}
//======================================================================================================================
//     Return the type of electrode
/*
 *  Returns the enum type of the electrode. This is used in the factory routine.
 *
 *  @return Returns an enum type, called   Electrode_Types_Enum
 */
Electrode_Types_Enum  Electrode_CSTR_MCMBAnode::electrodeType() const
{
    return CSTR_MCMB_ANODE_ET;
}
//======================================================================================================================
//  Setup the electrode
/*
 * @param ei    ELECTRODE_KEY_INPUT pointer object
 */
int
Electrode_CSTR_MCMBAnode::electrode_model_create(ELECTRODE_KEY_INPUT* ei)
{


    int flag = Electrode_CSTR::electrode_model_create(ei);
    if (flag != 0) {
        return flag;
    }

    init_sizes();


    ig_SolidLi_ = globalSpeciesIndex("Li_C6-bulk");
    if (ig_SolidLi_ < 0) {
        throw CanteraError("Electrode_CSTR_MCMBAnode::electrode_model_create()",
                           "I tried to find the species Li_C6-bulk but failed. May need to generalize the code now");
    }

    ig_SolidV_  = globalSpeciesIndex("V_C6-bulk");
    if (ig_SolidV_ < 0) {
        throw CanteraError("Electrode_CSTR_MCMBAnode::electrode_model_create()",
                           "I tried to find the species V_C6-bulk but failed. May need to generalize the code now");
    }

    ip_MCMB_ = globalPhaseIndex("MCMB_Interstitials_anode");
    if (ip_MCMB_ < 0) {
        throw CanteraError("Electrode_CSTR_MCMBAnode::electrode_model_create()",
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
        throw CanteraError("Electrode_CSTR::electrode_model_create()", "Relative Extent Rxn outside of bounds");
    }
    if (RelativeExtentRxn_final_ >  RelativeExtentRxn_RegionBoundaries_[1]) {
        throw CanteraError("Electrode_CSTR::electrode_model_create()", "Relative Extent Rxn outside of bounds");
    }
    xRegion_init_ = 0;
    xRegion_init_init_ = 0;
    xRegion_final_ = 0;
    xRegion_final_final_ = 0;



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
double Electrode_CSTR_MCMBAnode::calcRelativeExtentRxn_final() const
{
    double li_mf = moleFraction(ig_SolidLi_);
    return (1.0 - li_mf);
}
//======================================================================================================================
// Set the final state of the electrode using the relExtentRxn
/*
 *  This sets the state of the system, i.e., spmoles_final_[] for the solid phase
 *  components of the electrode using a single number.
 *
 *  It is virtual because there is no way to do this except by knowing about the system
 *
 *  It must be the case that  calcRelativeExtentRxn_final() and etStateFinal_fromRelativeExtentRxn()
 *  are inverses of one another. Note, this means that if the state of the system has more than one rank,
 *  then the other ranks are unperturbed by the round trip.
 *
 *  @param relExtentRxn  input of the relative extent of reaction.
 *                       Dimensionless -> equal to lithium vacancy mole fraction.
 */
void Electrode_CSTR_MCMBAnode::setState_relativeExtentRxn(double relExtentRxn)
{
    double li_mf =  1.0 - relExtentRxn;
    double sum =  spMoles_final_[ig_SolidV_] +  spMoles_final_[ig_SolidLi_];
    spMoles_final_[ig_SolidLi_] =         li_mf * sum;
    spMoles_final_[ig_SolidV_]  = (1.0 - li_mf) * sum;

    if (ip_MCMB_ < 0) {
        throw CanteraError("ddddd",
                           "ERROR");
    }
    /*
     *  Now Update all of the other numbers from spMoles_final_[]
     */
    updateState_Phase(ip_MCMB_);
}
//====================================================================================================================
} // End of namespace Cantera
//======================================================================================================================
