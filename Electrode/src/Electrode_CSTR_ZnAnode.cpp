/*
 * $Id: Electrode_CSTR_ZnAnode.cpp 571 2013-03-26 16:44:21Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"

#include "Electrode_CSTR_ZnAnode.h"

namespace Cantera
{
//======================================================================================================================
Electrode_CSTR_ZnAnode::Electrode_CSTR_ZnAnode() :
    Electrode_CSTR()
{
}
//======================================================================================================================
Electrode_CSTR_ZnAnode::Electrode_CSTR_ZnAnode(const Electrode_CSTR_ZnAnode& right) :
    Electrode_CSTR()
{
    *this = operator=(right);
}
//======================================================================================================================
// Assignment operator
/*
 *  @param right object to be copied
 */
Electrode_CSTR_ZnAnode& Electrode_CSTR_ZnAnode::operator=(const Electrode_CSTR_ZnAnode& right)
{
    // Check for self assignment.
    if (this == &right) {
        return *this;
    }

    Electrode_CSTR::operator=(right);

    return *this;
}
//======================================================================================================================
Electrode_CSTR_ZnAnode::~Electrode_CSTR_ZnAnode()
{
}
//======================================================================================================================
Electrode_Types_Enum  Electrode_CSTR_ZnAnode::electrodeType() const
{
    return CSTR_ZN_ANODE_ET;
}
//======================================================================================================================
int
Electrode_CSTR_ZnAnode::electrode_model_create(ELECTRODE_KEY_INPUT* ei)
{
    int flag = Electrode_CSTR::electrode_model_create(ei);
    if (flag != 0) {
        return flag;
    }

    init_sizes();

    zn_phase_index = globalPhaseIndex("Zn(S)");
    if (zn_phase_index < 0) {
        throw CanteraError("Electrode_CSTR_ZnAnode::electrode_model_create()",
                           "Phase Zn(S) not found.");
    }

    zn_sp_index = globalSpeciesIndex("Zn(S)", "Zn(S)");
    if (zn_sp_index < 0) {
        throw CanteraError("Electrode_CSTR_ZnAnode::electrode_model_create()",
                           "Species Zn(S) in phase Zn(S) not found.");
    }

    RelativeExtentRxn_RegionBoundaries_.resize(2);

    RelativeExtentRxn_RegionBoundaries_[1] = 1.0;
    RelativeExtentRxn_RegionBoundaries_[0] = 0.0;

    RelativeExtentRxn_init_ = 0.0;
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
double Electrode_CSTR_ZnAnode::calcRelativeExtentRxn_final() const
{
    // RelativeExtentRxn_NormalizationFactor_ contains the initial number of moles of solid.
    // For the Zn electrode that is the initial moles of Zn so we can directly use it to
    // calculate the extent of reaction.
    return (1.0 - spMoles_final_[zn_sp_index] / RelativeExtentRxn_NormalizationFactor_);
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
void Electrode_CSTR_ZnAnode::setStateFinal_fromRelativeExtentRxn(double relExtentRxn)
{
    spMoles_final_[zn_sp_index] = (1.0 - relExtentRxn) * RelativeExtentRxn_NormalizationFactor_;

    updatePhaseNumbers(zn_phase_index);
}

} // End of namespace Cantera
