/*
 * $Id: Electrode_Jacobian.cpp 496 2013-01-07 21:15:37Z hkmoffa $
 */

#include "Electrode_Jacobian.h"

#include <map>
#include <vector>

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif 
{
//====================================================================================================================
Electrode_Jacobian::Electrode_Jacobian(Electrode* elect) :
    electrode(elect),
    electrolytePhaseSpeciesStart(0),
    jac_dt(0.0),
    jac_t_init_init(0.0),
    tp_solnPhase(0),
    jac_numSubs_Max(0),
    jac_numSubs_Min(0),
    jac_energySource(0.0),
    jac_electrolytePhaseSource(0.0),
    jac_electronSource(0.0)
{
    size_t solnPhaseIndex = electrode->solnPhaseIndex();
    electrolytePhaseSpeciesStart = electrode->globalSpeciesIndex(electrode->solnPhaseIndex());
    tp_solnPhase = & electrode->phase(solnPhaseIndex);
    jac_lyteSpeciesSource.resize(tp_solnPhase->nSpecies(), 0.0);
}
//====================================================================================================================
Electrode_Jacobian::Electrode_Jacobian(const Electrode_Jacobian& right) :
    electrode(right.electrode),
    jac_dt(0.0),
    jac_t_init_init(0.0)
{
    operator=(right);
}
//====================================================================================================================
Electrode_Jacobian& Electrode_Jacobian::operator=(const Electrode_Jacobian& right)
{
    if (this == &right) {
        return *this;
    }

    electrode = right.electrode;
    electrolytePhaseSpeciesStart = right.electrolytePhaseSpeciesStart; 
    jacobian = right.jacobian;
    jac_centerpoint = right.jac_centerpoint;
    jac_dt = right.jac_dt;
    jac_t_init_init = right.jac_t_init_init;
    tp_solnPhase = right.tp_solnPhase;
    jac_numSubs_Max = right.jac_numSubs_Max;
    jac_numSubs_Min = right.jac_numSubs_Min;

    jac_energySource = right.jac_energySource;
    jac_electrolytePhaseSource = right.jac_electrolytePhaseSource;
    jac_electronSource = right.jac_electronSource;
    jac_lyteSpeciesSource = right.jac_lyteSpeciesSource;

    return *this;
}
//====================================================================================================================
Electrode_Jacobian::~Electrode_Jacobian()
{
}
//====================================================================================================================
std::string Electrode_Jacobian::dofsString(enum DOFS dd) const
{
    std::string ss;
    size_t nsp = tp_solnPhase->nSpecies();
    if (dd == SOLID_VOLTAGE) {
	ss = "SOLID_VOLTAGE";
    } else if (dd == LIQUID_VOLTAGE) {
	ss = "LYTE_VOLTAGE";
    } else if (dd == TEMPERATURE) {
	ss = "TEMPERATURE";
    } else if (dd == PRESSURE) {
	ss = "PRESSURE";
    } else if (dd >= SPECIES && dd < SPECIES + nsp) {
	size_t k = (size_t) (dd - SPECIES);
	ss = "C_" +  tp_solnPhase->speciesName(k);
    } else {
	throw CanteraError("Electrode_Jacobian::dofsString", "unknown DOF value" + int(dd)); 
    }
    return ss; 
}
//====================================================================================================================
void Electrode_Jacobian::add_entries_to_compute(const std::vector<DOF_SOURCE_PAIR> &entries)
{
    // Here we process a vector of entries adding each entry one at a time.
    std::vector<DOF_SOURCE_PAIR>::const_iterator dof_sources_it = entries.begin();
    std::vector<DOF_SOURCE_PAIR>::const_iterator dof_source_end = entries.end();
    for( ; dof_sources_it != dof_source_end; ++dof_sources_it) {
      add_entry_to_compute(*dof_sources_it);
    }
}
//====================================================================================================================
void Electrode_Jacobian::add_entry_to_compute(DOF_SOURCE_PAIR entry)
{
    // Here we check for an existing entry. If we don't find one, we add a zero entry in the jacobian map
    if (jacobian.find(entry) == jacobian.end()) {
      jacobian[entry] = 0.0;
    }
}
//====================================================================================================================
void Electrode_Jacobian::remove_entry_to_compute(DOF_SOURCE_PAIR entry)
{
    // Here we eliminate a jacobian entry.
    std::map< DOF_SOURCE_PAIR, double >::iterator entry_pos = jacobian.find(entry);
    if (entry_pos != jacobian.end()) {
      jacobian.erase(entry_pos);
    }
}
//====================================================================================================================
} // End of namespace
//----------------------------------------------------------------------------------------------------------------------------------
