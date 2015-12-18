/*
 * $Id: Electrode_Jacobian.cpp 496 2013-01-07 21:15:37Z hkmoffa $
 */

#include "Electrode_Jacobian.h"

#include <map>

namespace Cantera {

//====================================================================================================================
Electrode_Jacobian::Electrode_Jacobian(Electrode* elect) :
    electrode(elect)
{
    electrolytePhaseSpeciesStart = electrode->getGlobalSpeciesIndex(electrode->solnPhaseIndex());
}
//====================================================================================================================
Electrode_Jacobian::Electrode_Jacobian(const Electrode_Jacobian& right) :
    electrode(right.electrode)
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

    return *this;
}
//====================================================================================================================
Electrode_Jacobian::~Electrode_Jacobian()
{
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
}// End of namespace Cantera

