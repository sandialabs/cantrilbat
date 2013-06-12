/*
 * $Id: Electrode_Jacobian.cpp 496 2013-01-07 21:15:37Z hkmoffa $
 */

#include "Electrode_Jacobian.h"

namespace Cantera {

//====================================================================================================================
Electrode_Jacobian::Electrode_Jacobian(Electrode* elect) :
                electrode(elect)
{
}
//====================================================================================================================
Electrode_Jacobian::~Electrode_Jacobian()
{
}
//====================================================================================================================
void Electrode_Jacobian::add_entries_to_compute(const std::vector<DOF_SOURCE_PAIR> &entries)
{
  std::vector<DOF_SOURCE_PAIR>::const_iterator dof_sources_it = entries.begin();
  std::vector<DOF_SOURCE_PAIR>::const_iterator dof_source_end = entries.end();
  for( ; dof_sources_it != dof_source_end; ++dof_sources_it )
  {
    add_entry_to_compute(*dof_sources_it);
  }
}
//====================================================================================================================
void Electrode_Jacobian::add_entry_to_compute(DOF_SOURCE_PAIR entry)
{
  if( jacobian.find(entry) == jacobian.end() )
  {
    jacobian[entry] = 0.0;
  }
}
//====================================================================================================================
void Electrode_Jacobian::remove_entry_to_compute(DOF_SOURCE_PAIR entry)
{
  std::map< DOF_SOURCE_PAIR, double >::iterator entry_pos = jacobian.find(entry);
  if( entry_pos != jacobian.end() )
  {
    jacobian.erase(entry_pos);
  }
}
//====================================================================================================================
}// End of namespace Cantera

