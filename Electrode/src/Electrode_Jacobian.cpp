/*
 * $Id: Electrode_Jacobian.cpp 496 2013-01-07 21:15:37Z hkmoffa $
 */

#include "Electrode_Jacobian.h"

namespace Cantera {

//====================================================================================================================
Electrode_Jacobian::Electrode_Jacobian(Electrode* elect, const std::vector<DOF_SOURCE_PAIR> & entries_to_compute) :
                printLvl_(0),
                electrode(elect)
{
  std::map<DOF_SOURCE_PAIR, bool>::iterator dof_sources_it = entries_to_compute.begin();
  std::map<DOF_SOURCE_PAIR, bool>::const_iterator dof_source_end = entries_to_compute.end();
  for( ; dof_sources_it != dof_source_end; ++dof_sources_it )
  {
    jacobian[*dof_sources_it] = 0.0;
  }
}
//====================================================================================================================
Electrode_Jacobian::Electrode_Jacobian(const Electrode_Jacobian& right) :
                printLvl_(0),
                electrode(right.electrode),
                jacobian(right.jacobian)
{
    operator=(right);
}
//====================================================================================================================
Electrode_Jacobian::~Electrode_Jacobian()
{
}
//======================================================================================================================
Electrode_Jacobian& Electrode_Jacobian::operator=(const Electrode_Jacobian& right)
{
    if (this == &right) {
        return *this;
    }
    electrode = right.electrode;

    printLvl_ = right.printLvl_;

    jacobian = right.jacobian;

    return *this;
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
//======================================================================================================================

