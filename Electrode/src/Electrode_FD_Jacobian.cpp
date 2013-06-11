/*
 * Electrode_FD_Jacobian.cpp
 *
 *  Created on: Jun 10, 2013
 *      Author: vebruni
 */

#include "Electrode_FD_Jacobian.h"

namespace Cantera {

//====================================================================================================================
Electrode_FD_Jacobian::Electrode_FD_Jacobian(Electrode * elect, const std::vector<DOF_SOURCE_PAIR> & entries_to_compute)
  : Electrode_Jacobian(elect, entries_to_compute)
{
}
//====================================================================================================================
Electrode_FD_Jacobian::~Electrode_FD_Jacobian()
{
}
//====================================================================================================================
void Electrode_FD_Jacobian::compute_jacobian(const std::vector<double> & centerpoint, const double dt)
{
  std::list<DOFS>::iterator dof_it = dofs_to_fd.begin();
  for( ; dof_it != dofs_to_fd.end(); ++dof_it )
  {
    for( int negative=0; negative<2; ++negative)
    {
      // Perturb by -1^negative * base_delta and compute
    }
  }
}
//====================================================================================================================
void Electrode_FD_Jacobian::add_entry_to_compute(DOF_SOURCE_PAIR entry)
{
  if( jacobian.find(entry) == jacobian.end() )
  {
    jacobian[entry] = 0.0;
    if( std::find(dofs_to_fd.begin(), dofs_to_fd.end(), entry.first) == dofs_to_fd.end() )
    {
      dofs_to_fd.push_back(entry.first);
      num_sources_using_dof[entry.first] = 1;
    }
    else
    {
      ++num_sources_using_dof[entry.first];
    }
  }
}
//====================================================================================================================
void Electrode_FD_Jacobian::remove_entry_to_compute(DOF_SOURCE_PAIR entry)
{
  std::map< DOF_SOURCE_PAIR, double >::iterator entry_pos = jacobian.find(entry);
  if( entry_pos != jacobian.end() )
  {
    jacobian.erase(entry_pos);
    --num_sources_using_dof[entry.first];
    if( num_sources_using_dof[entry.first] == 0 )
    {
      dofs_to_fd.erase( std::find(dofs_to_fd.begin(), dofs_to_fd.end(), entry.first) );
    }
  }
}
//====================================================================================================================
void Electrode_FD_Jacobian::run_electrode_integration(const std::vector<double> & dof_values, const double dt)
{
  electrode->setState_TP(dof_values[TEMPERATURE], dof_values[PRESSURE]);
  electrode->setVoltages(dof_values[SOLID_VOLTAGE], dof_values[LIQUID_VOLTAGE]);
  electrode->setFinalStateFromInit();
  electrode->setElectrolyteMoleNumbers(&dof_values[SPECIES], true);
  electrode->integrate(dt);
}
//====================================================================================================================

} // namespace Cantera
