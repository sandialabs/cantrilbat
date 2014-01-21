/*
 * Electrode_FD_Jacobian.cpp
 *
 *  Created on: Jun 10, 2013
 *      Author: vebruni
 */

#include "Electrode_FD_Jacobian.h"

namespace Cantera {

//====================================================================================================================
Electrode_FD_Jacobian::Electrode_FD_Jacobian(Electrode * elect)
  : Electrode_Jacobian(elect)
{
}
//====================================================================================================================
Electrode_FD_Jacobian::~Electrode_FD_Jacobian()
{
}
//====================================================================================================================
void Electrode_FD_Jacobian::compute_jacobian(const std::vector<double> & centerpoint, const double dt)
{
  // Temporary storage used to do the centered difference calculation
  std::vector<double> speciesSources[2];
  speciesSources[0].resize(electrode->nSpecies());
  speciesSources[1].resize(electrode->nSpecies());
  double energySource[2];
  double electrolytePhaseSource[2];
  double individualSpeciesSource[2];

  int electrolytePhaseIndex = electrode->solnPhaseIndex();
  int electronIndex = electrode->kSpecElectron();

  std::vector<double> perturbed_point = centerpoint;

  const double base_delta = 1.e-5;

  // Iterate over the list of dofs that need to be finite differenced with respect to
  // For each dof store the source term values with the dof perturbed +- its centerpoint value.
  // These are then used to set the Jacobian entries using a centered difference formula
  std::list<DOFS>::iterator dof_it = dofs_to_fd.begin();
  for( ; dof_it != dofs_to_fd.end(); ++dof_it )
  {
    // Determine how far to perturb the dof using base_delta as a relative perturbation magnitude
    double dof_delta = base_delta;
    dof_delta *= (std::abs(perturbed_point[*dof_it]) > 0.) ? std::abs(perturbed_point[*dof_it]) : 1.0;

    for( int negative=0; negative<2; ++negative)
    {
      // Perturb by -1^negative * base_delta and compute
      perturbed_point[*dof_it] += dof_delta * std::pow(-1.0, negative);
      run_electrode_integration(perturbed_point, dt);

      // Store off the source term results in temporary storage
      energySource[negative] = electrode->getSourceTerm(Cantera::ENTHALPY_SOURCE);
      electrolytePhaseSource[negative] = electrode->getSourceTerm(Cantera::ELECTROLYTE_PHASE_SOURCE);
      electrode->integratedSourceTerm(&speciesSources[negative][0]);

      perturbed_point[*dof_it] = centerpoint[*dof_it];
    }

    // Use set_jacobian_entry to store the various sensitivities based on a centered difference.
    set_jacobian_entry( DOF_SOURCE_PAIR(*dof_it, ENTHALPY_SOURCE), energySource, dof_delta);
    set_jacobian_entry( DOF_SOURCE_PAIR(*dof_it, ELECTROLYTE_PHASE_SOURCE), electrolytePhaseSource, dof_delta);
    individualSpeciesSource[0] = speciesSources[0][electronIndex];
    individualSpeciesSource[1] = speciesSources[1][electronIndex];
    set_jacobian_entry( DOF_SOURCE_PAIR(*dof_it, CURRENT_SOURCE), individualSpeciesSource, dof_delta);
    for(int sp=0; sp < electrode->numSolnPhaseSpecies(); ++sp)
    {
      int idx = sp + electrolytePhaseSpeciesStart;
      individualSpeciesSource[0] = speciesSources[0][idx];
      individualSpeciesSource[1] = speciesSources[1][idx];
      set_jacobian_entry( DOF_SOURCE_PAIR(*dof_it, (SOURCES)(SPECIES_SOURCE + sp)), individualSpeciesSource, dof_delta);
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
void Electrode_FD_Jacobian::run_electrode_integration(const std::vector<double> & dof_values, double dt)
{
  electrode->revertToInitialTime(true);
  electrode->setState_TP(dof_values[TEMPERATURE], dof_values[PRESSURE]);
  electrode->setVoltages(dof_values[SOLID_VOLTAGE], dof_values[LIQUID_VOLTAGE]);
  electrode->setFinalStateFromInit();
  electrode->setElectrolyteMoleNumbers(&dof_values[SPECIES], true);
  electrode->integrate(dt);
}
//====================================================================================================================
void Electrode_FD_Jacobian::set_jacobian_entry(const DOF_SOURCE_PAIR & entry, double source_values[2], double delta)
{
  if(jacobian.find(entry) != jacobian.end())
  {
    jacobian[entry] = (source_values[0] - source_values[1]) / (2*delta);
  }
}
//====================================================================================================================

} // namespace Cantera
