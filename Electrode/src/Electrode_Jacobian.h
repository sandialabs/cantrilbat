/*
 * $Id: Electrode_Jacobian.h
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_JACOBIAN_H
#define _ELECTRODE_JACOBIAN_H

#include "Electrode.h"

namespace Cantera {

  enum SOURCES
  {
    CURRENT_SOURCE,
    ELECTROLYTE_PHASE_SOURCE,
    ENTHALPY_SOURCE,
    SPECIES_SOURCE,
    MAX_SOURCE
  };
  enum DOFS
  {
    SOLID_VOLTAGE,
    LIQUID_VOLTAGE,
    TEMPERATURE,
    PRESSURE,
    SPECIES,
    MAX_DOF
  };

class Electrode_Jacobian {

public:
  typedef std::pair<DOFS, SOURCES> DOF_SOURCE_PAIR ;

  Electrode_Jacobian(Electrode* elect, const std::vector<DOF_SOURCE_PAIR> & entries_to_compute);

  virtual ~Electrode_Jacobian();

  // Compute the Jacobian at the point specified by centerpoint where
  // centerpoint[DOFS] = dof_value
  // centerpoint.size() == (MAX_DOF + n_species - 1) (i.e. a value for all possible dofs must be specified even if
  // Jacobian entries are not being computed for that dof)
  // Species mole fractions are specified starting at centerpoint[SPECIES] and should be in the order expected by
  // electrode->setElectrolyteMoleNumbers
  virtual void compute_jacobian(const std::vector<double> & centerpoint, const double dt) = 0;

  double get_jacobian_value(DOF_SOURCE_PAIR dof_source_pair) const { return jacobian[dof_source_pair]; }

  virtual void add_entry_to_compute(DOF_SOURCE_PAIR entry);
  virtual void remove_entry_to_compute(DOF_SOURCE_PAIR entry);

protected:
  Electrode* const electrode;

  // Store the desired Jacobian contributions as a map from [dof, source] -> result
  std::map< DOF_SOURCE_PAIR, double > jacobian;

private:
  // Disallow copies
  Electrode_Jacobian(const Electrode_Jacobian& right);
  Electrode_Jacobian& operator=(const Electrode_Jacobian& right);
};

}
#endif
/*****************************************************************************/

