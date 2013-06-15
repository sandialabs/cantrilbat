/*
 * Electrode_FD_Jacobian.h
 *
 *  Created on: Jun 10, 2013
 *      Author: vebruni
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef ELECTRODE_FD_JACOBIAN_H_
#define ELECTRODE_FD_JACOBIAN_H_

#include <list>
#include <map>
#include "Electrode_Jacobian.h"

namespace Cantera {

// Compute electrode Jacobian by finite differences
class Electrode_FD_Jacobian : public Electrode_Jacobian {
public:
  Electrode_FD_Jacobian(Electrode* elect);

  virtual ~Electrode_FD_Jacobian();

  virtual void compute_jacobian(const std::vector<double> & centerpoint, const double dt);

  virtual void add_entry_to_compute(DOF_SOURCE_PAIR entry);
  virtual void remove_entry_to_compute(DOF_SOURCE_PAIR entry);

protected:
  std::list<DOFS> dofs_to_fd;
  std::map<DOFS, int> num_sources_using_dof;

private:
  Electrode_FD_Jacobian(const Electrode_FD_Jacobian& right);
  Electrode_FD_Jacobian& operator=(const Electrode_FD_Jacobian& right);

  void run_electrode_integration(const std::vector<double> & dof_values, double dt);
  void set_jacobian_entry(const DOF_SOURCE_PAIR & entry, double source_values[2], double delta);
};

} // namespace Cantera

#endif /* ELECTRODE_FD_JACOBIAN_H_ */
