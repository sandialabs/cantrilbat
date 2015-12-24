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

#include "Electrode_Jacobian.h"

#include <list>
#include <map>

namespace Cantera {

//! Implementation of the Electrode_Jacobian interface that uses numerical finite differencing to calculate source term sensitivity values.
class Electrode_FD_Jacobian : public Electrode_Jacobian {
public:
    Electrode_FD_Jacobian(Electrode* elect, double baseDelta);

    Electrode_FD_Jacobian(const Electrode_FD_Jacobian& right);

    Electrode_FD_Jacobian& operator=(const Electrode_FD_Jacobian& right);

  virtual ~Electrode_FD_Jacobian();

  virtual void compute_jacobian(const std::vector<double> & centerpoint, const double dt);

  virtual void add_entry_to_compute(DOF_SOURCE_PAIR entry);
  virtual void remove_entry_to_compute(DOF_SOURCE_PAIR entry);

  //! @copydoc Electrode_Jacobian::print_jacobian
  virtual void print_jacobian(int indentSpaces = 0) const;

protected:
  // These variables are used to determine which dofs we need to finite
  // difference with respect to.
  std::list<DOFS> dofs_to_fd;
  std::map<DOFS, int> num_sources_using_dof;

  //! Relative delta to use to create the jacobian
  double base_delta;

private:

  // This helper function handles running a single electrode->integrate call at
  // a given set of dof values.
  /*!
   *  @return              Returns the number of subintegration steps
   */
  int run_electrode_integration(const std::vector<double> & dof_values, double dt);

  // This helper function sets a specific jacobian entry using a centered difference formula
  // based on the specified source values and delta in the dof.
  void set_jacobian_entry(const DOF_SOURCE_PAIR & entry, double source_values[2], double delta);
};

} // namespace Cantera

#endif /* ELECTRODE_FD_JACOBIAN_H_ */
