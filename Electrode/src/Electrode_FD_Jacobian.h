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


namespace Cantera {

// Compute electrode Jacobian by finite differences
class Electrode_FD_Jacobian : public Electrode_Jacobian {
  public:
    void compute_jacobian();

  protected:
};

} // namespace Cantera

#endif /* ELECTRODE_FD_JACOBIAN_H_ */
