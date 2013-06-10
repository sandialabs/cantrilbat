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

#include <map>

namespace Cantera {

  enum SOURCES
  {
    CURRENT,
    SPECIES,
    ENTHALPY
  };
  enum DOFS
  {
    SOLID_VOLTAGE,
    LIQUID_VOLTAGE,
    SPECIES,
    TEMPERATURE,
    PRESSURE
  };

class Electrode_Jacobian {

public:
    typedef std::pair<DOFS, SOURCES> DOF_SOURCE_PAIR ;

    Electrode_Jacobian(Electrode* elect, const std::map<DOF_SOURCE_PAIR, bool> & dofs);

    virtual ~Electrode_Jacobian();

    Electrode_Jacobian(const Electrode_Jacobian& right);

    Electrode_Jacobian& operator=(const Electrode_Jacobian& right);

    //! Print level for input to vcs routines within this object
    /*!
     *    This is a public member so that it can be manipulated
     */
    int printLvl_;

    void compute_jacobian() = 0;

    double get_jacobian_value(DOF_SOURCE_PAIR dof_source_pair) { return jacobian[dof_source_pair].second; }

protected:

    Electrode* electrode;

    // Store the desired Jacobian contributions as a map from [dof, source] -> [ calculate?, result]
    std::map<DOF_SOURCE_PAIR, std::pair<bool, double> > jacobian;
};

}
#endif
/*****************************************************************************/

