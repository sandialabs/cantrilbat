/*
 * $Id: Electrode_Jacobian.cpp 496 2013-01-07 21:15:37Z hkmoffa $
 */

#include "Electrode_Jacobian.h"

namespace Cantera {

//====================================================================================================================
Electrode_Jacobian::Electrode_Jacobian(Electrode* elect, const std::map<DOF_SOURCE_PAIR, bool> & dof_source_pairs) :
                printLvl_(0),
                electrode(elect)
{
  std::map<DOF_SOURCE_PAIR, bool>::const_iterator dof_source_end = dof_source_pairs.end();
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
}// End of namespace Cantera
//======================================================================================================================

