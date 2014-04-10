/**
 * @file Electrode_Exception.cpp
 *
 */


#include "Electrode_Exception.h"
#include <cmath>

namespace Cantera
{

//==============================================================================
Electrode_Error::Electrode_Error(const std::string &proc, const std::string &msg) :
    CanteraError("Electrode_Error: " + proc, msg)
{
}
//==============================================================================
Electrode_Error::~Electrode_Error() throw()
{
}
//==============================================================================
Electrode_Error::Electrode_Error() :
    CanteraError()
{
}
//==============================================================================

bool doubleEqual(double a1, double a2, double atol, int digits)
{
    double denom = fabs(a1) + fabs(a2) + fabs(atol);
    double diff = fabs(a1 - a2);
    double rel = diff * 2.0 / denom;
    double tol = pow(10.0, -digits);
    if (rel < tol) {
        return true;
    }
    return false;
}

}
