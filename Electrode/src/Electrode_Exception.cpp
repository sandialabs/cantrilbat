/**
 * @file Electrode_Exception.cpp
 *
 */
#include "Electrode_Exception.h"
#include <cmath>
//----------------------------------------------------------------------------------------------------------------------------------
namespace Cantera
{
//==================================================================================================================================
Electrode_Error::Electrode_Error(const std::string &proc, const std::string &msg) :
    Cantera::CanteraError("Electrode_Error: " + proc, msg)
{
}
//==================================================================================================================================
Electrode_Error::~Electrode_Error() throw()
{
}
//==================================================================================================================================
Electrode_Error::Electrode_Error() :
    Cantera::CanteraError()
{
}
//==================================================================================================================================
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
//==================================================================================================================================
bool doubleVectorEqual(const std::vector<double>& a1, const std::vector<double>& a2, double atol, int digits)
{
    size_t j = a1.size();
    size_t j2 = a2.size();
    if (j2 != j) {
        return false; 
    } 
    for (size_t i = 0; i < j; ++i) {
       if (! doubleEqual(a1[i], a2[i], atol, digits)) {
           return false;
       }
    }
    return true;
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
