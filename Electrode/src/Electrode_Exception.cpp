/**
 * @file Electrode_Exception.cpp
 *
 */
#include "Electrode_Exception.h"
#include "Electrode.h"
#include "Electrode_Factory.h"
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

void Electrode_Warning(const Electrode& e,  const std::string &procedure, const std::string &msg)
{
    int eDom = e.electrodeDomainNumber_;
    int eCell = e.electrodeCellNumber_;
    int printLvl_ = e.printLvl_;

    if (printLvl_) {
    Electrode_Types_Enum etype = e.electrodeType();
    std::string estring = Electrode_Types_Enum_to_string(etype);
    std::string pmsg = "Electrode_Warning: dom: " + int2str(eDom) + " Cell: " + int2str(eCell) + ": " + estring + ":" + procedure;
    std::cerr << pmsg << " " << msg << std::endl; 
   }
}
//==================================================================================================================================
void ESModel_Warning(const std::string &procedure, const std::string &msg)
{
    int printLvl_ = 1;
    if (printLvl_) {
    std::string pmsg = "ESModel Warning: " + procedure;
    std::cerr << pmsg << " " << msg << std::endl; 
    }
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
