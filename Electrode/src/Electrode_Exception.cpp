/**
 * @file Electrode_Exception.cpp
 *
 */
#include "Electrode_Exception.h"
#include "Electrode.h"
#include "Electrode_Factory.h"
#include <cmath>
//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
Electrode_Error::Electrode_Error(const std::string &proc, const std::string &msg) :
    ZZCantera::CanteraError("Electrode_Error: " + proc, msg)
{
}
//==================================================================================================================================
Electrode_Error::~Electrode_Error() throw()
{
}
//==================================================================================================================================
Electrode_Error::Electrode_Error() :
    ZZCantera::CanteraError()
{
}

}
//----------------------------------------------------------------------------------------------------------------------------------
namespace esmodel 
{
//==================================================================================================================================
void Electrode_Warning(const ZZCantera::Electrode& e,  const std::string &procedure, const std::string &msg)
{
    int eDom = e.electrodeDomainNumber_;
    int eCell = e.electrodeCellNumber_;
    int printLvl_ = e.printLvl_;

    if (printLvl_) {
    ZZCantera::Electrode_Types_Enum etype = e.electrodeType();
    std::string estring = Electrode_Types_Enum_to_string(etype);
    std::string pmsg = "Electrode_Warning: dom: " + ZZCantera::int2str(eDom) + " Cell: " + ZZCantera::int2str(eCell) + ": " + estring + ":" + procedure;
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
bool doubleEqualNoAtol(double a1, double a2, int digits)
{
    double atol = 0.5 * (fabs(a1) + fabs(a2));
    if (atol == 0.0) {
        return true;
    }
    atol = atol * pow(10.0, -digits) + 1.0E-300;
    return doubleEqual(a1, a2, atol, digits);
}
//==================================================================================================================================
bool doubleEqual(double a1, double a2, double atol, int digits)
{
    double diff = fabs(a1 - a2);
    if (diff < atol) {
        return true;
    }
    double rel = diff * 2.0 / (fabs(a1) + fabs(a2) + fabs(atol));
    if (digits == 6) {
        if (rel < 4.0E-6) {
            return true;
        }
    } else {
        double tol = 4.0 * pow(10.0, -digits);
        if (rel < tol) {
            return true;
        }
    }
    return false;
}
//==================================================================================================================================
bool doubleVectorEqualNoAtol(const std::vector<double>& a1, const std::vector<double>& a2, int digits)
{
    size_t j = a1.size();
    size_t j2 = a2.size();
    if (j2 != j) {
        return false; 
    } 
    for (size_t i = 0; i < j; ++i) {
       if (!doubleEqualNoAtol(a1[i], a2[i], digits)) {
           return false;
       }
    }
    return true;
}
//==================================================================================================================================
bool doubleVectorEqual(const std::vector<double>& a1, const std::vector<double>& a2, double atol, int digits)
{
    AssertTrace(atol > 0.0);
    AssertTrace(digits < 16);
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
