/**
 * @file Electrode_Exception.cpp
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "Electrode_Exception.h"
#include "Electrode.h"
#include "Electrode_Factory.h"
#include <cmath>

#ifndef MAX
//! Define a fast max operator for doubles
#define MAX(x,y) (( (x) > (y) ) ? (x) : (y))
#endif

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
Electrode_Error::Electrode_Error(const std::string &proc, const std::string &msg) :
    ZZCantera::ZuzaxError("Electrode_Error: " + proc, msg)
{
}
//==================================================================================================================================
Electrode_Error::Electrode_Error(const std::string &proc, const char* fmt, ...) :
    ZZCantera::ZuzaxError()
{
    procedure_ = "Electrode_Error: " + proc;    
    msg_.resize(1024);
    va_list args;
    va_start(args, fmt);
    char* sbuf = const_cast<char*>(msg_.data());
#ifdef _MSC_VER
    int n = _vsnprintf(sbuf, 1023, fmt, args);
#else
    int n = vsnprintf(sbuf, 1023, fmt, args);
#endif 
    if (n < 1024) {
        // if n is negative, we just go ahead and put a zero at the end of the buffer and write anyway
        va_end(args);
        sbuf[1023] = '\0';
        msg_.resize(n+1);
    } else {
        int sze = n + 1;
        msg_.resize(sze);
        va_start(args, fmt);
        sbuf = const_cast<char*>(msg_.data());
#ifdef _MSC_VER
        n = _vsnprintf(sbuf, sze-1, fmt, args);
#else
        n = vsnprintf(sbuf, sze-1, fmt, args);
#endif
        // Negative n is not trapped. We just print anyway.
        va_end(args); 
        sbuf[sze-1] = '\0';
    }
    save();
}
//==================================================================================================================================
Electrode_Error::~Electrode_Error() throw()
{
}
//==================================================================================================================================
Electrode_Error::Electrode_Error() :
    ZZCantera::ZuzaxError()
{
}
//==================================================================================================================================
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
double l0norm(const std::vector<double>& v1, const std::vector<double>& v2, const std::vector<double>& atolVec, const double rtol)
{
    double max0 = 0.0, ee;
    for (size_t k = 0; k < v1.size(); k++) {
        ee = fabs(v1[k] - v2[k]) / MAX(rtol * MAX(fabs(v1[k]), fabs(v2[k])), atolVec[k]);
        if (ee > max0) {
            max0 = ee;
        }
    }
    return max0;
}
//==================================================================================================================================
double relv(double a, double b, double atol)
{
    if (a == 0.0 && b == 0.0) {
        return 0.0;
    }
    double denom = MAX(fabs(a), fabs(b));
    if (denom < atol) {
        denom = MAX(atol, 1.0E-300);
    }
    return fabs((a - b)/ denom);
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
