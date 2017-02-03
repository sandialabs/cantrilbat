/**
 * @file m1d_exception.cpp
 *
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "m1d_exception.h"
#include "m1d_app.h"
#include "m1d_globals.h"

#include <cmath>
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
m1d_Error::m1d_Error(const std::string& proc, const std::string& msg)
{
    m1d::app()->addError(proc, msg);
}
//==================================================================================================================================
m1d_Error::m1d_Error(const std::string& proc, const char* fmt, ...)
{
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
    m1d::app()->addError(proc, msg_);
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
}
//----------------------------------------------------------------------------------------------------------------------------------
