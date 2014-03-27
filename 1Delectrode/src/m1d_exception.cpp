/**
 * @file m1d_exception.cpp
 *
 */

/*
 *  $Id: m1d_exception.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "m1d_exception.h"
#include "m1d_app.h"
#include "m1d_globals.h"

#include <cmath>

namespace m1d
{

m1d_Error::m1d_Error(const std::string &proc, const std::string &msg)
{
  m1d::app()->addError(proc, msg);
}
//======================================================================================================
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
//==============================================================================
}
