/**
 * @file mdp_stringUtils.cpp
 *   Declarations for various string utilities.
 */
/*
 * Copywrite 2015 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "mdp_stringUtils.h"

#include <cctype>
#include <cstdio>

namespace mdpUtil {
//============================================================================================================
bool LowerCaseStringEquals(const std::string& a, const std::string& b)
{
     size_t sz = a.size();
     if (b.size() != sz) {
         return false;
     }
     for (size_t i = 0; i < sz; ++i) {
         if (tolower(a[i]) != tolower(b[i])) {
             return false;
         }
     }
     return true;
}
//============================================================================================================
std::string fp2str(const double x, const char* const fmt)
{
    char buf[64];
    int n = snprintf(buf, 63, fmt, x);
    if (n > 0) {
        buf[63] = '\0';
        return std::string(buf);
    }
    return std::string(" ");
}
//============================================================================================================
std::string int2str(const int n, const char* const fmt)
{
    char buf[32];
    int m = snprintf(buf, 30, fmt, n);
    if (m > 0) {
        buf[29] = '\0';
        return std::string(buf);
    }
    return std::string(" ");
}
//============================================================================================================

}
