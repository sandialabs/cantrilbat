/**
 * @file mdp_stringUtils.h
 *   Declarations for various string utilities.
 */
/*
 * Copywrite 2015 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef MDP_STRINGUTILS_H
#define MDP_STRINGUTILS_H

#include <string>

namespace mdpUtil 
{

    //!  Compares two strings to see if they are equal ignoring the difference in case
    /*!
     *   @param[in]  a  first string
     *   @param[in]  b  second string
     *
     *   @return Returns true if they are equal, false if they are unequal
     */
    bool LowerCaseStringEquals(const std::string& a, const std::string& b);

    //! Convert a double into a c++ string
    /*!
     *  @param[in]  x     double to be converted
     *  @param[in]  fmt   Format to be used (printf style)
     *
     *  @return     Returns the string
     */
    std::string fp2str(const double x, const char* const fmt = "%g");

    //!  Convert an int to a string using a format converter
    /*!
     *  @param[in]  n     int to be converted
     *  @param[in]  fmt   format converter for an int int the printf command
     *
     *  @return     Returns the string
     */
    std::string int2str(const int n, const char* const fmt ="%d");

}

#endif
