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

    //!  Compares two strings to see if they are equal ignoring the case
    /*!
     *   @param[in]  a  first string
     *   @param[in]  b  second string
     *
     *   @return Returns true if they are equal
     */
    bool LowerCaseStringEquals(const std::string& a, const std::string& b);

}

#endif
