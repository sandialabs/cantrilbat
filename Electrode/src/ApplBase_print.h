/**
 *  @file ApplBase_print.h
 *
 */
/*
 * $Id: ApplBase_print.h 571 2013-03-26 16:44:21Z hkmoffa $
 *
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */


#ifndef APPLBASE_PRINT_H
#define APPLBASE_PRINT_H

#include <map>
#include <string>

namespace ca_ab
{
//! Print a single character multiple times
/*!
 *  @param[in]                c                   character to be printed
 *  @param[in]                nTimes              multiple of times
 */
extern void print_char(const char c, const int nTimes);

//! Print a boolean using "yes" or "no"
/*!
 *  @param[in]                b                   value of bool
 */
extern void print_bool(const bool b);

//! Print an int using a width
/*!
 *  @param[in]                i                   value of int
 *  @param[in]                w                   value of width
 */
extern void pr_if(const int i, const int w);

//! Print a string using a certain amount of width
/*!
 *  Padding is placed on the left of the string.
 *
 *  @param[in]                s                   value of the string
 *  @param[in]                w                   value of the width
 */
extern void pr_sf(const std::string s, const int w);

//! Print a string in a fixed space, using left justification, and cropping
/*!
 *  @param[in]                s                   String
 *  @param[in]                w                   width of the space
 *  @param[in]                crop                boolean indicating you want to crop. Defaults to false.
 */
extern void pr_sf_lj(const std::string s, const int w, const int crop = 0);

//!  print with a fixed precision and width in fixed format.
/*!
 *  @param[in]                d                   value of the double to be printed
 *  @param[in]                w                   width
 *  @param[in]                p                   precision
 */
extern void pr_df(const double d, const int w, const int p);

//! Print with a fixed precision and a variable width in a
//! fixed format (i.e., non scientific notation)
/*!
 *  @param[in]                d                   value of the double to be printed
 *  @param[in]                p                   precision
 */
extern void pr_dfp(const double d, const int p);

//!  print with a fixed width and prescision in scientific format
/*!
 *  @param[in]                d                   value of the double to be printed
 *  @param[in]                w                   width
 *  @param[in]                p                   precision
 */
extern void pr_de(const double d, const int w, const int p);

//! Print a double using a best representation for the width and precision
/*!
 *  @param[in]                d                   value of the double to be printed
 *  @param[in]                w                   width
 *  @param[in]                p                   precision
 */
extern void pr_dg(const double d, const int w, const int p);

//! Indent according to sectional level
/*!
 *  @param[in]                i                   level of sectional indentation
 */
extern void dnt(const int i);

//! Print the value of a map<string,double> to standard out
/*!
 *  The map is printed all on one line.
 *
 *  @param[in]                m                   Map to be printed
 *  @param[in]                prefix              string to be printed before the values
 */
extern void print_map(const std::map<std::string,double>& m,  const std::string& prefix);

}
#endif
