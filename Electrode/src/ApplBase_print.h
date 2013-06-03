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
extern void print_char(const char c, const int nTimes);

extern void print_bool(const bool b);

extern void pr_if(const int i, const int w);

extern void pr_sf(const std::string s, const int w);

//! Print a string in a fixed space, using left justification, and cropping
/*!
 *  @param s       String
 *  @param w       width of the space
 *  @param p       boolean indicating you want to crop
 */
extern void pr_sf_lj(const std::string s, const int w, const int crop = 0);
extern void pr_df(const double d, const int w, const int p);
extern void pr_dfp(const double d, const int p);
extern void pr_de(const double d, const int w, const int p);
extern void pr_dg(const double d, const int w, const int p);
extern void dnt(const int i);

extern void print_map(const std::map<std::string,double>& m,
                      const std::string& prefix);

}
#endif
