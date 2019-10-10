/**
 *  @file sort.h
 *       Declarations for various sort routines.
 *
 */

#ifndef MDP_SORT_H
#define MDP_SORT_H

#include <vector>

namespace mdpUtil
{

//! Given two arrays x and y, sort the (x,y) pairs by the x
//! values. This version is for floating-point x, and integer y.
/*!
 *   @param[in,out]           x   Reference to  Vector of values of x.
 *   @param[in,out]           y   Reference to vector of values of y. These are sort-ordered according to x.
 */
extern void heapsort(std::vector<double>& x, std::vector<int>& y);

//! Given two arrays x and y, sort the (x,y) pairs by the x
//! values. This version is for floating-point x, and floating-point y.
/*!
 *   @param[in,out]           x   Reference to  Vector of values of x.
 *   @param[in,out]           y   Reference to vector of values of y. These are sort-ordered according to x.
 *
 */
extern void heapsort(std::vector<double>& x,  std::vector<double>& y);

template<typename T >
void heapsortT(std::vector<double>& x,  std::vector< T >& y);

}

#endif
