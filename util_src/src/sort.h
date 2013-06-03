/**
 *  @file sort.h
 *
 *  $Id: sort.h 485 2010-06-10 19:55:25Z hkmoffa $
 */

#ifndef MDP_SORT_H
#define MDP_SORT_H

#include <vector>

namespace mdpUtil {

    /// Given two arrays x and y, sort the (x,y) pairs by the x
    /// values. This version is for floating-point x, and integer y.
    void heapsort(std::vector<double> & x, std::vector<int> & y);

    /// Given two arrays x and y, sort the (x,y) pairs by the x
    /// values. This version is for floating-point x, and
    /// floating-point y.
    void heapsort(std::vector<double> &x,  std::vector<double> &y);
}

#endif
