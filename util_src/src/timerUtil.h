/**
 * @file timerUtil.h
 *    Declarations for a factory funcitons that return a clock timer based on the operating system
 *   (see \ref mdpUtil::timerUtil).
 */
/*
 * Copyright 2014 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef MDP_TIMERUTIL_H
#define MDP_TIMERUTIl_H

#include "timer.h"
#include "wallClock.h"
#include "clockID.h"

/*
 *  We take advantage of the _POSIX_VERSION define to tell whether the system is posix compliant.
 */

namespace mdpUtil
{
//===================================================================================================================================

timer* newTimer()
{
#ifdef _POSIX_VERSION
   return new clockID();
#else
   return new wallClock();
#endif
};

//===================================================================================================================================
timer* newCPUTimer()
{
#ifdef _POSIX_VERSION
   return new clockID(CLOCK_PROCESS_CPUTIME_ID);
#else
   return new wallClock();
#endif
};
//===================================================================================================================================
timer* newWallTimer()
{
#ifdef _POSIX_VERSION
   return new clockID(CLOCK_REALTIME);
#else
   return new wallClock();
#endif
};
//===================================================================================================================================
}
#endif
