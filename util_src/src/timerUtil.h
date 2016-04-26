/**
 * @file timer.h
 *    Declarations for a simple class that implements a clock timer
 *   (see \ref mdpUtil::clockID).
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
