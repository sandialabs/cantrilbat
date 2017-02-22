/**
 * @file mdp_clockID.cpp
 *    Definitions for a simple class that implements an Ansi C wall clock and interval timer
 *     (see \ref mdpUtil::clockID).
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "mdp_clockID.h"
#include <climits>
#include <stdexcept>
#include <cstdio>

// These need unistd.h, which is a POSIX standard only
//#define HAVE_POSIX
#ifdef  HAVE_UNISTD_H
#include <unistd.h>
#include <errno.h>
#endif


// turn on debugging
#ifndef DEBUG_MODE
#define DEBUG_MODE
#endif
//----------------------------------------------------------------------------------------------------------------------------------
namespace mdpUtil
{

#ifdef MDP_MAKE_CLOCK_GETTIME
#include <mach/mach_time.h>
#define ORWL_NANO (1.0E-9)
#define ORWL_GIGA (1000000000)
static double orwl_timebase = 0.0;
static uint64_t orwl_timestart = 0;
#endif



//  This is only true if we are on a posix system
#ifdef _POSIX_VERSION
//==================================================================================================================================
clockID::clockID(clockid_t cType) :
    timer(),
    clockType_(cType)
{
#ifdef  _POSIX_TIMERS 
#ifndef _POSIX_MONOTONIC_CLOCK
    if (cType == CLOCK_MONOTONIC) {
	 throw std::logic_error("Clock type, CLOCK_MONOTONIC, not supported. Use CLOCK_REALTIME\n");
    }
#endif
#ifndef _POSIX_CPUTIME
    if (cType == CLOCK_PROCESS_CPUTIME_ID) {
	 throw std::logic_error("Clock type, CLOCK_PROCESS_CPUTIME_ID, not supported. Use CLOCK_REALTIME\n");
    }
#endif
#ifndef _POSIX_THREAD_CPUTIME 
    if (cType == CLOCK_THREAD_CPUTIME_ID ) {
	 throw std::logic_error("Clock type, CLOCK_THREAD_CPUTIME_ID, not supported. Use CLOCK_REALTIME\n");
    }
#endif
#endif
#ifdef DEBUG_MODE
    int err = clock_gettime(clockType_, &startLastTicks_);
    if (err) {
#ifdef HAVE_UNISTD_H
	if (errno == EINVAL) {
	    throw std::logic_error("Clock type not supported\n");
	}
#else
	throw std::logic_error("Clock type not supported\n");
#endif
    }
#endif
}
//===============================================================================================================================
clockID::clockID(const clockID& right) :
    timer(right),
    currNumTicks_(right.currNumTicks_),
    startLastTicks_ (right.startLastTicks_)
{
}
//===============================================================================================================================
clockID& clockID::operator=(const clockID& right)
{
    if (this == &right) return *this;
    timer::operator=(right);
    startLastTicks_ = right.startLastTicks_;
    return *this;
}
//===============================================================================================================================
clockID::~clockID()
{
}
//===============================================================================================================================
void clockID::startTime()
{  
#ifdef DEBUG_MODE
    if (running_) {
	throw std::logic_error("Clock is already running");
    }
#endif
    running_ = true;
    clock_gettime(clockType_, &startLastTicks_);
}
//==================================================================================================================================
double clockID::stopTime()
{  
    clock_gettime(clockType_, &currNumTicks_);
    if (!running_) {
	throw std::logic_error("Clock wasn't already running");
    }
    running_ = false;
    double secsLast = (currNumTicks_.tv_sec - startLastTicks_.tv_sec) + 1.0E-9 * (currNumTicks_.tv_nsec - startLastTicks_.tv_nsec);
    storredSeconds_ += secsLast;
    startLastTicks_ = currNumTicks_;
    return secsLast;
}
//==================================================================================================================================
void clockID::clear()
{
    storredSeconds_ = 0.0;
}
//==================================================================================================================================
#ifdef ZZ_MAKE_CLOCK_GETTIME
int  clockID::clock_gettime(clockid_t clk_id, struct timespec* const sp)
{
    // currently this is for the Darwin system only

    if (!orwl_timestart) {
        mach_timebase_info_data_t tb = { 0};
        mach_timebase_info(&tb);
        orwl_timebase = tb.numer;
        orwl_timebase /= tb.denom;
        orwl_timestart = mach_absolute_time();
    }
    double diff = (mach_absolute_time() - orwl_timestart ) * orwl_timebase;
    sp->tv_sec  = diff * ORWL_NANO;
    sp->tv_nsec = diff - (sp->tv_sec * ORWL_GIGA);

    return 0;
}
#endif
//==================================================================================================================================
#endif
}
