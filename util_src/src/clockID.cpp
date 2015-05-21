/**
 * @file clockID.cpp
 *    Definitions for a simple class that implements an Ansi C wall clock and interval timer
 *     (see \ref mdpUtil::clockID).
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "clockID.h"
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

namespace mdpUtil
{

//==================================================================================================================================
clockID::clockID(clockid_t cType) :
    clockType_(cType),
    running_(false),
    storredSeconds_(0.0)
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
    currNumTicks_(right.currNumTicks_)
{
    operator=(right);
}
//===============================================================================================================================
clockID& clockID::operator=(const clockID& right)
{
    if (this == &right) return *this;
    running_      = right.running_;
    startLastTicks_ = right.startLastTicks_;
    storredSeconds_ = right.storredSeconds_;
    return *this;
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
double clockID::reportTime() const
{
    return storredSeconds_;
}
//==================================================================================================================================
void clockID::clear()
{
    storredSeconds_ = 0.0;
}
//==================================================================================================================================
}
