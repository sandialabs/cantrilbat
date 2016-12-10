/**
 * @file mdp_wallClock.cpp
 *    Definitions for a simple class that implements an Ansi C wall clock and interval timer
 *     (see \ref mdpUtil::wallClock).
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "mdp_wallClock.h"
#include <climits>
#include <stdexcept>
#ifndef DEBUG_MODE
#define DEBUG_MODE
#endif
namespace mdpUtil
{
//===============================================================================================================================
const double wallClock::s_inv_clocks_per_sec(1.0/(double)CLOCKS_PER_SEC);
//
//  This isn't correct. Really should use LONG_MAX
//
//   const double wallClock::s_clockTickWidth_((double)(1L<<((int)sizeof(clock_t)*8-2))*4./(double)CLOCKS_PER_SEC);
//===============================================================================================================================
wallClock::wallClock() :
    timer(),
    lastNumTicks_(clock()),
    startLastTicks_(lastNumTicks_),
    startTicksWC_(lastNumTicks_),
    storredSecondsWC_(0.0)
{
#ifdef DEBUG_MODE
    if (sizeof(clock_t) != sizeof(long int)) {
	throw std::logic_error("wrong assumptions about clock_t");
    }
    clock_t tt= - 45;
    tt -= 55;
    if (tt != -100) {
	throw std::logic_error("wrong assumptions about clock_t");
    }
#endif
}
//===============================================================================================================================
wallClock::wallClock(const wallClock& right) :
    timer(),
    currNumTicks_(right.currNumTicks_)
{
    operator=(right);
}
//===============================================================================================================================
wallClock& wallClock::operator=(const wallClock& right)
{
    if (this == &right) return *this;
    timer::operator=(right);
    lastNumTicks_ = right.lastNumTicks_;
    startLastTicks_ = right.startLastTicks_;
    startTicksWC_ = right.startTicksWC_;
    storredSecondsWC_ = right.storredSecondsWC_;
    return *this;
}
//===============================================================================================================================
void wallClock::startWC()
{
    lastNumTicks_ = clock();
    startTicksWC_ = startLastTicks_ = lastNumTicks_;
    storredSecondsWC_ = 0.0;
    running_ = true;
}
//===============================================================================================================================
void wallClock::startTime()
{  
 #ifdef DEBUG_MODE
    if (running_) {
	throw std::logic_error("Clock is already running");
    }
#endif
    running_ = true;
    startLastTicks_ = lastNumTicks_ = clock();
}
//===============================================================================================================================
double wallClock::stopTime()
{
    currNumTicks_ = clock();
    if (!running_) {
	throw std::logic_error("Clock wasn't already running");
    }
    running_ = false;
    double secsLast = 0.0;
    if (currNumTicks_ < lastNumTicks_) {
	secsLast = ((LONG_MAX - startLastTicks_) * s_inv_clocks_per_sec) + s_inv_clocks_per_sec;
	storredSeconds_ += secsLast;  
	startLastTicks_ = (clock_t) 0L;
	storredSecondsWC_ += ((LONG_MAX - startTicksWC_) * s_inv_clocks_per_sec) + s_inv_clocks_per_sec;
	startTicksWC_ = (clock_t) 0L;
    } 
    secsLast += (currNumTicks_ - startLastTicks_) * s_inv_clocks_per_sec;
    storredSeconds_ += secsLast;
    startLastTicks_ = lastNumTicks_ = currNumTicks_;
    return secsLast;
}
//===============================================================================================================================
double wallClock::reportTime() const
{
    return storredSeconds_;
}
//===============================================================================================================================
void wallClock::checkWC()
{   
    currNumTicks_ = clock();
    if (currNumTicks_ < lastNumTicks_) {
	storredSecondsWC_ += ((LONG_MAX - startTicksWC_) * s_inv_clocks_per_sec) + s_inv_clocks_per_sec;
	startTicksWC_ = (clock_t) 0L;
    }
    lastNumTicks_ = currNumTicks_;
}
//===============================================================================================================================
double wallClock::secondsWC()
{
    currNumTicks_ = clock();
    if (currNumTicks_ < lastNumTicks_) {
	storredSecondsWC_ += ((LONG_MAX - startTicksWC_) * s_inv_clocks_per_sec) + s_inv_clocks_per_sec;
	startTicksWC_ = (clock_t) 0L;
    }
    lastNumTicks_ = currNumTicks_;
    return storredSecondsWC_ + (currNumTicks_ - startTicksWC_) * s_inv_clocks_per_sec;
}
//===============================================================================================================================
void wallClock::clear()
{
    storredSeconds_ = 0.0;
    storredSecondsWC_ = 0.0;
}
//===============================================================================================================================
}
