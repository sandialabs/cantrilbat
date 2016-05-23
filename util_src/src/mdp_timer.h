/**
 * @file mdp_timer.h
 *    Declarations for a simple class that implements a clock timer
 *   (see \ref mdpUtil::timer).
 */
/*
 * Copyright 2014 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "config.h"

#ifndef MDP_TIMER_H
#define MDP_TIMER_H

#ifdef HAVE_UNISTD_H
#include "unistd.h"
#endif

namespace mdpUtil
{
//===================================================================================================================================
//! The class provides the wall clock timer in seconds for POSIX compliant systems
/*!
 *  This routine relies on the POSIX clock_gettime() routine, which
 *  includes four types of clocks.  There are several types of clocks. See man mage
 *  for clock_gettime() for more information.
 *
 *     CLOCK_REALTIME
 *     CLOCK_MONOTONIC
 *     CLOCK_MONOTONIC_RAW 
 *     CLOCK_PROCESS_CPUTIME_ID
 *     CLOCK_THREAD_CPUTIME_ID
 *
 *  An example of how to use the timer is given below. timeToDoCalcs
 *  contains the processor clock time calculated for the operation.
 *
 *  This clock may be started and stopped repeatedly. The time for each incremental
 *  interval will be reported.
 * 
 *  If this routine is used, the linker library -lrt must be added into the linker line.
 *
 *  @code
 *    clockID wc(CLOCK_PROCESS_CPUTIME_ID);
 *    wc.startTime();
 *    do_calculations_();
 *    wc.stopTime();
 *    double timeToDoCalcs = wc.reportTime();
 *  @endcode
 *
 *  In general, the process to be timed can be measured in nanoseconds with
 *  this implementation of the clock.
 *
 * @ingroup globalUtilFuncs
 *
 */
class timer
{
public:

    //! Constructor
    /*!
     *
     */
    timer() :
        running_(false),
        storredSeconds_(0.0)
    {
    }

    //! Copy constructor
    /*!
     *   Creates a clock from another clock
     *
     * @param right  Item to be copied.
     */
    timer(const timer& right) :
    running_(right.running_),
    storredSeconds_(right.storredSeconds_)
    {
    }

    //! Assign the timer info from one timer to another
    /*!
     *  @param right object to be copied
     *
     *  @return returns a reference to the current object
     */
    timer& operator=(const timer& right)
    {
	if (&right == this) return *this;
        running_ = right.running_;
	storredSeconds_ = right.storredSeconds_;	
        return *this;
    }

    virtual ~timer()
    {
    }

    //! Starts the timer on a new interval to be timed.
    virtual void startTime() = 0;

    //! Stop the timer
    /*!
     *      The timer is stopped. The incremental time is added into the total time for the clock
     *
     *  @return Returns the incremental time from the last increment
     */
    virtual double stopTime() = 0;

    //! Reports the time 
    /*!
     *  @return Returns the time (secs) spent in all intervals timed by this clock
     */
    virtual double reportTime() const
    {
        return storredSeconds_;
    } 

    //! Clears all of the timing information 
    virtual void clear()
    {
        storredSeconds_ = 0;
    }

protected:
    //! This is a flag stating whether interval timing is turned on
    bool running_;

    //! Storred seconds for calculating the timing of the operations
    double storredSeconds_;
 
};
//===================================================================================================================================


}
#endif
