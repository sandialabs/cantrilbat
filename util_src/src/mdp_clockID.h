/**
 * @file mdp_clockID.h
 *    Declarations for a simple class that implements a POSIX compliant  clock timer
 *   (see \ref mdpUtil::clockID).
 */
/*
 * Copyright 2014 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "mdp_timer.h"

#ifndef MDP_CLOCKID_H
#define MDP_CLOCKID_H

#include <ctime>

#ifdef  DOXYGEN_SHOULD_IGNORE_THIS
#ifndef _POSIX_VERSION
//! fake version number to get doxygen to compile this
#define  _POSIX_VERSION 12
#endif
#endif

#ifdef _POSIX_VERSION

namespace mdpUtil
{

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
class clockID : public timer
{
public:

    //! Constructor
    /*!
     *  @param cType   The type of clock to be used. We default to CLOCK_PROCESS_CPUTIME_ID
     *                 which measure the CPUTIME of the current process. This does not
     *                 measure the real time.
     */
    clockID(clockid_t cType = CLOCK_PROCESS_CPUTIME_ID);

    //! Copy constructor
    /*!
     *   Creates a clock from another clock
     *
     * @param right  Item to be copied.
     */
    clockID(const clockID& right);

    //! Assign the timer info from one timer to another
    /*!
     *  @param right object to be copied
     *
     *  @return returns a reference to the current object
     */
    clockID& operator=(const clockID& right);

    virtual ~clockID();
   
    //! Starts the timer on a new interval to be timed.
    virtual void startTime();

    //! Stop the timer
    /*!
     *      The timer is stopped. The incremental time is added into the total time for the clock
     *
     *  @return Returns the incremental time from the last increment
     */
    virtual double stopTime();

    //! Clears all of the timing information 
    virtual void clear();
  
private:

    //! Clock type, See the man page for clock_gettime() for possible values
    clockid_t clockType_;

    //! Counters the value of the number of ticks from the current call.
    struct timespec currNumTicks_;

    //! Counter containing the value of the number of ticks from the last start call
    struct timespec startLastTicks_;
};
}
#endif
#endif
