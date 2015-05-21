/**
 * @file clockWC.h
 *    Declarations for a simple class that implements an Ansi C wall clock timer
 *   (see \ref mdpUtil::clockWC).
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef MDP_CLOCKWC_H
#define MDP_CLOCKWC_H

#include <ctime>

namespace mdpUtil
{

//! The class provides the wall clock timer in seconds
/*!
 *  This routine relies on the ANSI C routine, clock(), for
 *  its basic operation. Therefore, it should be fairly
 *  portable.
 *
 * The clock will rollover if the calculation is long enough.
 * The wraparound time is roughly 72 minutes for a 32 bit system.
 * This object senses that by seeing if the raw tick counter is
 * has decreased from the last time. If it senses a wraparound has
 * occurred, it increments an internal counter to account for this.
 * Therefore, for long calculations, this object must be called
 * at regular intervals for the seconds timer to be accurate.
 *
 * An example of how to use the timer is given below. timeToDoCalcs
 * contains the wall clock time calculated for the operation.
 *
 *
 *  @code
 *   clockWC wc;
 *   do_hefty_calculations_atLeastgreaterThanAMillisecond();
 *   double timeToDoCalcs = wc.secondsWC();
 *  @endcode
 *
 *  In general, the process to be timed must take more than a millisecond
 *  for this clock to enough of a significant resolution to be
 *  accurate.
 *
 * @ingroup globalUtilFuncs
 *
 */
class clockWC
{
public:
    //! Constructor
    /*!
     * This also serves to initialize the ticks within the object
     */
    clockWC();

    clockWC(const clockWC& right);

    clockWC&  operator=(const clockWC& right);
   
    void startWC();

    //! Resets the internal counters and returns the wall clock time in seconds
    void startTime();

    void restartTime();

    double stopTime();

    double reportTime() const;

    void checkWC();

    //! Returns the wall clock time in seconds since the last reset.
    /*!
     *  Returns system cpu and wall clock time in seconds. This is a strictly
     *  Ansi C timer, since clock() is defined as an Ansi C function. On some
     *  machines clock() returns type unsigned long (HP) and on others (SUN)
     *  it returns type long. An attempt to recover the actual time for clocks
     *  which have rolled over is made also. However, it only works if this
     *  function is called fairly regularily during the solution procedure.
     */
    double secondsWC();

private:

    //! Internal constant containing clock ticks per second
    static const double s_inv_clocks_per_sec;

    //! Internal constant containing the total number of ticks per rollover.
    //static const double s_clockTickWidth_;

    //! Counters the value of the number of ticks from the current call.
    clock_t currNumTicks_;

    //! Counters the value of the number of ticks from the last call.
    clock_t lastNumTicks_;

    //! This is a flag stating whether interval timing is turned on
    bool running_;

    //! Counter containing the value of the number of ticks from the last start call
    clock_t startLastTicks_;

    //! Counter containing the value of the number of ticks from the first call to start();
    clock_t startTicksWC_;

    //! Storred seconds for calculating the timing of the operations
    /*!
     *  This handles clock rollovers
     */
    double storredSeconds_;

    //! Storred seconds for calculating the timing of the wall clock
    /*!
     *  This handles clock rollovers
     */
    double storredSecondsWC_;

    void clear();
};
}
#endif
