/*
 * FDJacobianTests.cpp
 *
 */
#include "../../../config.h"
#include "gtest/gtest.h"
#include "clockWC.h"
#include "clockID.h"

using namespace mdpUtil;
using namespace testing;

clockWC ptick;

double timeSpent(int n)
{
    clockWC tick;

    tick.startTime();
    //  Note, this is system time -> doesn't count
    sleep(n);
    double ff = 1.0;
    for (int i = 2; i < 88888880; i++) {
	ff *= (i - 1.0) / (i - 0.98);
    }
    FILE* ofp = fopen("testerOut.txt", "w");
    fprintf(ofp, "ff = %g\n", ff);
    fclose(ofp);
    tick.stopTime();
   
    return tick.reportTime();
}

double timeSpentID(int n, clockid_t ctype)
{
    clockID tick(ctype);
    tick.startTime();
    //  Note, this is system time -> doesn't count
    sleep(n);
    double ff = 1.0;
    for (int i = 2; i < 88888880; i++) {
	ff *= (i - 1.0) / (i - 0.98);
    }
    FILE* ofp = fopen("testerOut.txt", "w");
    fprintf(ofp, "ff = %g\n", ff);
    fclose(ofp);
    tick.stopTime();
    return tick.reportTime();
}

double do_work(int n)
{
    double ff = 1.0;
    for (int i = 2; i < 88888880; i++)
    {
	ff *= (i - 1.0) / (i - 0.98);
    }
    return ff;
}

void time_work(int m)
{



}



TEST(basicFunction, timeSpent)
{
    EXPECT_NEAR(1.2, timeSpent(5) , 1.);
}

//
// the clockID function is only defined on posix systems
//
#ifdef HAVE_UNISTD_H
TEST(basicFunction, timeSpentID_cpu)
{
    EXPECT_NEAR(1.2, timeSpentID(5, CLOCK_PROCESS_CPUTIME_ID) , 1.);
}
TEST(basicFunction, timeSpentID_real)
{
    EXPECT_NEAR(6.2, timeSpentID(5, CLOCK_REALTIME) , 1.);
}
#endif

int
main(int argc, char** argv)
{
    printf("Running main() from timerTests.cpp \n");
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
