/**
 * @file m1d_globals.cpp
 * This file contains globally defined functions and data
 * that are not part of any class.
 */

/*
 *  $Id: m1d_globals.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "m1d_globals.h"
#include "m1d_app.h"

#include <cstdlib>
#include <cmath>

namespace m1d
{

int s_printLvl_TimingInformation = 1;

int s_printLvl_DebugTables = 1;
//======================================================================================================================
// By default predictor_corrector printing is turned on, at least to the printLvl_ level.
//======================================================================================================================

void readEnvironmentalVariables() {

     char *PC_PRINTING = getenv("M1D_PRINT_TIMINGINFORMATION");
     if (PC_PRINTING) {
        if (PC_PRINTING[0] == 'f' || PC_PRINTING[0] == 'F' || PC_PRINTING[0] == '0') {
           m1d::s_printLvl_TimingInformation = 0;
        } else  if (PC_PRINTING[0] == 't' || PC_PRINTING[0] == 'T' || PC_PRINTING[0] == '1') {
           m1d::s_printLvl_TimingInformation = 1;
        }
     }
     PC_PRINTING = getenv("M1D_PRINT_DEBUGTABLES");
     if (PC_PRINTING) {
        if (PC_PRINTING[0] == 'f' || PC_PRINTING[0] == 'F' || PC_PRINTING[0] == '0') {
           m1d::s_printLvl_DebugTables = 0;
        } else  if (PC_PRINTING[0] == 't' || PC_PRINTING[0] == 'T' || PC_PRINTING[0] == '1') {
           m1d::s_printLvl_DebugTables = 1;
        }
     }

}


//==============================================================================
// Return a pointer to the application object
Appl*
app()
{
  return Appl::Instance();
}
//==============================================================================
void
setError(std::string r, std::string msg)
{
  app()->addError(r, msg);
}
//==============================================================================
void
showErrors(std::ostream& f)
{
  app()->showErrors(f);
}
//==============================================================================
void
popError()
{
  app()->popError();
}
//==============================================================================
std::string
intToString(const int p)
{
  char buf[64];
  sprintf(buf, "%d", p);
  return std::string(buf);
}
//==============================================================================
 bool checkDblAgree(double d1, double d2, double rtol, double atol)
{
    double dfabs = fabs(d1 - d2);
    double denom = 0.5 * (fabs(d1) + fabs(d2)) * rtol;
    denom = std::max(denom, atol);
    if (dfabs < denom) return true;
    return false;
}
//==============================================================================

}
