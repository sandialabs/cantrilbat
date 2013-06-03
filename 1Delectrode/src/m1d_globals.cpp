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

namespace m1d
{

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
}
