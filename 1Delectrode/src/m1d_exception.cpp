/**
 * @file m1d_exception.cpp
 *
 */

/*
 *  $Id: m1d_exception.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "m1d_exception.h"
#include "m1d_app.h"
#include "m1d_globals.h"

namespace m1d
{

m1d_Error::m1d_Error(const std::string &proc, const std::string &msg)
{
  m1d::app()->addError(proc, msg);
}
//==============================================================================
}
